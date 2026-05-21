# =============================================================================
# RNA Motif Discovery Explorer
# Script: motif_detection.R
# Purpose: k-mer extraction, motif ranking, enrichment scoring, significance
# =============================================================================

library(Biostrings)
library(dplyr)
library(tibble)
library(stringr)
library(universalmotif)

# -----------------------------------------------------------------------------
# extract_kmers
# Extract all k-mers of a given length from a single sequence string.
# Returns a character vector of k-mers.
# -----------------------------------------------------------------------------
extract_kmers <- function(seq, k) {
  n <- nchar(seq)
  if (n < k) return(character(0))
  sapply(seq_len(n - k + 1), function(i) substr(seq, i, i + k - 1))
}

normalize_rna_text <- function(seq) {
  gsub("T", "U", toupper(seq))
}

find_overlapping_positions <- function(seq, pattern) {
  seq <- normalize_rna_text(seq)
  pattern <- normalize_rna_text(pattern)
  k <- nchar(pattern)
  n <- nchar(seq)
  if (n < k) return(integer(0))
  starts <- seq_len(n - k + 1)
  candidates <- vapply(
    starts,
    function(i) substr(seq, i, i + k - 1),
    character(1)
  )
  unname(starts[candidates == pattern])
}

# -----------------------------------------------------------------------------
# count_kmers_dataset
# For each k-mer in the full dataset, count:
#   - total occurrences across all sequences
#   - number of sequences containing it (sequence frequency)
# Returns a tibble sorted by total_count descending.
# -----------------------------------------------------------------------------
count_kmers_dataset <- function(rna_stringset, k) {
  seqs_char <- as.character(rna_stringset)
  n_seqs    <- length(seqs_char)

  # Collect all k-mers with their source sequence index
  all_kmers <- lapply(seq_along(seqs_char), function(i) {
    km <- extract_kmers(seqs_char[i], k)
    if (length(km) == 0) return(NULL)
    data.frame(kmer = km, seq_idx = i, stringsAsFactors = FALSE)
  })
  all_kmers <- bind_rows(all_kmers)

  if (nrow(all_kmers) == 0) {
    return(tibble(kmer = character(0), total_count = integer(0),
                  seq_count = integer(0), seq_frequency = numeric(0)))
  }

  # Remove k-mers with ambiguous nucleotides (N)
  all_kmers <- all_kmers[!grepl("N", all_kmers$kmer), ]

  if (nrow(all_kmers) == 0) {
    return(tibble(kmer = character(0), total_count = integer(0),
                  seq_count = integer(0), seq_frequency = numeric(0)))
  }

  summary <- all_kmers %>%
    group_by(kmer) %>%
    summarise(
      total_count   = n(),
      seq_count     = n_distinct(seq_idx),
      .groups       = "drop"
    ) %>%
    mutate(
      seq_frequency = seq_count / n_seqs
    ) %>%
    arrange(desc(total_count))

  summary
}

# -----------------------------------------------------------------------------
# calculate_enrichment_scores
# Computes enrichment relative to background (uniform RNA nucleotide model).
# Expected frequency = (0.25)^k for a random k-mer under uniform assumption.
# Enrichment score = observed_freq / expected_freq.
# -----------------------------------------------------------------------------
calculate_enrichment_scores <- function(kmer_counts, total_positions) {
  if (nrow(kmer_counts) == 0) {
    return(kmer_counts %>%
             mutate(observed_freq = numeric(0), expected_freq = numeric(0),
                    expected_count = numeric(0), enrichment_score = numeric(0)))
  }

  k <- nchar(kmer_counts$kmer[1])

  # Background: uniform model p(A)=p(U)=p(G)=p(C)=0.25
  expected_freq  <- 0.25^k
  expected_count <- expected_freq * total_positions

  kmer_counts %>%
    mutate(
      observed_freq  = total_count / total_positions,
      expected_freq  = expected_freq,
      expected_count = expected_count,
      enrichment_score = log2((observed_freq + 1e-10) / (expected_freq + 1e-10))
    )
}

# -----------------------------------------------------------------------------
# calculate_pvalues
# Binomial test for each k-mer: is observed count significantly higher
# than expected under the background model?
# Uses Bonferroni correction for multiple testing.
# -----------------------------------------------------------------------------
calculate_pvalues <- function(kmer_enriched, total_positions) {
  if (nrow(kmer_enriched) == 0) {
    return(kmer_enriched %>%
             mutate(p_value = numeric(0), p_adjusted = numeric(0),
                    significant = logical(0), neg_log10_p = numeric(0)))
  }

  k             <- nchar(kmer_enriched$kmer[1])
  expected_prob <- 0.25^k
  n_tests       <- nrow(kmer_enriched)

  pvals <- mapply(function(obs, total) {
    tryCatch(
      binom.test(obs, total, p = expected_prob, alternative = "greater")$p.value,
      error = function(e) 1.0
    )
  }, kmer_enriched$total_count, rep(total_positions, nrow(kmer_enriched)))

  kmer_enriched %>%
    mutate(
      p_value       = pvals,
      p_adjusted    = p.adjust(pvals, method = "bonferroni"),
      significant   = p_adjusted < 0.05,
      neg_log10_p   = -log10(p_value + 1e-300)
    )
}

# -----------------------------------------------------------------------------
# compute_positional_distribution
# For each of the top N motifs, compute how their occurrences are distributed
# across relative positions (0–1 normalized) in the sequence collection.
# Returns a tibble with kmer, position (absolute), and relative_position.
# -----------------------------------------------------------------------------
compute_positional_distribution <- function(rna_stringset, top_kmers,
                                            n_bins = 20) {
  seqs_char <- normalize_rna_text(as.character(rna_stringset))
  top_kmers <- unique(normalize_rna_text(as.character(top_kmers)))
  top_kmers <- top_kmers[!is.na(top_kmers) & nzchar(top_kmers)]

  if (length(seqs_char) == 0 || length(top_kmers) == 0) {
    return(tibble(kmer = character(0), abs_position = integer(0),
                  seq_length = integer(0), relative_position = numeric(0),
                  position_bin = integer(0)))
  }

  pos_data <- lapply(top_kmers, function(km) {
    k <- nchar(km)
    pos_list <- lapply(seq_along(seqs_char), function(i) {
      seq_len_i <- unname(nchar(seqs_char[i]))
      if (seq_len_i < k) return(NULL)
      hits <- find_overlapping_positions(seqs_char[i], km)
      if (length(hits) == 0) return(NULL)
      data.frame(
        kmer              = unname(km),
        abs_position      = as.integer(hits),
        seq_length        = seq_len_i,
        relative_position = as.integer(hits) / seq_len_i,
        stringsAsFactors  = FALSE
      )
    })
    bind_rows(pos_list)
  })

  pos_df <- bind_rows(pos_data)

  # Bin relative positions
  if (nrow(pos_df) > 0) {
    pos_df <- pos_df %>%
      mutate(
        position_bin = cut(relative_position,
                           breaks = seq(0, 1, length.out = n_bins + 1),
                           include.lowest = TRUE, labels = FALSE)
      )
  }

  pos_df
}

# -----------------------------------------------------------------------------
# build_pwm_from_kmers
# Given a set of k-mer strings (same length), build a position weight matrix
# compatible with universalmotif.
# -----------------------------------------------------------------------------
build_pwm_from_kmers <- function(kmer_vec) {
  k    <- nchar(kmer_vec[1])
  nts  <- c("A", "U", "G", "C")

  # Count matrix (4 x k)
  counts <- matrix(0, nrow = 4, ncol = k, dimnames = list(nts, NULL))
  for (km in kmer_vec) {
    chars <- strsplit(km, "")[[1]]
    for (pos in seq_len(k)) {
      nt <- chars[pos]
      if (nt %in% nts) counts[nt, pos] <- counts[nt, pos] + 1
    }
  }

  # Pseudo-count normalization
  freq_matrix <- (counts + 0.1) / (colSums(counts) + 0.4)

  tryCatch({
    motif <- universalmotif::create_motif(freq_matrix,
                                          alphabet = "RNA",
                                          type = "PPM",
                                          name = paste(kmer_vec[1], "...",
                                                        sep = ""))
    motif
  }, error = function(e) NULL)
}

# -----------------------------------------------------------------------------
# rank_motifs
# Combined ranking: frequency rank + enrichment rank + significance rank.
# Lower rank_score = better motif.
# -----------------------------------------------------------------------------
rank_motifs <- function(motif_table) {
  if (nrow(motif_table) == 0) return(motif_table)

  motif_table %>%
    mutate(
      rank_freq        = rank(-total_count,    ties.method = "min"),
      rank_enrichment  = rank(-enrichment_score, ties.method = "min"),
      rank_significance = rank(p_adjusted,     ties.method = "min"),
      combined_rank    = rank_freq + rank_enrichment + rank_significance
    ) %>%
    arrange(combined_rank)
}

# -----------------------------------------------------------------------------
# motif_discovery_pipeline
# Full motif discovery: k-mer counts -> enrichment -> p-values -> ranking
# -> positional distribution -> PWMs for top motifs
# Parameters:
#   rna_stringset  : preprocessed RNAStringSet
#   k              : motif length (3–15)
#   min_freq       : minimum sequence frequency threshold (0–1)
#   min_count      : minimum total occurrence count
#   top_n          : number of top motifs to analyse in depth
# -----------------------------------------------------------------------------
motif_discovery_pipeline <- function(rna_stringset,
                                      k           = 6,
                                      min_freq    = 0.1,
                                      min_count   = 3,
                                      top_n       = 20,
                                      n_pos_bins  = 20,
                                      verbose     = TRUE) {
  log <- character(0)
  log_msg <- function(msg) {
    if (verbose) message(msg)
    log <<- c(log, msg)
  }

  log_msg(sprintf("=== Motif Discovery Pipeline (k=%d) ===", k))

  if (length(rna_stringset) == 0) {
    stop("No sequences supplied for motif discovery.")
  }

  seqs_char       <- as.character(rna_stringset)
  total_positions <- sum(pmax(nchar(seqs_char) - k + 1, 0))
  total_positions <- max(total_positions, 1)

  # Step 1: Count k-mers
  log_msg(sprintf("[1/5] Counting %d-mers across %d sequences...",
                  k, length(seqs_char)))
  kmer_counts <- count_kmers_dataset(rna_stringset, k)
  log_msg(sprintf("      Found %d unique %d-mers", nrow(kmer_counts), k))

  if (nrow(kmer_counts) == 0) {
    stop(sprintf("No valid %d-mers found in dataset.", k))
  }

  # Step 2: Filter
  log_msg("[2/5] Filtering by frequency and count thresholds...")
  kmer_filtered <- kmer_counts %>%
    filter(seq_frequency >= min_freq, total_count >= min_count)
  log_msg(sprintf("      %d k-mers pass filters", nrow(kmer_filtered)))

  if (nrow(kmer_filtered) == 0) {
    warning("No k-mers pass the current filters. Relaxing to top 20 by count.")
    kmer_filtered <- head(kmer_counts, 20)
  }

  # Step 3: Enrichment
  log_msg("[3/5] Calculating enrichment scores...")
  kmer_enriched <- calculate_enrichment_scores(kmer_filtered, total_positions)

  # Step 4: P-values
  log_msg("[4/5] Computing p-values (binomial test + Bonferroni)...")
  kmer_scored <- calculate_pvalues(kmer_enriched, total_positions)

  # Step 5: Rank
  log_msg("[5/5] Ranking motifs...")
  motif_ranked <- rank_motifs(kmer_scored)

  top_motifs   <- head(motif_ranked, top_n)
  top_kmer_vec <- top_motifs$kmer

  # Positional distribution for top motifs
  log_msg("      Computing positional distributions for top motifs...")
  pos_dist <- tryCatch(
    compute_positional_distribution(rna_stringset, top_kmer_vec, n_pos_bins),
    error = function(e) {
      warning(sprintf("Positional distribution failed: %s", conditionMessage(e)))
      data.frame()
    }
  )

  # Build PWMs for top motifs (group similar kmers)
  log_msg("      Building PWMs for top motifs...")
  pwm_list <- lapply(seq_len(min(5, nrow(top_motifs))), function(i) {
    km <- top_motifs$kmer[i]
    build_pwm_from_kmers(c(km))
  })
  pwm_list <- Filter(Negate(is.null), pwm_list)

  log_msg("=== Motif Discovery Complete ===")

  list(
    motif_table      = motif_ranked,
    top_motifs       = top_motifs,
    positional_dist  = pos_dist,
    pwm_list         = pwm_list,
    params = list(
      k          = k,
      min_freq   = min_freq,
      min_count  = min_count,
	      top_n      = top_n,
	      total_kmers_found = nrow(kmer_counts),
	      kmers_after_filter = nrow(kmer_filtered),
	      significant_motifs = sum(motif_ranked$significant, na.rm = TRUE)
	    ),
    log = log
  )
}
