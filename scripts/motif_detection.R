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
# Uses Biostrings::oligonucleotideFrequency for fast counting.
# Returns a tibble sorted by total_count descending.
# -----------------------------------------------------------------------------
count_kmers_dataset <- function(rna_stringset, k) {
  seqs_char <- as.character(rna_stringset)
  n_seqs    <- length(seqs_char)

  if (n_seqs == 0 || k <= 0) {
    return(tibble(kmer = character(0), total_count = integer(0),
                  seq_count = integer(0), seq_frequency = numeric(0)))
  }

  counts <- oligonucleotideFrequency(rna_stringset,
                                     width = k,
                                     step = 1,
                                     as.prob = FALSE,
                                     with.labels = TRUE)
  counts <- as.matrix(counts)
  kmer_names <- colnames(counts)

  if (is.null(kmer_names) || length(kmer_names) == 0) {
    return(tibble(kmer = character(0), total_count = integer(0),
                  seq_count = integer(0), seq_frequency = numeric(0)))
  }

  valid <- grepl("^[AUGC]+$", kmer_names)
  if (!any(valid)) {
    return(tibble(kmer = character(0), total_count = integer(0),
                  seq_count = integer(0), seq_frequency = numeric(0)))
  }

  counts <- counts[, valid, drop = FALSE]
  kmer_names <- colnames(counts)

  total_count <- colSums(counts)
  seq_count   <- colSums(counts > 0)

  tibble(
    kmer = kmer_names,
    total_count = as.integer(total_count),
    seq_count = as.integer(seq_count),
    seq_frequency = seq_count / n_seqs
  ) %>%
    arrange(desc(total_count))
}

calculate_mono_background_probabilities <- function(rna_stringset, kmers) {
  seqs <- normalize_rna_text(as.character(rna_stringset))
  total_nt <- sum(nchar(seqs))
  if (total_nt == 0) {
    stop("No nucleotides available to calculate mononucleotide background.")
  }

  freqs <- c(
    A = sum(str_count(seqs, "A")) / total_nt,
    U = sum(str_count(seqs, "U")) / total_nt,
    G = sum(str_count(seqs, "G")) / total_nt,
    C = sum(str_count(seqs, "C")) / total_nt
  )
  freqs[freqs <= 0] <- 1e-10

  vapply(kmers, function(km) {
    nts <- strsplit(km, "")[[1]]
    prod(freqs[nts])
  }, numeric(1), USE.NAMES = FALSE)
}

merge_control_counts <- function(fg_counts, control_counts) {
  if (is.null(control_counts) || nrow(control_counts) == 0) {
    return(fg_counts)
  }
  merged <- full_join(fg_counts, control_counts, by = "kmer",
                      suffix = c("", "_ctrl"))
  merged %>%
    mutate(
      total_count_ctrl = ifelse(is.na(total_count_ctrl), 0L, total_count_ctrl),
      seq_count_ctrl   = ifelse(is.na(seq_count_ctrl), 0L, seq_count_ctrl),
      seq_frequency_ctrl = ifelse(is.na(seq_frequency_ctrl), 0, seq_frequency_ctrl)
    )
}

# -----------------------------------------------------------------------------
# calculate_enrichment_scores
# Computes enrichment relative to a chosen background model.
# Background options: uniform, mononucleotide, control.
# -----------------------------------------------------------------------------
calculate_enrichment_scores <- function(kmer_counts, total_positions,
                                        background_model = c("uniform",
                                                             "mononucleotide",
                                                             "control"),
                                        rna_stringset = NULL,
                                        control_total_positions = NULL) {
  if (nrow(kmer_counts) == 0) {
    return(kmer_counts %>%
             mutate(observed_freq = numeric(0), expected_freq = numeric(0),
                    expected_count = numeric(0), enrichment_score = numeric(0),
                    control_freq = numeric(0)))
  }

  background_model <- match.arg(background_model)
  k <- nchar(kmer_counts$kmer[1])
  eps <- 1e-10

  expected_freq <- switch(background_model,
    uniform = rep(0.25^k, nrow(kmer_counts)),
    mononucleotide = {
      if (is.null(rna_stringset)) {
        stop("rna_stringset is required for mononucleotide background model.")
      }
      calculate_mono_background_probabilities(rna_stringset, kmer_counts$kmer)
    },
    control = {
      if (is.null(control_total_positions) || !"total_count_ctrl" %in% colnames(kmer_counts)) {
        stop("Control background requires control dataset counts and total positions.")
      }
      control_total_positions <- max(control_total_positions, 1)
      pmax(kmer_counts$total_count_ctrl / control_total_positions, eps)
    }
  )

  observed_freq <- kmer_counts$total_count / total_positions
  control_freq <- if ("total_count_ctrl" %in% colnames(kmer_counts) &&
                       !is.null(control_total_positions)) {
    kmer_counts$total_count_ctrl / max(control_total_positions, 1)
  } else {
    rep(NA_real_, nrow(kmer_counts))
  }

  kmer_counts %>%
    mutate(
      observed_freq  = observed_freq,
      expected_freq  = expected_freq,
      expected_count = expected_freq * total_positions,
      enrichment_score = log2((observed_freq + eps) / (expected_freq + eps)),
      control_freq   = control_freq
    )
}

# -----------------------------------------------------------------------------
# calculate_pvalues
# Binomial test for a single dataset, Fisher's exact test when control is supplied.
# Uses Bonferroni or FDR correction.
# -----------------------------------------------------------------------------
calculate_pvalues <- function(kmer_enriched, total_positions,
                              p_adjust_method = c("bonferroni", "fdr"),
                              control_total_positions = NULL) {
  if (nrow(kmer_enriched) == 0) {
    return(kmer_enriched %>%
             mutate(p_value = numeric(0), p_adjusted = numeric(0),
                    significant = logical(0), neg_log10_p = numeric(0)))
  }

  p_adjust_method <- match.arg(p_adjust_method)
  use_control <- !is.null(control_total_positions) &&
    "total_count_ctrl" %in% colnames(kmer_enriched)

  if (use_control) {
    control_total_positions <- max(control_total_positions, 1)
    pvals <- mapply(function(obs, ctrl_obs) {
      tryCatch(
        fisher.test(matrix(c(obs,
                             total_positions - obs,
                             ctrl_obs,
                             control_total_positions - ctrl_obs), nrow = 2),
                    alternative = "greater")$p.value,
        error = function(e) 1.0
      )
    },
    kmer_enriched$total_count,
    kmer_enriched$total_count_ctrl,
    SIMPLIFY = TRUE)
  } else {
    pvals <- mapply(function(obs, p) {
      tryCatch(
        binom.test(obs, total_positions, p = p, alternative = "greater")$p.value,
        error = function(e) 1.0
      )
    },
    kmer_enriched$total_count,
    kmer_enriched$expected_freq,
    SIMPLIFY = TRUE)
  }

  kmer_enriched %>%
    mutate(
      p_value     = pvals,
      p_adjusted  = p.adjust(pvals, method = p_adjust_method),
      significant = p_adjusted < 0.05,
      neg_log10_p = -log10(p_value + 1e-300)
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
# -> positional distribution -> PWMs for top motifs.
# Supports optional control dataset comparison and background model selection.
# -----------------------------------------------------------------------------
motif_discovery_pipeline <- function(rna_stringset,
                                      k           = 6,
                                      min_freq    = 0.1,
                                      min_count   = 3,
                                      top_n       = 20,
                                      n_pos_bins  = 20,
                                      background_model = c("uniform",
                                                           "mononucleotide",
                                                           "control"),
                                      p_adjust_method = c("bonferroni", "fdr"),
                                      control_rna_stringset = NULL,
                                      verbose     = TRUE) {
  log <- character(0)
  log_msg <- function(msg) {
    if (verbose) message(msg)
    log <<- c(log, msg)
  }

  background_model <- match.arg(background_model)
  p_adjust_method <- match.arg(p_adjust_method)

  log_msg(sprintf("=== Motif Discovery Pipeline (k=%d) ===", k))

  if (length(rna_stringset) == 0) {
    stop("No sequences supplied for motif discovery.")
  }

  seqs_char       <- as.character(rna_stringset)
  total_positions <- sum(pmax(nchar(seqs_char) - k + 1, 0))
  total_positions <- max(total_positions, 1)

  control_total_positions <- NULL
  control_counts <- NULL
  if (!is.null(control_rna_stringset)) {
    control_seqs_char <- as.character(control_rna_stringset)
    control_total_positions <- sum(pmax(nchar(control_seqs_char) - k + 1, 0))
    control_total_positions <- max(control_total_positions, 1)
    control_counts <- count_kmers_dataset(control_rna_stringset, k)
  }

  if (background_model == "control" && is.null(control_rna_stringset)) {
    warning("Control background model selected but no control dataset supplied; falling back to uniform background.")
    background_model <- "uniform"
  }

  log_msg(sprintf("[1/5] Counting %d-mers across %d sequences...", k, length(seqs_char)))
  fg_counts <- count_kmers_dataset(rna_stringset, k)
  if (!is.null(control_counts)) {
    kmer_counts <- merge_control_counts(fg_counts, control_counts)
  } else {
    kmer_counts <- fg_counts
  }
  log_msg(sprintf("      Found %d unique %d-mers", nrow(kmer_counts), k))

  if (nrow(kmer_counts) == 0) {
    stop(sprintf("No valid %d-mers found in dataset.", k))
  }

  log_msg("[2/5] Filtering by frequency and count thresholds...")
  kmer_filtered <- kmer_counts %>%
    filter(seq_frequency >= min_freq, total_count >= min_count)
  log_msg(sprintf("      %d k-mers pass filters", nrow(kmer_filtered)))

  if (nrow(kmer_filtered) == 0) {
    warning("No k-mers pass the current filters. Relaxing to top 20 by count.")
    kmer_filtered <- head(kmer_counts, 20)
  }

  log_msg("[3/5] Calculating enrichment scores...")
  kmer_enriched <- calculate_enrichment_scores(
    kmer_filtered,
    total_positions,
    background_model = background_model,
    rna_stringset = rna_stringset,
    control_total_positions = control_total_positions
  )

  log_msg("[4/5] Computing p-values...")
  kmer_scored <- calculate_pvalues(
    kmer_enriched,
    total_positions,
    p_adjust_method = p_adjust_method,
    control_total_positions = control_total_positions
  )

  log_msg("[5/5] Ranking motifs...")
  motif_ranked <- rank_motifs(kmer_scored)

  top_motifs   <- head(motif_ranked, top_n)
  top_kmer_vec <- top_motifs$kmer

  log_msg("      Computing positional distributions for top motifs...")
  pos_dist <- tryCatch(
    compute_positional_distribution(rna_stringset, top_kmer_vec, n_pos_bins),
    error = function(e) {
      warning(sprintf("Positional distribution failed: %s", conditionMessage(e)))
      data.frame()
    }
  )

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
      k                 = k,
      min_freq          = min_freq,
      min_count         = min_count,
      top_n             = top_n,
      background_model  = background_model,
      p_adjust_method   = p_adjust_method,
      control_dataset   = !is.null(control_rna_stringset),
      total_kmers_found = nrow(kmer_counts),
      kmers_after_filter= nrow(kmer_filtered),
      significant_motifs= sum(motif_ranked$significant, na.rm = TRUE)
    ),
    log = log
  )
}
