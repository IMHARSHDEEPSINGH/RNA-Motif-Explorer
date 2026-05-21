# =============================================================================
# RNA Motif Discovery Explorer
# Script: preprocessing.R
# Purpose: FASTA ingestion, validation, cleaning, and summary statistics
# =============================================================================

library(Biostrings)
library(stringr)
library(dplyr)
library(tibble)

# -----------------------------------------------------------------------------
# read_fasta_file
# Reads a FASTA file and returns a named character vector of sequences.
# Accepts .fasta, .fa, or .txt extensions.
# -----------------------------------------------------------------------------
read_fasta_file <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("File not found: %s", filepath))
  }
  ext <- tolower(tools::file_ext(filepath))
  if (!ext %in% c("fasta", "fa", "txt")) {
    stop(sprintf("Unsupported file extension: .%s. Use .fasta, .fa, or .txt", ext))
  }
  tryCatch({
    seqs <- readRNAStringSet(filepath)
    if (length(seqs) == 0) {
      stop("No sequences found in the file.")
    }
    seqs
  }, error = function(e) {
    stop(sprintf("Error reading FASTA: %s", conditionMessage(e)))
  })
}

# -----------------------------------------------------------------------------
# validate_rna_sequences
# Checks for valid RNA characters (A, U, G, C, N).
# Returns a list with valid sequences and a report of removed entries.
# -----------------------------------------------------------------------------
validate_rna_sequences <- function(rna_stringset) {
  n_input      <- length(rna_stringset)
  seqs_char    <- as.character(rna_stringset)
  names_char   <- names(rna_stringset)

  # Allowed characters: A U G C N (case-insensitive handled after toupper)
  valid_pattern <- "^[AUGCN]+$"
  is_valid      <- grepl(valid_pattern, toupper(seqs_char))

  invalid_report <- tibble(
    name   = names_char[!is_valid],
    reason = "Contains invalid nucleotide characters"
  )

  valid_seqs <- rna_stringset[is_valid]

  list(
    valid_sequences  = valid_seqs,
    invalid_report   = invalid_report,
    n_input          = n_input,
    n_valid          = sum(is_valid),
    n_invalid        = sum(!is_valid)
  )
}

# -----------------------------------------------------------------------------
# clean_rna_sequences
# - Converts to uppercase
# - Replaces T with U (handles accidental DNA input)
# - Strips whitespace from sequence names
# - Removes sequences shorter than min_length
# -----------------------------------------------------------------------------
clean_rna_sequences <- function(rna_stringset, min_length = 10) {
  seqs_char  <- toupper(as.character(rna_stringset))
  seqs_char  <- gsub("T", "U", seqs_char)   # DNA -> RNA
  seqs_char  <- gsub("[^AUGCN]", "", seqs_char)  # strip remaining junk

  seq_names  <- str_trim(names(rna_stringset))

  # Filter by minimum length
  lengths    <- nchar(seqs_char)
  keep       <- lengths >= min_length

  cleaned <- RNAStringSet(seqs_char[keep])
  names(cleaned) <- seq_names[keep]

  removed <- tibble(
    name   = seq_names[!keep],
    length = lengths[!keep],
    reason = sprintf("Shorter than minimum length (%d nt)", min_length)
  )

  list(
    sequences     = cleaned,
    removed_short = removed,
    n_removed     = sum(!keep)
  )
}

# -----------------------------------------------------------------------------
# remove_duplicate_sequences
# Deduplicates by sequence content, retaining first occurrence.
# -----------------------------------------------------------------------------
remove_duplicate_sequences <- function(rna_stringset) {
  seqs_char <- as.character(rna_stringset)
  dup_mask  <- duplicated(seqs_char)

  unique_seqs <- rna_stringset[!dup_mask]

  dup_report <- tibble(
    name     = names(rna_stringset)[dup_mask],
    sequence = seqs_char[dup_mask]
  )

  list(
    sequences   = unique_seqs,
    duplicates  = dup_report,
    n_removed   = sum(dup_mask)
  )
}

# -----------------------------------------------------------------------------
# calculate_sequence_stats
# Computes per-sequence and dataset-level summary statistics.
# Returns a list with a per-sequence tibble and global summary.
# -----------------------------------------------------------------------------
calculate_sequence_stats <- function(rna_stringset) {
  seqs_char <- as.character(rna_stringset)
  seq_names <- names(rna_stringset)
  lengths   <- nchar(seqs_char)

  # Per-nucleotide counts
  count_nt <- function(seq, nt) str_count(seq, nt)

  per_seq <- tibble(
    name      = seq_names,
    length    = lengths,
    count_A   = count_nt(seqs_char, "A"),
    count_U   = count_nt(seqs_char, "U"),
    count_G   = count_nt(seqs_char, "G"),
    count_C   = count_nt(seqs_char, "C"),
    count_N   = count_nt(seqs_char, "N")
  ) %>%
    mutate(
      freq_A = count_A / length,
      freq_U = count_U / length,
      freq_G = count_G / length,
      freq_C = count_C / length,
      gc_content = (count_G + count_C) / length
    )

  summary_stats <- list(
    total_sequences  = length(seqs_char),
    mean_length      = mean(lengths),
    median_length    = median(lengths),
    min_length       = min(lengths),
    max_length       = max(lengths),
    sd_length        = sd(lengths),
    mean_gc          = mean(per_seq$gc_content),
    mean_freq_A      = mean(per_seq$freq_A),
    mean_freq_U      = mean(per_seq$freq_U),
    mean_freq_G      = mean(per_seq$freq_G),
    mean_freq_C      = mean(per_seq$freq_C),
    total_nucleotides = sum(lengths)
  )

  list(
    per_sequence  = per_seq,
    summary       = summary_stats
  )
}

# -----------------------------------------------------------------------------
# preprocess_pipeline
# Full pipeline: read -> validate -> clean -> deduplicate -> stats
# Returns a named list with all intermediate results for logging/reporting.
# -----------------------------------------------------------------------------
preprocess_pipeline <- function(filepath, min_length = 10, verbose = TRUE) {
  log <- character(0)
  log_msg <- function(msg) {
    if (verbose) message(msg)
    log <<- c(log, msg)
  }

  log_msg("=== RNA Preprocessing Pipeline Started ===")

  # Step 1: Read
  log_msg(sprintf("[1/5] Reading FASTA: %s", filepath))
  raw_seqs <- read_fasta_file(filepath)
  log_msg(sprintf("      Loaded %d sequences", length(raw_seqs)))

  # Step 2: Validate
  log_msg("[2/5] Validating sequences...")
  val_result <- validate_rna_sequences(raw_seqs)
  log_msg(sprintf("      Valid: %d | Invalid: %d",
                  val_result$n_valid, val_result$n_invalid))

  # Step 3: Clean
  log_msg("[3/5] Cleaning sequences...")
  clean_result <- clean_rna_sequences(val_result$valid_sequences,
                                       min_length = min_length)
  log_msg(sprintf("      Removed %d sequences below min length (%d nt)",
                  clean_result$n_removed, min_length))

  # Step 4: Deduplicate
  log_msg("[4/5] Removing duplicates...")
  dedup_result <- remove_duplicate_sequences(clean_result$sequences)
  log_msg(sprintf("      Removed %d duplicate sequences",
                  dedup_result$n_removed))

  # Step 5: Stats
  log_msg("[5/5] Calculating statistics...")
  final_seqs <- dedup_result$sequences
  stats      <- calculate_sequence_stats(final_seqs)
  log_msg(sprintf("      Final dataset: %d sequences, mean length %.1f nt",
                  stats$summary$total_sequences,
                  stats$summary$mean_length))

  log_msg("=== Preprocessing Complete ===")

  list(
    sequences       = final_seqs,
    stats           = stats,
    invalid_report  = val_result$invalid_report,
    removed_short   = clean_result$removed_short,
    duplicates      = dedup_result$duplicates,
    log             = log,
    params          = list(min_length = min_length, filepath = filepath)
  )
}
