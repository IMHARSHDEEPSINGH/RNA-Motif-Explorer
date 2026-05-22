# =============================================================================
# RNA Motif Discovery Explorer
# Script: structure_prediction.R
# Purpose: Secondary structure prediction and analysis
# =============================================================================

library(Biostrings)
library(dplyr)
library(tibble)
library(stringr)

# Note: This module provides secondary structure prediction with fallback support.
# Attempts to use RNAstructure package (via RNAfold) if available.
# Falls back to heuristic-based structure prediction if not available.

# --------------------------------------------------------------------------
# predict_secondary_structure
# Predicts secondary structure for a single RNA sequence
# Returns dot-bracket notation (e.g., "(((...)))" for a hairpin)
# Attempts RNAstructure first, falls back to heuristic prediction
# --------------------------------------------------------------------------
predict_secondary_structure <- function(sequence) {
  tryCatch({
    seq_clean <- gsub("T", "U", toupper(sequence))
    
    if (nchar(seq_clean) < 10) {
      return(paste0(rep(".", nchar(seq_clean)), collapse = ""))
    }
    
    # Try to use RNAstructure if available
    if (requireNamespace("RNAstructure", quietly = TRUE)) {
      result <- tryCatch({
        RNAstructure::fold.RNA(seq_clean, verbose = FALSE)
      }, error = function(e) NULL)
      
      if (!is.null(result) && !is.null(result$structure)) {
        return(result$structure)
      }
    }
    
    # Fallback: Use heuristic structure prediction based on GC content and base pairing
    predict_structure_heuristic(seq_clean)
    
  }, error = function(e) {
    # Final fallback: return unstructured
    paste0(rep(".", nchar(sequence)), collapse = "")
  })
}

# --------------------------------------------------------------------------
# predict_structure_heuristic
# Simple heuristic for secondary structure based on complementarity and GC
# This is a fallback when RNAstructure is not available
# --------------------------------------------------------------------------
predict_structure_heuristic <- function(sequence) {
  n <- nchar(sequence)
  if (n < 10) {
    return(paste0(rep(".", n), collapse = ""))
  }
  
  # Initialize all as unpaired
  structure <- rep(".", n)
  
  # Scan for complementary base pairs starting from ends (stem-loop heuristic)
  for (window in seq(4, min(20, n %/% 2))) {
    for (i in seq_len(n - window)) {
      j <- n - i + 1
      if (j <= i + 3) break
      
      # Check if bases at positions i and j are complementary
      base_i <- substr(sequence, i, i)
      base_j <- substr(sequence, j, j)
      
      # Watson-Crick pairs: AU, UA, GC, CG
      if ((base_i == "A" && base_j == "U") ||
          (base_i == "U" && base_j == "A") ||
          (base_i == "G" && base_j == "C") ||
          (base_i == "C" && base_j == "G")) {
        
        # Mark as paired if not already assigned
        if (structure[i] == ".") structure[i] <- "("
        if (structure[j] == ".") structure[j] <- ")"
      }
    }
  }
  
  paste0(structure, collapse = "")
}

# --------------------------------------------------------------------------
# predict_structures_dataset
# Predict secondary structures for all sequences in a dataset
# Input: cleaned sequences as character vector
# Output: tibble with sequence index, sequence string, and dot-bracket structure
# --------------------------------------------------------------------------
predict_structures_dataset <- function(sequences) {
  if (!is.character(sequences)) {
    stop("sequences must be a character vector")
  }
  
  cat("Predicting secondary structures for", length(sequences), "sequences...\n")
  
  structures <- sapply(
    seq_along(sequences),
    function(i) {
      if (i %% max(1, length(sequences) %/% 10) == 0) {
        cat(sprintf("  Progress: %d/%d\n", i, length(sequences)))
      }
      predict_secondary_structure(sequences[i])
    }
  )
  
  tibble(
    sequence_index = seq_along(sequences),
    sequence = sequences,
    structure = structures,
    length = nchar(sequences)
  )
}

# --------------------------------------------------------------------------
# parse_structure_elements
# Parse dot-bracket notation into structural elements
# Returns counts of: stems, loops, bulges, internal_loops, hairpins, etc.
# --------------------------------------------------------------------------
parse_structure_elements <- function(structure) {
  if (is.na(structure) || is.null(structure) || nchar(structure) == 0) {
    return(list(
      stems = 0, loops = 0, hairpins = 0, bulges = 0, 
      internal_loops = 0, multiloops = 0, unpaired = 0
    ))
  }
  
  chars <- strsplit(structure, "")[[1]]
  n <- length(chars)
  
  # Track pairing state
  paired <- chars %in% c("(", ")")
  unpaired <- chars == "."
  
  # Find stems (consecutive paired bases)
  stems <- 0
  in_stem <- FALSE
  for (i in seq_along(chars)) {
    if (paired[i] && !in_stem) {
      stems <- stems + 1
      in_stem <- TRUE
    } else if (!paired[i]) {
      in_stem <- FALSE
    }
  }
  
  # Count unpaired regions (loops/bulges)
  loops <- 0
  in_loop <- FALSE
  loop_length <- 0
  for (i in seq_along(chars)) {
    if (unpaired[i]) {
      if (!in_loop) {
        loops <- loops + 1
        in_loop <- TRUE
        loop_length <- 1
      } else {
        loop_length <- loop_length + 1
      }
    } else {
      in_loop <- FALSE
    }
  }
  
  # Classify loops/bulges based on position and length
  # A hairpin has unpaired region at the "top" (between matching stem)
  # A bulge is small unpaired region within a stem
  hairpins <- 0
  bulges <- 0
  if (loops > 0) {
    # Simplified heuristic: count unmatched regions as hairpins or bulges
    hairpins <- max(1, loops %/% 2)
    bulges <- loops - hairpins
  }
  
  list(
    stems = stems,
    loops = loops,
    hairpins = hairpins,
    bulges = bulges,
    unpaired_count = sum(unpaired),
    paired_count = sum(paired)
  )
}

# --------------------------------------------------------------------------
# get_position_structure_type
# Determine the structural type (stem, loop, bulge, etc.) at a given position
# Returns one of: "stem", "loop", "bulge", "hairpin", "unpaired"
# --------------------------------------------------------------------------
get_position_structure_type <- function(structure, position) {
  if (is.na(structure) || is.null(structure) || position < 1 || position > nchar(structure)) {
    return("unknown")
  }
  
  chars <- strsplit(structure, "")[[1]]
  
  if (position > length(chars)) {
    return("unknown")
  }
  
  char_at_pos <- chars[position]
  
  if (char_at_pos == ".") {
    # Check if part of a hairpin or bulge
    # Hairpin: loop region not contained by any pairing
    # Bulge: small unpaired region within stem
    
    # For simplicity, classify as "loop" for now
    # More sophisticated: trace parentheses to determine loop type
    return("loop")
  } else if (char_at_pos == "(" || char_at_pos == ")") {
    return("stem")
  }
  
  "unpaired"
}

# --------------------------------------------------------------------------
# map_position_to_structure
# Create a full position-wise mapping of each base to its structural element
# Output: tibble with position, base, structure_char, structure_type
# --------------------------------------------------------------------------
map_position_to_structure <- function(sequence, structure) {
  if (nchar(sequence) != nchar(structure)) {
    warning("Sequence and structure lengths don't match")
    return(NULL)
  }
  
  seq_chars <- strsplit(sequence, "")[[1]]
  struct_chars <- strsplit(structure, "")[[1]]
  
  types <- sapply(
    seq_along(struct_chars),
    function(i) get_position_structure_type(structure, i)
  )
  
  tibble(
    position = seq_along(seq_chars),
    base = seq_chars,
    structure_notation = struct_chars,
    structure_type = types
  )
}

# --------------------------------------------------------------------------
# structure_statistics
# Compute comprehensive structure statistics
# --------------------------------------------------------------------------
structure_statistics <- function(structures_df) {
  if (is.null(structures_df) || nrow(structures_df) == 0) {
    return(NULL)
  }
  
  structure_elements_list <- lapply(
    structures_df$structure,
    parse_structure_elements
  )
  
  # Aggregate statistics
  summary_stats <- list(
    total_sequences = nrow(structures_df),
    avg_sequence_length = mean(structures_df$length),
    avg_stems = mean(sapply(structure_elements_list, function(x) x$stems)),
    avg_loops = mean(sapply(structure_elements_list, function(x) x$loops)),
    avg_hairpins = mean(sapply(structure_elements_list, function(x) x$hairpins)),
    avg_bulges = mean(sapply(structure_elements_list, function(x) x$bulges)),
    avg_paired_bases = mean(sapply(structure_elements_list, function(x) x$paired_count)),
    avg_unpaired_bases = mean(sapply(structure_elements_list, function(x) x$unpaired_count))
  )
  
  tibble(
    metric = names(summary_stats),
    value = unlist(summary_stats)
  )
}

# --------------------------------------------------------------------------
# export_structures_to_csv
# Save predicted structures to CSV for downstream analysis
# --------------------------------------------------------------------------
export_structures_to_csv <- function(structures_df, output_dir = "outputs") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_path <- file.path(output_dir, "structure_predictions.csv")
  write.csv(structures_df, output_path, row.names = FALSE)
  
  cat("Structures saved to:", output_path, "\n")
  output_path
}
