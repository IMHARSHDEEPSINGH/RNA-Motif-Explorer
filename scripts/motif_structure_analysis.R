# =============================================================================
# RNA Motif Discovery Explorer
# Script: motif_structure_analysis.R
# Purpose: Map motifs to secondary structure elements and compute enrichment
# =============================================================================

library(dplyr)
library(tibble)
library(stringr)

# --------------------------------------------------------------------------
# find_motif_positions_in_sequence
# Find all positions where a motif occurs in a sequence
# Returns a vector of starting positions (1-indexed)
# --------------------------------------------------------------------------
find_motif_positions_in_sequence <- function(sequence, motif) {
  seq_upper <- toupper(gsub("T", "U", sequence))
  motif_upper <- toupper(gsub("T", "U", motif))
  
  k <- nchar(motif_upper)
  n <- nchar(seq_upper)
  
  if (n < k) return(integer(0))
  
  positions <- integer(0)
  for (i in seq_len(n - k + 1)) {
    if (substr(seq_upper, i, i + k - 1) == motif_upper) {
      positions <- c(positions, i)
    }
  }
  
  positions
}

# --------------------------------------------------------------------------
# classify_motif_positions
# For each motif occurrence, determine its structural context
# Returns a tibble with position and structure_type
# --------------------------------------------------------------------------
classify_motif_positions <- function(sequence, motif, position_structure_map) {
  if (is.null(position_structure_map) || nrow(position_structure_map) == 0) {
    return(NULL)
  }
  
  positions <- find_motif_positions_in_sequence(sequence, motif)
  k <- nchar(motif)
  
  if (length(positions) == 0) {
    return(NULL)
  }
  
  # For each position, get the structure types of all bases in the motif
  classifications <- lapply(positions, function(start_pos) {
    end_pos <- start_pos + k - 1
    
    # Get structure types for all positions covered by this motif
    covered_positions <- which(
      position_structure_map$position >= start_pos & 
      position_structure_map$position <= end_pos
    )
    
    if (length(covered_positions) == 0) {
      return(list(position = start_pos, structure_type = "unknown", 
                  structure_types_covered = NA))
    }
    
    structure_types <- position_structure_map$structure_type[covered_positions]
    
    # Determine primary structure type (most common)
    primary_type <- names(sort(table(structure_types), decreasing = TRUE))[1]
    
    list(
      position = start_pos,
      motif = motif,
      sequence_position = start_pos,
      structure_type = primary_type,
      structure_types_covered = paste(structure_types, collapse = ",")
    )
  })
  
  do.call(rbind, lapply(classifications, as.data.frame, stringsAsFactors = FALSE))
}

# --------------------------------------------------------------------------
# map_motifs_to_structures_dataset
# Create comprehensive motif-structure mapping across entire dataset
# Input: sequences, structures_df from structure_prediction, motif_table from discovery
# Output: tibble with each motif occurrence and its structural context
# --------------------------------------------------------------------------
map_motifs_to_structures_dataset <- function(sequences, structures_df, motif_table) {
  if (is.null(structures_df) || is.null(motif_table) || nrow(motif_table) == 0) {
    return(NULL)
  }

  motif_column <- if ("kmer" %in% names(motif_table)) {
    "kmer"
  } else if ("motif" %in% names(motif_table)) {
    "motif"
  } else {
    stop("motif_table must contain a 'kmer' or 'motif' column")
  }

  # Get unique motifs to analyze
  motifs <- unique(motif_table[[motif_column]])

  # Source structure prediction to get helper functions
  # (Assuming it's already loaded)

  all_mappings <- list()

  for (seq_idx in seq_along(sequences)) {
    sequence <- sequences[seq_idx]
    structure <- structures_df$structure[seq_idx]
    
    # Create position-to-structure map for this sequence
    pos_struct_map <- map_position_to_structure(sequence, structure)
    
    # Find motif occurrences in this sequence
    for (motif in motifs) {
      classifications <- classify_motif_positions(sequence, motif, pos_struct_map)
      
      if (!is.null(classifications)) {
        classifications$sequence_idx <- seq_idx
        classifications$sequence <- sequence
        classifications$structure <- structure
        all_mappings[[paste(seq_idx, motif, sep = "_")]] <- classifications
      }
    }
  }
  
  if (length(all_mappings) == 0) {
    return(NULL)
  }
  
  result <- do.call(rbind, all_mappings)
  rownames(result) <- NULL
  
  as_tibble(result)
}

# --------------------------------------------------------------------------
# calculate_motif_structure_enrichment
# Compute enrichment scores for each motif in each structural context
# Returns Fisher's exact test p-values and log-odds ratios
# --------------------------------------------------------------------------
calculate_motif_structure_enrichment <- function(motif_structure_map, motif_table) {
  if (is.null(motif_structure_map) || nrow(motif_structure_map) == 0) {
    return(NULL)
  }

  motif_column <- if ("kmer" %in% names(motif_table)) {
    "kmer"
  } else if ("motif" %in% names(motif_table)) {
    "motif"
  } else {
    stop("motif_table must contain a 'kmer' or 'motif' column")
  }

  motifs <- unique(motif_table[[motif_column]])
  structure_types <- unique(motif_structure_map$structure_type)
  structure_types <- structure_types[!is.na(structure_types)]

  enrichment_results <- list()

  for (motif in motifs) {
    motif_rows <- which(motif_structure_map$motif == motif)
    
    if (length(motif_rows) == 0) next
    
    for (struct_type in structure_types) {
      # Create 2x2 contingency table
      # a: motif in this structure type
      # b: motif NOT in this structure type
      # c: other motifs in this structure type
      # d: other motifs NOT in this structure type
      
      in_motif_in_struct <- sum(
        motif_structure_map$motif[motif_rows] == motif & 
        motif_structure_map$structure_type[motif_rows] == struct_type
      )
      
      in_motif_not_struct <- sum(
        motif_structure_map$motif[motif_rows] == motif & 
        motif_structure_map$structure_type[motif_rows] != struct_type
      )
      
      not_motif_in_struct <- sum(
        motif_structure_map$motif != motif & 
        motif_structure_map$structure_type == struct_type
      )
      
      not_motif_not_struct <- sum(
        motif_structure_map$motif != motif & 
        motif_structure_map$structure_type != struct_type
      )
      
      # Fisher's exact test
      contingency_table <- matrix(
        c(in_motif_in_struct, in_motif_not_struct,
          not_motif_in_struct, not_motif_not_struct),
        nrow = 2
      )
      
      test_result <- fisher.test(contingency_table)
      p_value <- test_result$p.value
      
      # Calculate log odds ratio
      log_odds <- log((in_motif_in_struct + 0.5) / (in_motif_not_struct + 0.5)) -
                  log((not_motif_in_struct + 0.5) / (not_motif_not_struct + 0.5))
      
      # Calculate observed vs expected proportion
      observed_in_struct <- in_motif_in_struct / (in_motif_in_struct + in_motif_not_struct + 1)
      expected_in_struct <- (in_motif_in_struct + not_motif_in_struct) / nrow(motif_structure_map)
      
      enrichment_results[[paste(motif, struct_type, sep = "_")]] <- list(
        motif = motif,
        structure_type = struct_type,
        count_in_structure = in_motif_in_struct,
        count_not_in_structure = in_motif_not_struct,
        p_value = p_value,
        log_odds_ratio = log_odds,
        observed_proportion = observed_in_struct,
        expected_proportion = expected_in_struct
      )
    }
  }
  
  if (length(enrichment_results) == 0) {
    return(NULL)
  }
  
  result_df <- do.call(rbind, lapply(enrichment_results, as.data.frame))
  rownames(result_df) <- NULL
  
  as_tibble(result_df) %>%
    mutate(
      padj = p.adjust(p_value, method = "BH"),
      significance = ifelse(padj < 0.05, "**", ifelse(p_value < 0.05, "*", ""))
    )
}

# --------------------------------------------------------------------------
# get_motif_structure_distribution
# Summarize which structures contain each motif
# Returns: motif -> structure type -> count and percentage
# --------------------------------------------------------------------------
get_motif_structure_distribution <- function(motif_structure_map) {
  if (is.null(motif_structure_map) || nrow(motif_structure_map) == 0) {
    return(NULL)
  }
  
  motif_structure_map %>%
    filter(!is.na(structure_type)) %>%
    group_by(motif, structure_type) %>%
    summarise(
      count = n(),
      .groups = "drop"
    ) %>%
    group_by(motif) %>%
    mutate(
      total = sum(count),
      percentage = round(100 * count / total, 2)
    ) %>%
    ungroup()
}

# --------------------------------------------------------------------------
# compute_structure_correlation_metrics
# Calculate correlation between motif frequency and structural features
# Returns: correlation matrix motif x structure_type
# --------------------------------------------------------------------------
compute_structure_correlation_metrics <- function(motif_structure_map, structures_df) {
  if (is.null(motif_structure_map) || nrow(motif_structure_map) == 0) {
    return(NULL)
  }
  
  # Build a presence/absence matrix: sequences x motifs x structure types
  motifs <- unique(motif_structure_map$motif)
  struct_types <- unique(motif_structure_map$structure_type)
  struct_types <- struct_types[!is.na(struct_types)]
  sequences <- unique(motif_structure_map$sequence_idx)
  
  correlations <- list()
  
  for (motif in motifs) {
    motif_data <- motif_structure_map %>% filter(motif == !!motif)
    
    for (struct_type in struct_types) {
      # For each sequence, does this motif appear in this structure type?
      presence_by_seq <- sapply(sequences, function(seq_id) {
        sum(motif_data$sequence_idx == seq_id & motif_data$structure_type == struct_type) > 0
      })
      
      # For each sequence, how much of this structure type is present?
      struct_content_by_seq <- sapply(sequences, function(seq_id) {
        if (seq_id > nrow(structures_df)) return(0)
        structure <- structures_df$structure[seq_id]
        if (is.na(structure)) return(0)
        
        # Count unpaired bases in structure for loop context
        sum(gregexpr("\\.", structure)[[1]] > 0) / nchar(structure)
      })
      
      # Calculate point-biserial correlation
      if (length(unique(presence_by_seq)) > 1 && sd(struct_content_by_seq) > 0) {
        corr <- cor(as.numeric(presence_by_seq), struct_content_by_seq, method = "pearson")
      } else {
        corr <- 0
      }
      
      correlations[[paste(motif, struct_type, sep = "_")]] <- list(
        motif = motif,
        structure_type = struct_type,
        correlation = corr
      )
    }
  }
  
  if (length(correlations) == 0) {
    return(NULL)
  }
  
  result <- do.call(rbind, lapply(correlations, as.data.frame))
  rownames(result) <- NULL
  
  as_tibble(result) %>%
    arrange(desc(abs(correlation)))
}

# --------------------------------------------------------------------------
# export_motif_structure_analysis
# Save all analysis results to CSV files
# --------------------------------------------------------------------------
export_motif_structure_analysis <- function(
  motif_structure_map, 
  enrichment_results, 
  distribution_results,
  correlation_results,
  output_dir = "outputs"
) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results_list <- list()
  
  if (!is.null(motif_structure_map)) {
    path <- file.path(output_dir, "motif_structure_mapping.csv")
    write.csv(motif_structure_map, path, row.names = FALSE)
    results_list[["mapping"]] <- path
  }
  
  if (!is.null(enrichment_results)) {
    path <- file.path(output_dir, "structure_enrichment.csv")
    write.csv(enrichment_results, path, row.names = FALSE)
    results_list[["enrichment"]] <- path
  }
  
  if (!is.null(distribution_results)) {
    path <- file.path(output_dir, "motif_structure_distribution.csv")
    write.csv(distribution_results, path, row.names = FALSE)
    results_list[["distribution"]] <- path
  }
  
  if (!is.null(correlation_results)) {
    path <- file.path(output_dir, "structure_correlation.csv")
    write.csv(correlation_results, path, row.names = FALSE)
    results_list[["correlation"]] <- path
  }
  
  invisible(results_list)
}
