# =============================================================================
# Tests for Secondary Structure Prediction Module
# =============================================================================

library(testthat)
source("../../scripts/structure_prediction.R")
source("../../scripts/motif_structure_analysis.R")

context("Structure Prediction Functions")

# Test data
test_sequence_short <- "AUGCGUACG"
test_sequence_long <- paste0(
  "GGUAUUAGCGCCAGAUGUAAUCCGUAACCCGC",
  "AUGGAAACGGGCACUGAACAUCGAUGUAAGUA"
)

test_sequences <- c(
  "AUGCGUACGAUGCUAGCUAGCU",
  "GCUAGCUAGCUAGCUAGCUAGC",
  "AUGAUGAUGAUGAUGAUGAUGA"
)

# ---- Tests for predict_secondary_structure ----
test_that("predict_secondary_structure returns valid dot-bracket notation", {
  structure <- predict_secondary_structure(test_sequence_long)
  
  # Check that output is character
  expect_is(structure, "character")
  
  # Check length matches input
  expect_equal(nchar(structure), nchar(test_sequence_long))
  
  # Check valid characters only
  expect_true(all(grepl("^[().]*$", structure)))
})

test_that("predict_secondary_structure handles short sequences", {
  structure <- predict_secondary_structure(test_sequence_short)
  
  # Should return a structure of same length
  expect_equal(nchar(structure), nchar(test_sequence_short))
})

test_that("predict_secondary_structure handles very short sequences", {
  structure <- predict_secondary_structure("AU")
  
  expect_equal(nchar(structure), 2)
  expect_true(grepl("^[().]*$", structure))
})

# ---- Tests for parse_structure_elements ----
test_that("parse_structure_elements counts stems correctly", {
  # Simple hairpin structure: ((())) = 3 stems
  structure <- "(((...)))"
  elements <- parse_structure_elements(structure)
  
  expect_is(elements, "list")
  expect_true(elements$stems > 0)
  expect_equal(elements$paired_count, 6)
  expect_equal(elements$unpaired_count, 3)
})

test_that("parse_structure_elements handles all dots (unstructured)", {
  structure <- "........."
  elements <- parse_structure_elements(structure)
  
  expect_equal(elements$stems, 0)
  expect_equal(elements$paired_count, 0)
  expect_equal(elements$unpaired_count, 9)
})

test_that("parse_structure_elements handles NA/NULL", {
  elements_na <- parse_structure_elements(NA)
  elements_null <- parse_structure_elements(NULL)
  
  expect_equal(elements_na$stems, 0)
  expect_equal(elements_null$stems, 0)
})

# ---- Tests for get_position_structure_type ----
test_that("get_position_structure_type identifies stems and loops", {
  structure <- "(((...)))"
  
  # Position 1 should be stem (opening paren)
  type_stem <- get_position_structure_type(structure, 1)
  expect_equal(type_stem, "stem")
  
  # Position 5 should be loop (dot)
  type_loop <- get_position_structure_type(structure, 5)
  expect_equal(type_loop, "loop")
})

test_that("get_position_structure_type handles invalid positions", {
  structure <- "(((...)))"
  
  # Out of bounds
  type_invalid <- get_position_structure_type(structure, 100)
  expect_equal(type_invalid, "unknown")
  
  # Negative
  type_neg <- get_position_structure_type(structure, -1)
  expect_equal(type_neg, "unknown")
})

# ---- Tests for map_position_to_structure ----
test_that("map_position_to_structure creates valid mapping", {
  sequence <- "AUGCUAGCUA"
  structure <- "(((...)))"
  
  mapping <- map_position_to_structure(sequence, structure)
  
  expect_is(mapping, "data.frame")
  expect_equal(nrow(mapping), nchar(sequence))
  expect_true(all(c("position", "base", "structure_notation", 
                   "structure_type") %in% names(mapping)))
})

test_that("map_position_to_structure handles length mismatch", {
  sequence <- "AUGCUAGCUA"
  structure <- "(((.)))"  # Different length
  
  expect_warning(
    result <- map_position_to_structure(sequence, structure),
    "lengths don't match"
  )
})

# ---- Tests for predict_structures_dataset ----
test_that("predict_structures_dataset processes multiple sequences", {
  structures_df <- predict_structures_dataset(test_sequences)
  
  expect_is(structures_df, "data.frame")
  expect_equal(nrow(structures_df), length(test_sequences))
  expect_true(all(c("sequence_index", "sequence", "structure", 
                    "length") %in% names(structures_df)))
})

test_that("predict_structures_dataset returns valid structures", {
  structures_df <- predict_structures_dataset(test_sequences)
  
  # All structures should match sequence length
  for (i in seq_len(nrow(structures_df))) {
    expect_equal(
      nchar(structures_df$sequence[i]),
      nchar(structures_df$structure[i])
    )
  }
})

# ---- Tests for structure_statistics ----
test_that("structure_statistics computes valid metrics", {
  structures_df <- predict_structures_dataset(test_sequences)
  stats <- structure_statistics(structures_df)
  
  expect_is(stats, "data.frame")
  expect_true("metric" %in% names(stats))
  expect_true("value" %in% names(stats))
  expect_true(nrow(stats) > 0)
})

test_that("structure_statistics handles NULL input", {
  stats <- structure_statistics(NULL)
  expect_null(stats)
  
  stats_empty <- structure_statistics(data.frame())
  expect_null(stats_empty)
})

# ---- Tests for Motif-Structure Mapping ----
test_that("find_motif_positions_in_sequence finds all occurrences", {
  sequence <- "AUGAUGAUG"
  motif <- "AUG"
  
  positions <- find_motif_positions_in_sequence(sequence, motif)
  
  expect_is(positions, "integer")
  expect_equal(length(positions), 3)  # Should find AUG three times
  expect_equal(positions, c(1, 4, 7))
})

test_that("find_motif_positions_in_sequence handles T/U conversion", {
  sequence_dna <- "ATGATGATG"
  sequence_rna <- "AUGAUGAUG"
  motif <- "AUG"
  
  positions_dna <- find_motif_positions_in_sequence(sequence_dna, motif)
  positions_rna <- find_motif_positions_in_sequence(sequence_rna, motif)
  
  # Both should find same positions
  expect_equal(positions_dna, positions_rna)
})

test_that("find_motif_positions_in_sequence handles missing motif", {
  sequence <- "AUGCGUACG"
  motif <- "GGGG"
  
  positions <- find_motif_positions_in_sequence(sequence, motif)
  
  expect_equal(length(positions), 0)
})

context("Structure Analysis Integration")

test_that("export_structures_to_csv creates valid CSV", {
  structures_df <- predict_structures_dataset(test_sequences)
  output_dir <- tempdir()
  
  result_path <- export_structures_to_csv(structures_df, output_dir)
  
  expect_true(file.exists(result_path))
  expect_true(grepl("structure_predictions.csv", result_path))
  
  # Verify CSV can be read back
  read_back <- read.csv(result_path)
  expect_equal(nrow(read_back), nrow(structures_df))
})

test_that("export_structures_to_csv creates directory if needed", {
  structures_df <- predict_structures_dataset(test_sequences)
  output_dir <- file.path(tempdir(), "new_test_dir", "nested")
  
  result_path <- export_structures_to_csv(structures_df, output_dir)
  
  expect_true(file.exists(result_path))
  unlink(file.path(tempdir(), "new_test_dir"), recursive = TRUE)
})

context("Edge Cases and Error Handling")

test_that("Structure functions handle empty sequences gracefully", {
  expect_error(
    predict_structures_dataset(character(0)),
    NA  # Should not error
  )
})

test_that("Structure functions handle sequences with N bases", {
  sequence_with_n <- "AUGCNAUGCN"
  structure <- predict_secondary_structure(sequence_with_n)
  
  # Should complete without error
  expect_equal(nchar(structure), nchar(sequence_with_n))
})

test_that("get_position_structure_type is robust to edge cases", {
  # Empty structure
  type <- get_position_structure_type("", 1)
  expect_equal(type, "unknown")
  
  # Single character
  type_single <- get_position_structure_type(".", 1)
  expect_equal(type_single, "loop")
})
