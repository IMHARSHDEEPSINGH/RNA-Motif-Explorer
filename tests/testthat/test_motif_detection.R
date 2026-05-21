library(testthat)
library(Biostrings)

source(file.path("..", "..", "scripts", "motif_detection.R"))

# Use a small dataset to verify Biostrings-based k-mer counting and scoring.
test_that("count_kmers_dataset returns correct counts", {
  seqs <- RNAStringSet(c("AUGCAUGC", "AUGCUA"))
  result <- count_kmers_dataset(seqs, 3)

  expect_true("AUG" %in% result$kmer)
  expect_equal(result$total_count[result$kmer == "AUG"], 3)
  expect_equal(result$seq_count[result$kmer == "AUG"], 2)
})

test_that("motif_discovery_pipeline supports control dataset comparisons", {
  fg <- RNAStringSet(c("AUGCGA", "AUGCGA"))
  ctrl <- RNAStringSet(c("AUGCGU", "AUGCGU"))

  result <- motif_discovery_pipeline(
    rna_stringset         = fg,
    k                     = 3,
    min_freq              = 0.1,
    min_count             = 1,
    top_n                 = 5,
    background_model      = "control",
    p_adjust_method       = "fdr",
    control_rna_stringset = ctrl,
    verbose               = FALSE
  )

  expect_true(result$params$control_dataset)
  expect_equal(result$params$background_model, "control")
  expect_equal(result$params$p_adjust_method, "fdr")
  expect_true("total_count_ctrl" %in% names(result$motif_table))
  expect_true(all(!is.na(result$motif_table$p_adjusted)))
})

test_that("motif_discovery_pipeline supports mononucleotide background and FDR", {
  seqs <- RNAStringSet(c("AUGCAUGC", "AUGCAUGC"))
  result <- motif_discovery_pipeline(
    rna_stringset    = seqs,
    k                = 3,
    min_freq         = 0,
    min_count        = 1,
    top_n            = 5,
    background_model = "mononucleotide",
    p_adjust_method  = "fdr",
    verbose          = FALSE
  )

  expect_equal(result$params$background_model, "mononucleotide")
  expect_equal(result$params$p_adjust_method, "fdr")
  expect_true(all(result$motif_table$p_adjusted <= 1))
})
