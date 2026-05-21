library(testthat)
library(Biostrings)

source(file.path("..", "..", "scripts", "preprocessing.R"))

test_that("preprocess_pipeline reads and cleans FASTA correctly", {
  fasta <- tempfile(fileext = ".fasta")
  cat(">seq1\nAUGCAUGC\n>seq2\nAUGNXX\n>seq3\nAUGC\n", file = fasta)

  result <- preprocess_pipeline(fasta, min_length = 4, verbose = FALSE)
  expect_equal(result$stats$summary$total_sequences, 2)
  expect_equal(length(result$sequences), 2)
  expect_true(all(grepl("^[AUGC]+$", as.character(result$sequences))))
  expect_false(any(grepl("N", as.character(result$sequences))))
})
