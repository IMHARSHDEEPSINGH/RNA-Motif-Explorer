#!/usr/bin/env Rscript

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop("Please install testthat first: install.packages('testthat')")
}

testthat::test_dir("tests/testthat")
