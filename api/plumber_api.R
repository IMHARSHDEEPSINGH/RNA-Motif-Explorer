# =============================================================================
# RNA Motif Discovery Explorer
# File: api/plumber_api.R
# Purpose: REST API layer via plumber — exposes backend pipeline as endpoints.
#          Deploy with: plumber::plumb("api/plumber_api.R")$run(port = 8080)
#          Designed for React frontend or cloud integration.
# =============================================================================

library(plumber)
library(jsonlite)

# Source modules (adjust paths relative to project root when deployed)
source("scripts/preprocessing.R")
source("scripts/motif_detection.R")
source("scripts/visualization.R")

#* @apiTitle RNA Motif Discovery API
#* @apiDescription RESTful API for RNA motif discovery and enrichment analysis.
#* @apiVersion 1.0.0

# ---- CORS (required for React/browser frontends) ----------------------------
#* @filter cors
function(req, res) {
  res$setHeader("Access-Control-Allow-Origin",  "*")
  res$setHeader("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
  res$setHeader("Access-Control-Allow-Headers",
                "Content-Type, Authorization")
  if (req$REQUEST_METHOD == "OPTIONS") {
    res$status <- 200
    return(list())
  }
  plumber::forward()
}

# -----------------------------------------------------------------------------
#* Health check endpoint
#* @get /health
#* @serializer json
function() {
  list(
    status    = "ok",
    service   = "RNA Motif Discovery API",
    version   = "1.0.0",
    timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
  )
}

# -----------------------------------------------------------------------------
#* Upload and preprocess a FASTA file
#* @post /preprocess
#* @param fasta_file:file FASTA file upload
#* @param min_length:int Minimum sequence length (default 10)
#* @serializer json
function(req, res, fasta_file, min_length = 10) {
  min_length <- as.integer(min_length)

  # Write temp file
  tmp <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp))

  tryCatch({
    writeBin(fasta_file$content, tmp)
    result <- preprocess_pipeline(tmp, min_length = min_length,
                                   verbose = FALSE)
    list(
      success          = TRUE,
      total_sequences  = result$stats$summary$total_sequences,
      mean_length      = result$stats$summary$mean_length,
      median_length    = result$stats$summary$median_length,
      min_length       = result$stats$summary$min_length,
      max_length       = result$stats$summary$max_length,
      mean_gc          = result$stats$summary$mean_gc,
      total_nucleotides = result$stats$summary$total_nucleotides,
      n_invalid        = nrow(result$invalid_report),
      n_duplicates     = nrow(result$duplicates),
      per_sequence     = result$stats$per_sequence
    )
  }, error = function(e) {
    res$status <- 400
    list(success = FALSE, error = conditionMessage(e))
  })
}

# -----------------------------------------------------------------------------
#* Discover motifs in preprocessed sequences (JSON body)
#* @post /discover
#* @param sequences:[str] Array of RNA sequence strings
#* @param control_sequences:[str] Optional control/background RNA sequences
#* @param k:int Motif length (default 6)
#* @param min_freq:double Min sequence frequency (default 0.05)
#* @param min_count:int Min occurrence count (default 2)
#* @param top_n:int Top N motifs to return (default 20)
#* @param background_model:str Background model: uniform, mononucleotide, control
#* @param p_adjust_method:str P-value correction method: bonferroni, fdr
#* @serializer json
function(req, res, sequences, control_sequences = NULL,
         k = 6, min_freq = 0.05, min_count = 2, top_n = 20,
         background_model = "uniform", p_adjust_method = "bonferroni") {
  k                <- as.integer(k)
  min_freq         <- as.double(min_freq)
  min_count        <- as.integer(min_count)
  top_n            <- as.integer(top_n)
  background_model <- tolower(background_model)
  p_adjust_method  <- tolower(p_adjust_method)

  if (length(sequences) == 0) {
    res$status <- 400
    return(list(success = FALSE, error = "No sequences provided."))
  }

  tryCatch({
    cleaned_sequences <- gsub("T", "U", toupper(sequences))
    rna_ss <- RNAStringSet(cleaned_sequences)

    control_ss <- NULL
    if (!is.null(control_sequences) && length(control_sequences) > 0) {
      cleaned_control <- gsub("T", "U", toupper(control_sequences))
      control_ss <- RNAStringSet(cleaned_control)
    }

    result <- motif_discovery_pipeline(
      rna_stringset         = rna_ss,
      k                     = k,
      min_freq              = min_freq,
      min_count             = min_count,
      top_n                 = top_n,
      background_model      = background_model,
      p_adjust_method       = p_adjust_method,
      control_rna_stringset = control_ss,
      verbose               = FALSE
    )

    list(
      success          = TRUE,
      total_kmers      = result$params$total_kmers_found,
      filtered         = result$params$kmers_after_filter,
      background_model = result$params$background_model,
      p_adjust_method  = result$params$p_adjust_method,
      control_dataset  = result$params$control_dataset,
      top_motifs       = as.list(result$top_motifs),
      motif_table      = as.list(result$motif_table)
    )
  }, error = function(e) {
    res$status <- 500
    list(success = FALSE, error = conditionMessage(e))
  })
}

# -----------------------------------------------------------------------------
#* Get top motif k-mer strings only (lightweight endpoint)
#* @post /top_kmers
#* @param sequences:[str] Array of RNA sequences
#* @param control_sequences:[str] Optional control/background RNA sequences
#* @param k:int Motif length
#* @param top_n:int Number of top motifs
#* @param background_model:str Background model: uniform, mononucleotide, control
#* @param p_adjust_method:str P-value correction: bonferroni, fdr
#* @serializer json
function(req, res, sequences, control_sequences = NULL,
         k = 6, top_n = 10,
         background_model = "uniform", p_adjust_method = "bonferroni") {
  k                <- as.integer(k)
  top_n            <- as.integer(top_n)
  background_model <- tolower(background_model)
  p_adjust_method  <- tolower(p_adjust_method)

  if (length(sequences) == 0) {
    res$status <- 400
    return(list(success = FALSE, error = "No sequences provided."))
  }

  tryCatch({
    cleaned_sequences <- gsub("T", "U", toupper(sequences))
    rna_ss <- RNAStringSet(cleaned_sequences)

    control_ss <- NULL
    if (!is.null(control_sequences) && length(control_sequences) > 0) {
      cleaned_control <- gsub("T", "U", toupper(control_sequences))
      control_ss <- RNAStringSet(cleaned_control)
    }

    result <- motif_discovery_pipeline(
      rna_stringset         = rna_ss,
      k                     = k,
      min_freq              = 0,
      min_count             = 1,
      top_n                 = top_n,
      background_model      = background_model,
      p_adjust_method       = p_adjust_method,
      control_rna_stringset = control_ss,
      verbose               = FALSE
    )

    list(
      success   = TRUE,
      top_kmers = head(result$motif_table$kmer, top_n)
    )
  }, error = function(e) {
    res$status <- 500
    list(success = FALSE, error = conditionMessage(e))
  })
}

# -----------------------------------------------------------------------------
#* Get supported motif lengths
#* @get /config
#* @serializer json
function() {
  list(
    min_k             = 3,
    max_k             = 15,
    default_k         = 6,
    default_min_freq  = 0.05,
    default_min_count = 2,
    default_top_n     = 20,
    supported_formats = c(".fasta", ".fa", ".txt"),
    background_models = c("uniform", "mononucleotide", "control"),
    p_adjust_methods  = c("bonferroni", "fdr")
  )
}
