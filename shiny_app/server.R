# =============================================================================
# RNA Motif Discovery Explorer
# File: shiny_app/server.R
# Purpose: Shiny server logic — reactive pipeline orchestration
# =============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(shinyjs)
library(dplyr)

# Source all modules
source("../scripts/preprocessing.R")
source("../scripts/motif_detection.R")
source("../scripts/visualization.R")
source("../scripts/report_generation.R")

server <- function(input, output, session) {

  # --------------------------------------------------------------------------
  # Reactive values store
  # --------------------------------------------------------------------------
  rv <- reactiveValues(
    preproc_result  = NULL,
    motif_result    = NULL,
    status_msg      = "Ready. Upload a FASTA file to begin.",
    pdf_path        = NULL
  )

  set_status <- function(msg, type = "info") {
    prefix <- switch(type,
      info    = "[INFO] ",
      success = "[OK]   ",
      warning = "[WARN] ",
      error   = "[ERR]  ",
      "[INFO] "
    )
    rv$status_msg <- paste0(prefix, msg)
  }

  # --------------------------------------------------------------------------
  # Step 1: Preprocessing
  # --------------------------------------------------------------------------
  observeEvent(input$run_preprocess, {
    req(input$fasta_file)
    set_status("Preprocessing sequences...", "info")

    withProgress(message = "Preprocessing FASTA...", value = 0, {
      tryCatch({
        incProgress(0.3, detail = "Reading file...")
        result <- preprocess_pipeline(
          filepath   = input$fasta_file$datapath,
          min_length = input$min_seq_length,
          verbose    = FALSE
        )
        incProgress(0.6, detail = "Computing statistics...")
        rv$preproc_result <- result
        rv$motif_result   <- NULL  # reset motif results on new upload

        incProgress(1.0, detail = "Done!")
        set_status(sprintf("Preprocessing complete: %d sequences ready.",
                           result$stats$summary$total_sequences), "success")

        # Update motif selector (empty until discovery runs)
        updateSelectInput(session, "logo_motif_select", choices = NULL)
      }, error = function(e) {
        set_status(paste("Preprocessing failed:", conditionMessage(e)), "error")
        rv$preproc_result <- NULL
      })
    })
  })

  # --------------------------------------------------------------------------
  # Step 2: Motif Discovery
  # --------------------------------------------------------------------------
  observeEvent(input$run_discovery, {
    req(rv$preproc_result)
    set_status("Running motif discovery...", "info")

    withProgress(message = "Discovering motifs...", value = 0, {
      tryCatch({
        incProgress(0.2, detail = "Counting k-mers...")
        result <- motif_discovery_pipeline(
          rna_stringset = rv$preproc_result$sequences,
          k             = input$motif_length,
          min_freq      = input$min_freq,
          min_count     = input$min_count,
          top_n         = input$top_n_motifs,
          n_pos_bins    = input$pos_bins,
          verbose       = FALSE
        )
        incProgress(0.8, detail = "Finalising ranks...")
        rv$motif_result <- result

        # Update logo selector with top motifs
        top_choices <- result$top_motifs$kmer
        updateSelectInput(session, "logo_motif_select",
                          choices = top_choices,
                          selected = top_choices[1])

        incProgress(1.0, detail = "Done!")
        set_status(
          sprintf("Discovery complete: %d motifs found, %d significant (k=%d).",
                  nrow(result$motif_table),
                  sum(result$motif_table$significant, na.rm = TRUE),
                  input$motif_length),
          "success"
        )
      }, error = function(e) {
        set_status(paste("Discovery failed:", conditionMessage(e)), "error")
      })
    })
  })

  # --------------------------------------------------------------------------
  # Value Boxes
  # --------------------------------------------------------------------------
  output$vbox_sequences <- renderValueBox({
    n <- if (!is.null(rv$preproc_result))
      rv$preproc_result$stats$summary$total_sequences else "—"
    valueBox(n, "Sequences", icon = icon("dna"),
             color = "blue")
  })

  output$vbox_mean_length <- renderValueBox({
    v <- if (!is.null(rv$preproc_result))
      round(rv$preproc_result$stats$summary$mean_length, 1) else "—"
    valueBox(v, "Mean Length (nt)", icon = icon("ruler-horizontal"),
             color = "teal")
  })

  output$vbox_motifs_found <- renderValueBox({
    n <- if (!is.null(rv$motif_result))
      nrow(rv$motif_result$motif_table) else "—"
    valueBox(n, "Unique Motifs", icon = icon("search"),
             color = "green")
  })

  output$vbox_significant <- renderValueBox({
    n <- if (!is.null(rv$motif_result))
      sum(rv$motif_result$motif_table$significant, na.rm = TRUE) else "—"
    valueBox(n, "Significant", icon = icon("star"),
             color = "red")
  })

  # --------------------------------------------------------------------------
  # Status output
  # --------------------------------------------------------------------------
  output$status_message <- renderText({ rv$status_msg })

  output$preprocess_log <- renderText({
    req(rv$preproc_result)
    paste(rv$preproc_result$log, collapse = "\n")
  })

  # --------------------------------------------------------------------------
  # Overview Plots
  # --------------------------------------------------------------------------
  output$plot_length_dist <- renderPlotly({
    req(rv$preproc_result)
    plot_sequence_length_distribution(rv$preproc_result$stats$per_sequence)
  })

  output$plot_nt_freq <- renderPlotly({
    req(rv$preproc_result)
    plot_nucleotide_frequencies(rv$preproc_result$stats$per_sequence)
  })

  output$plot_gc <- renderPlotly({
    req(rv$preproc_result)
    plot_gc_content(rv$preproc_result$stats$per_sequence)
  })

  output$plot_summary_bar <- renderPlotly({
    req(rv$preproc_result, rv$motif_result)
    plot_dashboard_summary(rv$preproc_result$stats$summary,
                            rv$motif_result$params)
  })

  # --------------------------------------------------------------------------
  # Motif Logos Tab
  # --------------------------------------------------------------------------
  output$plot_logo_selected <- renderPlot({
    req(rv$motif_result, input$logo_motif_select)
    km <- input$logo_motif_select
    # Get top-k sequences containing this kmer
    seqs_char <- as.character(rv$preproc_result$sequences)
    hit_seqs  <- seqs_char[grepl(km, seqs_char, fixed = TRUE)]
    # Extract the aligned kmer context
    k         <- nchar(km)
    kmer_instances <- unlist(lapply(hit_seqs, function(s) {
      positions <- find_overlapping_positions(s, km)
      if (length(positions) == 0) return(NULL)
      sapply(as.integer(positions), function(p) substr(s, p, p + k - 1))
    }))
    kmer_instances <- unique(kmer_instances[!is.na(kmer_instances)])
    if (length(kmer_instances) == 0) kmer_instances <- c(km)
    plot_sequence_logo(kmer_instances, title = paste("Logo:", km))
  }, bg = "#0D1117")

  output$plot_logos_grid <- renderPlot({
    req(rv$motif_result)
    top_kmers <- head(rv$motif_result$top_motifs$kmer, 6)
    seqs_char <- as.character(rv$preproc_result$sequences)

    plots <- lapply(top_kmers, function(km) {
      k    <- nchar(km)
      hits <- unlist(lapply(seqs_char, function(s) {
        p <- find_overlapping_positions(s, km)
        if (length(p) == 0) return(NULL)
        sapply(as.integer(p), function(pos) substr(s, pos, pos + k - 1))
      }))
      hits <- unique(hits[!is.na(hits)])
      if (length(hits) == 0) hits <- c(km)
      plot_sequence_logo(hits, title = km)
    })
    plots <- Filter(Negate(is.null), plots)

    if (length(plots) > 0) {
      n_plots <- length(plots)
      n_cols  <- min(2, n_plots)
      n_rows  <- ceiling(n_plots / n_cols)

      # Combine with gridExtra
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        gridExtra::grid.arrange(grobs = plots,
                                 ncol = n_cols,
                                 nrow = n_rows)
      } else {
        print(plots[[1]])
      }
    }
  }, bg = "#0D1117")

  # --------------------------------------------------------------------------
  # Frequency Analysis Tab
  # --------------------------------------------------------------------------
  output$plot_freq_bar <- renderPlotly({
    req(rv$motif_result)
    plot_motif_frequency_bar(rv$motif_result$top_motifs,
                              n = input$n_top_display)
  })

  output$plot_volcano <- renderPlotly({
    req(rv$motif_result)
    plot_enrichment_volcano(rv$motif_result$motif_table,
                             top_n_labels = 10)
  })

	  output$table_motifs <- renderDT({
	    req(rv$motif_result)
    df <- rv$motif_result$motif_table %>%
      select(kmer, total_count, seq_count, seq_frequency,
             enrichment_score, p_value, p_adjusted,
             significant, combined_rank) %>%
      mutate(
        seq_frequency   = round(seq_frequency, 4),
        enrichment_score = round(enrichment_score, 3),
        p_value         = formatC(p_value, format = "e", digits = 2),
        p_adjusted      = formatC(p_adjusted, format = "e", digits = 2)
      )

    datatable(
      df,
      filter    = "top",
      rownames  = FALSE,
      options   = list(
        pageLength = 20,
        scrollX    = TRUE,
        dom        = "Bfrtip",
        columnDefs = list(
          list(className = "dt-center",
               targets   = seq_len(ncol(df)) - 1)
        )
      ),
      colnames = c("k-mer", "Total Count", "Seq Count",
                   "Seq Freq", "Enrichment", "p-value",
                   "adj. p", "Significant", "Rank"),
      class = "cell-border stripe"
    ) %>%
      formatStyle("significant",
                   backgroundColor = styleEqual(
                     c(TRUE, FALSE),
                     c("#0E2B1A", "#2B0E0E")
                   ))
  })

	  # --------------------------------------------------------------------------
	  # Enrichment Heatmap Tab
	  # --------------------------------------------------------------------------
	  current_positional_dist <- reactive({
	    req(rv$motif_result, rv$preproc_result)
	    n_motifs <- max(input$n_top_display, 8)
	    top_kmers <- unique(c(
	      head(rv$motif_result$top_motifs$kmer, n_motifs),
	      head(rv$motif_result$motif_table$kmer, n_motifs)
	    ))
	    pos_dist <- compute_positional_distribution(
	      rv$preproc_result$sequences,
	      top_kmers,
	      n_bins = input$pos_bins
	    )
	    rv$motif_result$positional_dist <- pos_dist
	    pos_dist
	  })

	  output$plot_pos_heatmap <- renderPlotly({
	    req(rv$motif_result)
	    top_kmers <- head(rv$motif_result$top_motifs$kmer, input$n_top_display)
	    plot_positional_heatmap(current_positional_dist(),
	                             top_kmers  = top_kmers,
	                             n_bins     = input$pos_bins)
	  })

	  output$plot_density <- renderPlotly({
	    req(rv$motif_result)
	    plot_motif_density(current_positional_dist(), top_n = 8)
	  })

  # --------------------------------------------------------------------------
  # Download Handlers — CSV
  # --------------------------------------------------------------------------
  output$dl_motif_csv <- downloadHandler(
    filename = function() {
      paste0("motif_table_k", input$motif_length, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content  = function(file) {
      req(rv$motif_result)
      export_motif_table_csv(rv$motif_result$motif_table, file)
    }
  )

  output$dl_stats_csv <- downloadHandler(
    filename = function() {
      paste0("sequence_stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content  = function(file) {
      req(rv$preproc_result)
      export_stats_csv(rv$preproc_result$stats$per_sequence, file)
    }
  )

  output$dl_pos_csv <- downloadHandler(
    filename = function() {
      paste0("positional_dist_k", input$motif_length, "_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content  = function(file) {
      req(rv$motif_result)
      export_positional_csv(rv$motif_result$positional_dist, file)
    }
  )

  # --------------------------------------------------------------------------
  # Download Handlers — Plots
  # --------------------------------------------------------------------------
  output$dl_plot_freq <- downloadHandler(
    filename = function() paste0("freq_bar_k", input$motif_length, ".png"),
    content  = function(file) {
      req(rv$motif_result)
      # Convert ggplotly back to ggplot for saving
      p <- ggplot(head(rv$motif_result$top_motifs, input$n_top_display),
                  aes(x = total_count,
                      y = reorder(kmer, total_count),
                      fill = enrichment_score)) +
        geom_col(alpha = 0.9) +
        scale_fill_gradient2(low = "#388BFD", mid = "#F0C27F",
                              high = "#FF7B72", midpoint = 0) +
        labs(title = "Top Motifs by Frequency", x = "Count", y = "k-mer") +
        THEME_RNA()
      ggsave(file, plot = p, width = 10, height = 7, dpi = 150,
             bg = "#0D1117")
    }
  )

  output$dl_plot_volcano <- downloadHandler(
    filename = function() paste0("volcano_k", input$motif_length, ".png"),
    content  = function(file) {
      req(rv$motif_result)
      df  <- rv$motif_result$motif_table
      p   <- ggplot(df, aes(x = enrichment_score, y = neg_log10_p,
                              color = significant)) +
        geom_point(alpha = 0.7, size = 2) +
        scale_color_manual(values = c("TRUE" = "#FF7B72",
                                       "FALSE" = "#388BFD")) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   color = "#8B949E") +
        labs(title = "Enrichment Volcano",
             x = "log2 Enrichment", y = "-log10(p)") +
        THEME_RNA()
      ggsave(file, plot = p, width = 8, height = 6, dpi = 150,
             bg = "#0D1117")
    }
  )

  output$dl_plot_logos <- downloadHandler(
    filename = function() paste0("logos_grid_k", input$motif_length, ".png"),
    content  = function(file) {
      req(rv$motif_result, rv$preproc_result)
      top_kmers <- head(rv$motif_result$top_motifs$kmer, 6)
      seqs_char <- as.character(rv$preproc_result$sequences)
      plots <- lapply(top_kmers, function(km) {
        k    <- nchar(km)
        hits <- unlist(lapply(seqs_char, function(s) {
          p <- find_overlapping_positions(s, km)
          if (length(p) == 0) return(NULL)
          sapply(as.integer(p), function(pos) substr(s, pos, pos + k - 1))
        }))
        hits <- unique(hits[!is.na(hits)])
        if (length(hits) == 0) hits <- c(km)
        plot_sequence_logo(hits, title = km)
      })
      plots <- Filter(Negate(is.null), plots)
      if (requireNamespace("gridExtra", quietly = TRUE)) {
        png(file, width = 1400, height = 900, bg = "#0D1117")
        gridExtra::grid.arrange(grobs = plots, ncol = 2)
        dev.off()
      }
    }
  )

  # --------------------------------------------------------------------------
  # PDF Report Generation
  # --------------------------------------------------------------------------
  observeEvent(input$generate_pdf, {
    req(rv$motif_result, rv$preproc_result)
    set_status("Generating PDF report...", "info")

    withProgress(message = "Rendering PDF...", value = 0.5, {
      tryCatch({
        out_dir  <- tempdir()
        params   <- c(rv$motif_result$params,
                       list(filepath   = input$fasta_file$name,
                            min_freq   = input$min_freq,
                            min_count  = input$min_count))
        pdf_path <- generate_rmd_report(
          summary_stats = rv$preproc_result$stats$summary,
          motif_table   = rv$motif_result$motif_table,
          top_motifs    = rv$motif_result$top_motifs,
          per_seq_stats = rv$preproc_result$stats$per_sequence,
          params_list   = params,
          out_dir       = out_dir
        )
        rv$pdf_path <- pdf_path
        incProgress(1.0)
        if (!is.null(pdf_path) && file.exists(pdf_path)) {
          set_status("PDF report ready for download.", "success")
        } else {
          set_status("PDF generation failed. Check rmarkdown/tinytex.", "warning")
        }
      }, error = function(e) {
        set_status(paste("PDF error:", conditionMessage(e)), "error")
      })
    })
  })

  output$pdf_download_ui <- renderUI({
    req(rv$pdf_path)
    if (!is.null(rv$pdf_path) && file.exists(rv$pdf_path)) {
      downloadButton("dl_pdf", "Download PDF Report",
                     class = "btn-danger btn-block")
    }
  })

  output$dl_pdf <- downloadHandler(
    filename = function() {
      paste0("rna_motif_report_",
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content  = function(file) {
      req(rv$pdf_path)
      file.copy(rv$pdf_path, file)
    }
  )

  # --------------------------------------------------------------------------
  # Parameters Summary
  # --------------------------------------------------------------------------
  output$params_summary <- renderText({
    lines <- c(
      "=== Analysis Parameters ===",
      sprintf("Input file     : %s",
              if (!is.null(input$fasta_file)) input$fasta_file$name
              else "Not loaded"),
      sprintf("Min seq length : %d nt",   input$min_seq_length),
      sprintf("Motif length   : k = %d",  input$motif_length),
      sprintf("Min seq freq   : %.2f",    input$min_freq),
      sprintf("Min count      : %d",      input$min_count),
      sprintf("Top N motifs   : %d",      input$top_n_motifs),
      sprintf("Position bins  : %d",      input$pos_bins),
      ""
    )
    if (!is.null(rv$preproc_result)) {
      s <- rv$preproc_result$stats$summary
      lines <- c(lines,
        "=== Dataset Statistics ===",
        sprintf("Total sequences    : %d", s$total_sequences),
        sprintf("Mean length        : %.1f nt", s$mean_length),
        sprintf("Median length      : %.1f nt", s$median_length),
        sprintf("GC content (mean)  : %.1f%%", s$mean_gc * 100),
        sprintf("Total nucleotides  : %s",
                format(s$total_nucleotides, big.mark = ","))
      )
    }
    if (!is.null(rv$motif_result)) {
      p <- rv$motif_result$params
      lines <- c(lines,
        "",
        "=== Motif Discovery ===",
        sprintf("Unique k-mers found   : %d", p$total_kmers_found),
        sprintf("k-mers after filter   : %d", p$kmers_after_filter),
        sprintf("Significant motifs    : %d",
                sum(rv$motif_result$motif_table$significant, na.rm = TRUE))
      )
    }
    paste(lines, collapse = "\n")
  })
}
