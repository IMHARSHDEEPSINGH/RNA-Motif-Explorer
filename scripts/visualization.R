# =============================================================================
# RNA Motif Discovery Explorer
# Script: visualization.R
# Purpose: All plot generation functions — logos, frequency, heatmaps, etc.
# =============================================================================

library(ggplot2)
library(ggseqlogo)
library(universalmotif)
library(plotly)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Global RNA color scheme (consistent across all plots)
RNA_COLORS <- c(
  A = "#E74C3C",   # red
  U = "#3498DB",   # blue
  G = "#2ECC71",   # green
  C = "#F39C12"    # orange
)

THEME_RNA <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.background  = element_rect(fill = "#0D1117", color = NA),
      panel.background = element_rect(fill = "#161B22", color = NA),
      panel.grid.major = element_line(color = "#21262D"),
      panel.grid.minor = element_blank(),
      text             = element_text(color = "#E6EDF3"),
      axis.text        = element_text(color = "#8B949E"),
      axis.title       = element_text(color = "#C9D1D9"),
      plot.title       = element_text(color = "#F0F6FC", face = "bold",
                                       size = 15),
      plot.subtitle    = element_text(color = "#8B949E", size = 11),
      legend.background = element_rect(fill = "#161B22", color = NA),
      legend.text       = element_text(color = "#C9D1D9"),
      legend.title      = element_text(color = "#E6EDF3"),
      strip.text        = element_text(color = "#E6EDF3", face = "bold")
    )
}

plotly_empty_rna <- function(title = "No data available") {
  plot_ly(x = numeric(0), y = numeric(0), type = "scatter", mode = "markers") %>%
    layout(
      title = list(text = title, font = list(color = "#8B949E")),
      paper_bgcolor = "#0D1117",
      plot_bgcolor = "#161B22",
      font = list(color = "#E6EDF3"),
      xaxis = list(visible = FALSE, zeroline = FALSE),
      yaxis = list(visible = FALSE, zeroline = FALSE)
    )
}

# -----------------------------------------------------------------------------
# plot_sequence_length_distribution
# Histogram of sequence lengths with density overlay.
# -----------------------------------------------------------------------------
plot_sequence_length_distribution <- function(per_seq_stats) {
  p <- ggplot(per_seq_stats, aes(x = length)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins   = 30,
                   fill   = "#388BFD",
                   color  = "#58A6FF",
                   alpha  = 0.8) +
    geom_density(color = "#FF7B72", linewidth = 1.2, adjust = 1.5) +
    labs(
      title    = "Sequence Length Distribution",
      subtitle = sprintf("n = %d sequences", nrow(per_seq_stats)),
      x        = "Sequence Length (nt)",
      y        = "Density"
    ) +
    THEME_RNA()

  ggplotly(p, tooltip = c("x", "y")) %>%
    layout(paper_bgcolor = "#0D1117", plot_bgcolor = "#161B22",
           font = list(color = "#E6EDF3"))
}

# -----------------------------------------------------------------------------
# plot_nucleotide_frequencies
# Grouped bar chart showing mean A/U/G/C frequencies.
# -----------------------------------------------------------------------------
plot_nucleotide_frequencies <- function(per_seq_stats) {
  freq_long <- per_seq_stats %>%
    select(name, freq_A, freq_U, freq_G, freq_C) %>%
    pivot_longer(cols = -name,
                 names_to  = "nucleotide",
                 values_to = "frequency") %>%
    mutate(nucleotide = gsub("freq_", "", nucleotide))

  mean_freqs <- freq_long %>%
    group_by(nucleotide) %>%
    summarise(mean_freq = mean(frequency), sd_freq = sd(frequency),
              .groups = "drop")

  p <- ggplot(mean_freqs, aes(x = nucleotide, y = mean_freq,
                               fill = nucleotide)) +
    geom_col(width = 0.6, alpha = 0.9) +
    geom_errorbar(aes(ymin = pmax(mean_freq - sd_freq, 0),
                      ymax = mean_freq + sd_freq),
                  width = 0.2, color = "#E6EDF3") +
    scale_fill_manual(values = RNA_COLORS, guide = "none") +
    labs(
      title    = "Mean Nucleotide Frequencies",
      subtitle = "Error bars = ±1 SD across sequences",
      x        = "Nucleotide",
      y        = "Mean Frequency"
    ) +
    THEME_RNA()

  ggplotly(p, tooltip = c("x", "y")) %>%
    layout(paper_bgcolor = "#0D1117", plot_bgcolor = "#161B22",
           font = list(color = "#E6EDF3"))
}

# -----------------------------------------------------------------------------
# plot_gc_content
# Violin + box plot of GC content distribution.
# -----------------------------------------------------------------------------
plot_gc_content <- function(per_seq_stats) {
  mean_gc <- mean(per_seq_stats$gc_content)
  p <- ggplot(per_seq_stats, aes(x = "Sequences", y = gc_content))

  if (nrow(per_seq_stats) < 2) {
    p <- p +
      geom_boxplot(width = 0.1, fill = "#58A6FF",
                   color = "#E6EDF3", outlier.color = "#FF7B72") +
      geom_jitter(width = 0.1, height = 0, size = 3,
                  color = "#58A6FF", alpha = 0.9)
  } else {
    p <- p +
      geom_violin(fill = "#388BFD", alpha = 0.6, color = NA) +
      geom_boxplot(width = 0.1, fill = "#58A6FF",
                   color = "#E6EDF3", outlier.color = "#FF7B72")
  }

  p <- p +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title    = "GC Content Distribution",
      subtitle = sprintf("Mean GC = %.1f%%", mean_gc * 100),
      x        = NULL,
      y        = "GC Content"
    ) +
    THEME_RNA()

  ggplotly(p, tooltip = c("y")) %>%
    layout(paper_bgcolor = "#0D1117", plot_bgcolor = "#161B22",
           font = list(color = "#E6EDF3"))
}

# -----------------------------------------------------------------------------
# plot_motif_frequency_bar
# Ranked bar chart of top N motifs by total count.
# -----------------------------------------------------------------------------
plot_motif_frequency_bar <- function(top_motifs, n = 20) {
  df <- head(top_motifs, n) %>%
    mutate(kmer = factor(kmer, levels = rev(kmer)))

  p <- ggplot(df, aes(x = total_count, y = kmer,
                       fill = enrichment_score,
                       text = paste0("k-mer: ", kmer,
                                     "<br>Count: ", total_count,
                                     "<br>Seq freq: ",
                                     scales::percent(seq_frequency, 0.1),
                                     "<br>Enrichment: ",
                                     round(enrichment_score, 2)))) +
    geom_col(alpha = 0.9) +
    scale_fill_gradient2(low = "#388BFD", mid = "#F0C27F",
                          high = "#FF7B72", midpoint = 0,
                          name = "log2\nEnrichment") +
    labs(
      title    = sprintf("Top %d Motifs by Frequency", n),
      subtitle = "Colored by enrichment score",
      x        = "Total Occurrences",
      y        = "k-mer"
    ) +
    THEME_RNA()

  ggplotly(p, tooltip = "text") %>%
    layout(paper_bgcolor = "#0D1117", plot_bgcolor = "#161B22",
           font = list(color = "#E6EDF3"))
}

# -----------------------------------------------------------------------------
# plot_enrichment_volcano
# Volcano plot: enrichment score vs. -log10(p-value).
# -----------------------------------------------------------------------------
plot_enrichment_volcano <- function(motif_table, top_n_labels = 10) {
  df <- motif_table %>%
    mutate(
      significant = p_adjusted < 0.05,
      label       = ifelse(combined_rank <= top_n_labels, kmer, NA_character_)
    )

  p <- ggplot(df, aes(x = enrichment_score, y = neg_log10_p,
                       color = significant,
                       text = paste0("k-mer: ", kmer,
                                     "<br>Enrichment: ",
                                     round(enrichment_score, 3),
                                     "<br>-log10(p): ",
                                     round(neg_log10_p, 3),
                                     "<br>Adj p: ",
                                     formatC(p_adjusted,
                                             format = "e", digits = 2)))) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("TRUE" = "#FF7B72", "FALSE" = "#388BFD"),
                        name = "Significant\n(adj. p<0.05)") +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "#8B949E", linewidth = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "#F0C27F", linewidth = 0.5) +
    labs(
      title    = "Enrichment Volcano Plot",
      subtitle = "Red = statistically significant (Bonferroni corrected)",
      x        = "log2 Enrichment Score",
      y        = "-log10(p-value)"
    ) +
    THEME_RNA()

  ggplotly(p, tooltip = "text") %>%
    layout(paper_bgcolor = "#0D1117", plot_bgcolor = "#161B22",
           font = list(color = "#E6EDF3"))
}

# -----------------------------------------------------------------------------
# plot_positional_heatmap
# Heatmap of motif density across sequence position bins.
# -----------------------------------------------------------------------------
plot_positional_heatmap <- function(pos_dist, top_kmers = NULL,
                                     n_bins = 20) {
  if (is.null(pos_dist) || nrow(pos_dist) == 0) {
    return(plotly_empty_rna("No positional data available"))
  }

  if (!is.null(top_kmers)) {
    pos_dist <- pos_dist %>% filter(kmer %in% top_kmers)
  }

  if (nrow(pos_dist) == 0) {
    return(plotly_empty_rna("No positional data for selected motifs"))
  }

  heat_data <- pos_dist %>%
    group_by(kmer, position_bin) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(kmer) %>%
    mutate(density = count / sum(count)) %>%
    ungroup()

  # Build full grid (fill zeros for missing bins)
  all_kmers <- unique(heat_data$kmer)
  all_bins  <- seq_len(n_bins)
  grid      <- expand.grid(kmer = all_kmers, position_bin = all_bins,
                            stringsAsFactors = FALSE)
  heat_full <- grid %>%
    left_join(heat_data, by = c("kmer", "position_bin")) %>%
    replace_na(list(density = 0, count = 0))

  mat       <- heat_full %>%
    pivot_wider(id_cols = kmer, names_from = position_bin,
                values_from = density, values_fill = 0)
  kmer_labels <- mat$kmer
  mat_vals    <- as.matrix(mat[, -1])

  bin_labels <- paste0(round(seq(0, 1 - 1/n_bins, by = 1/n_bins) * 100), "%")

  plot_ly(
    z          = mat_vals,
    x          = bin_labels,
    y          = kmer_labels,
    type       = "heatmap",
    colorscale = list(c(0, "#161B22"), c(0.5, "#388BFD"), c(1, "#FF7B72")),
    hovertemplate = "k-mer: %{y}<br>Position: %{x}<br>Density: %{z:.3f}<extra></extra>"
  ) %>%
    layout(
      title      = list(text = "Motif Positional Enrichment Heatmap",
                         font = list(color = "#F0F6FC")),
      xaxis      = list(title = "Relative Position in Sequence",
                         tickfont = list(color = "#8B949E"),
                         titlefont = list(color = "#C9D1D9")),
      yaxis      = list(title = "k-mer",
                         tickfont = list(color = "#8B949E"),
                         titlefont = list(color = "#C9D1D9")),
      paper_bgcolor = "#0D1117",
      plot_bgcolor  = "#161B22",
      font          = list(color = "#E6EDF3")
    )
}

# -----------------------------------------------------------------------------
# plot_motif_density
# Density curves of positional distribution for top motifs.
# -----------------------------------------------------------------------------
plot_motif_density <- function(pos_dist, top_n = 8) {
  if (is.null(pos_dist) || nrow(pos_dist) == 0) {
    return(plotly_empty_rna("No positional data available"))
  }

  top_kmers <- pos_dist %>%
    count(kmer) %>%
    top_n(top_n, n) %>%
    pull(kmer)

  df <- pos_dist %>% filter(kmer %in% top_kmers)

  if (nrow(df) == 0) {
    return(plotly_empty_rna("No positional data for selected motifs"))
  }

  p <- ggplot(df, aes(x = relative_position, color = kmer, fill = kmer)) +
    geom_density(alpha = 0.15, linewidth = 1) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
      title    = sprintf("Motif Positional Density (Top %d)", top_n),
      subtitle = "Relative position within sequence",
      x        = "Relative Position",
      y        = "Density",
      color    = "k-mer",
      fill     = "k-mer"
    ) +
    THEME_RNA()

  ggplotly(p) %>%
    layout(paper_bgcolor = "#0D1117", plot_bgcolor = "#161B22",
           font = list(color = "#E6EDF3"))
}

# -----------------------------------------------------------------------------
# plot_sequence_logo
# Generates a sequence logo for a given character vector of aligned k-mers.
# Returns a ggplot object.
# -----------------------------------------------------------------------------
plot_sequence_logo <- function(kmer_vec, title = "Sequence Logo") {
  if (length(kmer_vec) == 0) return(NULL)
  k <- nchar(kmer_vec[1])
  if (k == 0) return(NULL)

  tryCatch({
    p <- ggseqlogo(kmer_vec,
                   method = "prob",
                   col_scheme = make_col_scheme(
                     chars = c("A", "U", "G", "C"),
                     cols  = unname(RNA_COLORS)
                   )) +
      ggtitle(title) +
      THEME_RNA() +
      theme(
        axis.text.x  = element_text(size = 10),
        axis.text.y  = element_text(size = 9),
        plot.title   = element_text(size = 13)
      )
    p
  }, error = function(e) {
    ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Logo unavailable:", conditionMessage(e)),
               color = "#8B949E") +
      THEME_RNA()
  })
}

# -----------------------------------------------------------------------------
# plot_universalmotif_logo
# Uses universalmotif to render a PWM-based logo. Returns a ggplot.
# -----------------------------------------------------------------------------
plot_universalmotif_logo <- function(motif_obj, title = NULL) {
  if (is.null(motif_obj)) return(NULL)
  tryCatch({
    view_motifs(motif_obj, show.positions = TRUE) +
      THEME_RNA()
  }, error = function(e) NULL)
}

# -----------------------------------------------------------------------------
# plot_dashboard_summary
# Bar chart summarizing key dataset metrics for the Overview tab.
# -----------------------------------------------------------------------------
plot_dashboard_summary <- function(summary_stats, motif_params) {
  significant_motifs <- if (is.null(motif_params$significant_motifs)) {
    0
  } else {
    motif_params$significant_motifs
  }
  unique_motifs <- if (is.null(motif_params$kmers_after_filter)) {
    motif_params$total_kmers_found
  } else {
    motif_params$kmers_after_filter
  }

  metrics <- tibble(
    label = c("Total Sequences", "Mean Length (nt)",
               "Unique Motifs Found", "Significant Motifs"),
    value = c(
      summary_stats$total_sequences,
      round(summary_stats$mean_length, 1),
      unique_motifs,
      significant_motifs
    ),
    display_value = scales::comma(value)
  ) %>%
    mutate(
      label = factor(label, levels = rev(label)),
      plot_value = log10(value + 1)
    )

  p <- ggplot(metrics, aes(x = label, y = plot_value,
	                         fill = label,
	                         text = paste0(label, ": ", display_value))) +
    geom_col(width = 0.6, alpha = 0.9, show.legend = FALSE) +
    geom_text(aes(label = display_value),
              hjust = -0.15, color = "#E6EDF3", size = 4) +
    coord_flip() +
    scale_fill_manual(values = c(
      "Total Sequences" = "#2ECC71",
      "Mean Length (nt)" = "#388BFD",
      "Unique Motifs Found" = "#F0C27F",
      "Significant Motifs" = "#FF7B72"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "Dataset Summary", x = NULL, y = "Value (log10 scale)") +
    THEME_RNA()

  ggplotly(p, tooltip = "text") %>%
    layout(paper_bgcolor = "#0D1117", plot_bgcolor = "#161B22",
           font = list(color = "#E6EDF3"))
}
