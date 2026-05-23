# =============================================================================
# RNA Motif Discovery Explorer
# File: shiny_app/ui.R
# Purpose: Shiny dashboard user interface
# =============================================================================

library(shiny)
library(shinydashboard)
library(DT)
library(plotly)
library(shinyjs)

# Custom CSS for dark bioinformatics aesthetic
custom_css <- "
  /* Global background */
  body, .wrapper { background-color: #0D1117 !important; }
  .content-wrapper, .right-side { background-color: #0D1117 !important; }

  /* Sidebar */
  .main-sidebar, .left-side { background-color: #161B22 !important; }
  .sidebar-menu > li > a { color: #8B949E !important; }
  .sidebar-menu > li.active > a,
  .sidebar-menu > li > a:hover { color: #58A6FF !important;
    background-color: #21262D !important; }
  .sidebar-menu > li > a > .fa { color: #58A6FF !important; }

	  /* Header */
	  .main-header .navbar { background-color: #161B22 !important; }
	  .main-header .logo { background-color: #0D1117 !important;
	    color: #58A6FF !important; font-weight: 700; font-size: 14px; }
	  .main-header .logo:hover { background-color: #161B22 !important; }
	  .app-logo-icon { margin-right: 6px; color: #58A6FF; }

  /* Boxes */
  .box { background: #161B22 !important; border-top-color: #388BFD !important;
    border-radius: 8px; }
  .box-body { background: #161B22 !important; }
  .box-header { color: #E6EDF3 !important; }
  .box-title { color: #58A6FF !important; font-size: 14px; font-weight: 600; }
  .content { background-color: #0D1117 !important; }
  .tab-pane .row { margin-left: 0; margin-right: 0; }
  .tab-pane .box { margin-bottom: 20px; }

  /* Plotly sizing */
  .plotly.html-widget,
  .js-plotly-plot,
  .plot-container,
  .svg-container {
    width: 100% !important;
  }

  /* Value boxes */
  .small-box { border-radius: 8px !important; }
  .small-box:hover { transform: translateY(-2px);
    box-shadow: 0 4px 15px rgba(56,139,253,0.3) !important; }
  .small-box .inner h3 { font-size: 28px; font-weight: 700; }

  /* Tabs */
  .nav-tabs-custom { background: #161B22 !important; }
  .nav-tabs-custom > .nav-tabs > li.active > a { color: #58A6FF !important;
    background: #21262D !important; border-top-color: #388BFD !important; }
  .nav-tabs-custom > .nav-tabs > li > a { color: #8B949E !important; }
  .nav-tabs-custom > .tab-content { background: #161B22 !important; }

  /* Inputs */
  .form-control { background-color: #21262D !important; color: #E6EDF3 !important;
    border-color: #30363D !important; }
  .form-control:focus { border-color: #388BFD !important;
    box-shadow: 0 0 0 3px rgba(56,139,253,0.15) !important; }
  label { color: #C9D1D9 !important; font-size: 12px; font-weight: 500; }

  /* Sliders */
  .irs-bar, .irs-bar-edge { background: #388BFD !important;
    border-color: #388BFD !important; }
  .irs-single { background: #388BFD !important; }
  .irs-line { background: #30363D !important; }
  .irs-grid-text { color: #8B949E !important; }
  .irs-min, .irs-max { color: #8B949E !important; }

  /* Buttons */
  .btn-primary { background-color: #388BFD !important;
    border-color: #58A6FF !important; color: #0D1117 !important;
    font-weight: 600; }
  .btn-primary:hover { background-color: #58A6FF !important; }
  .btn-success { background-color: #2ECC71 !important;
    border-color: #2ECC71 !important; color: #0D1117 !important; }
  .btn-warning { background-color: #F0C27F !important;
    border-color: #F0C27F !important; color: #0D1117 !important; }
  .btn-danger  { background-color: #FF7B72 !important;
    border-color: #FF7B72 !important; }

  /* DataTables */
  .dataTables_wrapper { color: #E6EDF3 !important; }
  table.dataTable { background: #0D1117 !important; color: #E6EDF3 !important; }
  table.dataTable thead th { background: #21262D !important;
    color: #58A6FF !important; border-bottom-color: #30363D !important; }
  table.dataTable tbody tr { background: #161B22 !important; }
  table.dataTable tbody tr:hover { background: #21262D !important; }
  table.dataTable tbody td { border-color: #21262D !important; }
  .dataTables_filter input { background: #21262D !important;
    color: #E6EDF3 !important; border-color: #30363D !important; }
  .dataTables_length select { background: #21262D !important;
    color: #E6EDF3 !important; border-color: #30363D !important; }
  .paginate_button { color: #8B949E !important; }
  .paginate_button.current { background: #388BFD !important;
    color: #0D1117 !important; border-radius: 4px; }

  /* Upload area */
  .shiny-file-upload { border: 2px dashed #388BFD !important;
    background: #161B22 !important; border-radius: 8px;
    color: #8B949E !important; padding: 10px; }

  /* Progress bar */
  .progress-bar { background-color: #388BFD !important; }

  /* Alert boxes */
  .alert-info { background-color: #1C2B3A !important;
    border-color: #388BFD !important; color: #58A6FF !important; }
  .alert-warning { background-color: #2B1F0A !important;
    border-color: #F0C27F !important; color: #F0C27F !important; }
  .alert-danger { background-color: #2B0E0E !important;
    border-color: #FF7B72 !important; color: #FF7B72 !important; }
  .alert-success { background-color: #0E2B1A !important;
    border-color: #2ECC71 !important; color: #2ECC71 !important; }

  /* Status text */
  #status_message { font-family: 'Courier New', monospace;
    font-size: 12px; color: #2ECC71; padding: 8px;
    background: #0D1117; border-radius: 4px; }
  #preprocess_log { font-family: 'Courier New', monospace;
    font-size: 11px; color: #8B949E; white-space: pre-wrap;
    max-height: 200px; overflow-y: auto; }

  /* Divider */
  hr { border-color: #21262D; }

  /* Sidebar section headers */
  .sidebar-section-header { color: #58A6FF; font-size: 10px;
    font-weight: 700; letter-spacing: 1.5px; text-transform: uppercase;
    margin: 12px 8px 4px; padding-bottom: 4px;
    border-bottom: 1px solid #21262D; }
"

# =============================================================================
# UI Definition
# =============================================================================
ui <- dashboardPage(
  skin = "blue",

  # ---- Header ----
	  dashboardHeader(
	    title = tags$span(
	      icon("dna", class = "app-logo-icon"),
	      "RNA Motif Explorer"
	    ),
    titleWidth = 240
  ),

  # ---- Sidebar ----
  dashboardSidebar(
    width = 240,
    useShinyjs(),
    tags$style(custom_css),

    tags$div(class = "sidebar-section-header", "Input"),

    fileInput(
      inputId  = "fasta_file",
      label    = "Upload FASTA File",
      accept   = c(".fasta", ".fa", ".fna", ".txt"),
      placeholder = "Drop .fasta / .fa / .fna / .txt"
    ),

    numericInput(
      inputId = "min_seq_length",
      label   = "Min Sequence Length (nt)",
      value   = 10, min = 1, max = 500, step = 1
    ),

    actionButton("run_preprocess", "Preprocess Sequences",
                 icon  = icon("filter"),
                 class = "btn-primary btn-block",
                 style = "margin-bottom:8px;"),

    tags$div(class = "sidebar-section-header", "Motif Discovery"),

    sliderInput(
      inputId = "motif_length",
      label   = "Motif Length (k)",
      min = 3, max = 15, value = 6, step = 1
    ),

    sliderInput(
      inputId = "min_freq",
      label   = "Min Sequence Frequency",
      min = 0, max = 1, value = 0.05, step = 0.01
    ),

    sliderInput(
      inputId = "min_count",
      label   = "Min Occurrence Count",
      min = 1, max = 50, value = 2, step = 1
    ),

    sliderInput(
      inputId = "top_n_motifs",
      label   = "Top N Motifs",
      min = 5, max = 100, value = 20, step = 5
    ),

    actionButton("run_discovery", "Discover Motifs",
                 icon  = icon("search"),
                 class = "btn-success btn-block",
                 style = "margin-bottom:8px;"),

    fileInput(
      inputId   = "control_fasta_file",
      label     = "Optional Control FASTA",
      accept    = c(".fasta", ".fa", ".fna", ".txt"),
      placeholder = "Upload a control/background FASTA"
    ),

    helpText(
      "Upload a second FASTA to compare foreground and background motif counts.",
      "Select 'Control dataset' as the background model to enable group comparisons."
    ),

    selectInput(
      inputId = "background_model",
      label   = "Background Model",
      choices = c(
        "Uniform" = "uniform",
        "Mononucleotide composition" = "mononucleotide",
        "Control dataset" = "control"
      ),
      selected = "uniform"
    ),

    selectInput(
      inputId = "p_adjust_method",
      label   = "P-value correction",
      choices = c("Bonferroni" = "bonferroni", "FDR" = "fdr"),
      selected = "bonferroni"
    ),

    tags$div(class = "sidebar-section-header", "Visualization"),

    sliderInput(
      inputId = "n_top_display",
      label   = "Motifs in Plots",
      min = 5, max = 50, value = 15, step = 5
    ),

    sliderInput(
      inputId = "pos_bins",
      label   = "Position Bins (heatmap)",
      min = 5, max = 40, value = 20, step = 5
    ),

    tags$div(class = "sidebar-section-header", "Status"),

    verbatimTextOutput("status_message"),
    br()
  ),

  # ---- Body ----
  dashboardBody(

    # Value boxes row
    fluidRow(
      valueBoxOutput("vbox_sequences", width = 3),
      valueBoxOutput("vbox_mean_length", width = 3),
      valueBoxOutput("vbox_motifs_found", width = 3),
      valueBoxOutput("vbox_significant", width = 3)
    ),

    # Main tab panel
    tabBox(
      width = 12,
      id    = "main_tabs",

      # ----- TAB 1: Overview -----
      tabPanel(
        title = tagList(icon("chart-bar"), " Overview"),
        value = "tab_overview",

        fluidRow(
          box(
            width = 6, title = "Sequence Length Distribution",
            solidHeader = TRUE, collapsible = TRUE,
            plotlyOutput("plot_length_dist", height = "320px")
          ),
          box(
            width = 6, title = "Nucleotide Frequencies",
            solidHeader = TRUE, collapsible = TRUE,
            plotlyOutput("plot_nt_freq", height = "320px")
          )
        ),
        fluidRow(
          box(
            width = 4, title = "GC Content",
            solidHeader = TRUE, collapsible = TRUE,
            plotlyOutput("plot_gc", height = "280px")
          ),
          box(
            width = 8, title = "Dataset Summary",
            solidHeader = TRUE, collapsible = TRUE,
            plotlyOutput("plot_summary_bar", height = "280px")
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Preprocessing Log",
            solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
            verbatimTextOutput("preprocess_log")
          )
        )
      ),

      # ----- TAB 2: Motif Logos -----
      tabPanel(
        title = tagList(icon("dna"), " Motif Logos"),
        value = "tab_logos",

        fluidRow(
          box(
            width = 12,
            title = "Select Motif to Display Logo",
            solidHeader = TRUE,
            selectInput("logo_motif_select",
                        label = NULL,
                        choices = NULL,
                        width   = "400px"),
            plotOutput("plot_logo_selected", height = "250px")
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Top 6 Motif Logos",
            solidHeader = TRUE, collapsible = TRUE,
            plotOutput("plot_logos_grid", height = "500px")
          )
        )
      ),

      # ----- TAB 3: Frequency Analysis -----
      tabPanel(
        title = tagList(icon("signal"), " Frequency Analysis"),
        value = "tab_freq",

        fluidRow(
          box(
            width = 8, title = "Top Motifs by Frequency",
            solidHeader = TRUE, collapsible = TRUE,
            plotlyOutput("plot_freq_bar", height = "500px")
          ),
          box(
            width = 4, title = "Enrichment Volcano",
            solidHeader = TRUE, collapsible = TRUE,
            plotlyOutput("plot_volcano", height = "500px")
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Full Motif Table",
            solidHeader = TRUE, collapsible = TRUE,
            DTOutput("table_motifs")
          )
        )
      ),

      # ----- TAB 4: Enrichment Heatmaps -----
      tabPanel(
        title = tagList(icon("th"), " Enrichment Heatmaps"),
        value = "tab_heat",

        fluidRow(
          box(
            width = 12, title = "Positional Enrichment Heatmap",
            solidHeader = TRUE,
            plotlyOutput("plot_pos_heatmap", width = "100%", height = "500px")
          )
        ),
        fluidRow(
          box(
            width = 12, title = "Motif Positional Density",
            solidHeader = TRUE,
            plotlyOutput("plot_density", width = "100%", height = "380px")
          )
        )
      ),

      # ----- TAB 5: Downloads -----
      tabPanel(
        title = tagList(icon("download"), " Downloads"),
        value = "tab_downloads",

        fluidRow(
          box(
            width = 12,
            title = "Export Results",
            solidHeader = TRUE,

            fluidRow(
              column(4,
                tags$h4("CSV Reports", style = "color:#58A6FF;"),
                downloadButton("dl_motif_csv",  "Motif Table (.csv)",
                               class = "btn-primary btn-block",
                               style = "margin-bottom:8px;"),
                downloadButton("dl_stats_csv",  "Sequence Stats (.csv)",
                               class = "btn-primary btn-block",
                               style = "margin-bottom:8px;"),
                downloadButton("dl_pos_csv",    "Positional Data (.csv)",
                               class = "btn-primary btn-block")
              ),
              column(4,
                tags$h4("Plot Downloads", style = "color:#58A6FF;"),
                downloadButton("dl_plot_freq",   "Frequency Bar Plot (.png)",
                               class = "btn-warning btn-block",
                               style = "margin-bottom:8px;"),
                downloadButton("dl_plot_volcano","Volcano Plot (.png)",
                               class = "btn-warning btn-block",
                               style = "margin-bottom:8px;"),
                downloadButton("dl_plot_logos",  "Logos Grid (.png)",
                               class = "btn-warning btn-block")
              ),
              column(4,
                tags$h4("Full Report", style = "color:#58A6FF;"),
                tags$p("Generate a comprehensive PDF report including all results,
                         tables, and figures.",
                        style = "color:#8B949E; font-size:12px;"),
                actionButton("generate_pdf", "Generate PDF Report",
                             icon  = icon("file-pdf"),
                             class = "btn-danger btn-block",
                             style = "margin-bottom:8px;"),
                uiOutput("pdf_download_ui")
              )
            )
          )
        ),

        fluidRow(
          box(
            width = 12,
            title = "Analysis Parameters Used",
            solidHeader = TRUE, collapsible = TRUE,
            verbatimTextOutput("params_summary")
          )
        )
      )
    )
  )
)
