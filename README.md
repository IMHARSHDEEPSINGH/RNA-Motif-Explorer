# RNA Motif Discovery Explorer

> An interactive R/Shiny bioinformatics platform for discovering, analyzing, and
> visualizing recurring RNA motifs in large FASTA datasets.

---

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Usage Guide](#usage-guide)
6. [Pipeline Details](#pipeline-details)
7. [Output Files](#output-files)
8. [Deployment](#deployment)
9. [API Integration](#api-integration)
10. [Future Extensibility](#future-extensibility)
11. [Troubleshooting](#troubleshooting)

---

## Overview

The **RNA Motif Discovery Explorer** accepts large RNA sequence datasets in FASTA
format and runs a fully automated analysis pipeline:

- **Preprocessing** — validation, cleaning, deduplication, quality metrics
- **Motif Discovery** — k-mer extraction, frequency counting, enrichment scoring
- **Statistical Testing** — binomial tests, Bonferroni correction, combined ranking
- **Secondary Structure Prediction** — RNA folding, structure element mapping, motif enrichment by structure
- **Visualization** — sequence logos, frequency bars, volcano plots, positional
  heatmaps, structure plots, enrichment heatmaps
- **Export** — CSV tables, PNG plots, PDF reports

All of this is wrapped in a dark-themed **Shiny Dashboard** with interactive
Plotly plots, searchable DataTables, and one-click downloads.

---

## Architecture

```
RNA-Motif-Explorer/
├── data/
│   ├── input_FASTA/         # Raw FASTA uploads (sample included)
│   └── processed/           # Intermediate processed data
├── scripts/
│   ├── preprocessing.R      # FASTA reading, validation, cleaning, stats
│   ├── motif_detection.R    # k-mer counts, enrichment, p-values, PWMs
│   ├── structure_prediction.R    # Secondary structure prediction (RNAstructure)
│   ├── motif_structure_analysis.R # Motif-structure mapping & enrichment
│   ├── visualization.R      # All ggplot2 / plotly plot functions
│   └── report_generation.R  # CSV/PNG/PDF export utilities
├── shiny_app/
│   ├── ui.R                 # Shiny Dashboard UI definition
│   ├── server.R             # Reactive server logic
│   └── app.R                # Entry point
├── api/
│   └── plumber_api.R        # REST API (Plumber) for web integration
├── outputs/
│   ├── logos/               # Exported sequence logos
│   ├── plots/               # Exported plots
│   └── reports/             # CSV reports and PDF output
├── install_dependencies.R   # Package installer
├── Dockerfile               # Docker container definition
└── README.md
```

### Design Principles

- **Modularity** — Each script has a single responsibility. Backend functions
  can be tested independently.
- **Reactivity** — Shiny reactivity is isolated to `server.R`; all scientific
  logic lives in the `scripts/` modules.
- **Separation of concerns** — UI, server, data, and API layers are cleanly
  separated so any layer can be swapped (e.g., replace Shiny with React).
- **No hardcoded values** — All parameters flow from user inputs or function
  arguments.

---

## Installation

### Requirements

- R ≥ 4.2
- For PDF reports: a LaTeX distribution (TinyTeX recommended)

### Install R Packages

```r
# From the project root directory:
Rscript install_dependencies.R
```

This installs:

**CRAN:** shiny, shinydashboard, shinyjs, DT, plotly, ggplot2, ggseqlogo,
dplyr, tidyr, tibble, stringr, writexl, rmarkdown, knitr, gridExtra, scales,
RColorBrewer, plumber

**Bioconductor:** Biostrings, universalmotif

### TinyTeX (for PDF reports — optional)

```r
install.packages("tinytex")
tinytex::install_tinytex()
```

#### Pandoc requirement

PDF report generation also requires `pandoc` 1.12.3 or newer. Most RStudio installations include pandoc, but if you use plain R you may need to install it separately from https://pandoc.org.

## Testing

```r
Rscript install_dependencies.R
Rscript run_tests.R
```

---

## Quick Start

```r
# 1. Install dependencies
Rscript install_dependencies.R

# 2. Launch the Shiny app from project root
shiny::runApp("shiny_app/")

# 3. In the browser:
#    - Upload data/input_FASTA/sample_rna.fasta
#    - Click "Preprocess Sequences"
#    - Adjust motif length (default: 6)
#    - Click "Discover Motifs"
#    - Explore the tabs
```

---

## Usage Guide

### Sidebar Controls

| Control | Description | Default |
|---|---|---|
| Upload FASTA | `.fasta`, `.fa`, `.fna`, `.txt` files | — |
| Min Sequence Length | Filter sequences shorter than N nt | 10 |
| Motif Length (k) | k-mer length to extract (3–15) | 6 |
| Min Sequence Frequency | Motif must appear in at least this fraction of sequences | 0.05 |
| Min Occurrence Count | Minimum total occurrence count | 2 |
| Top N Motifs | Number of top-ranked motifs to analyse | 20 |
| Motifs in Plots | Number to show in frequency/heatmap plots | 15 |
| Position Bins | Resolution of positional heatmap | 20 |

### Tabs

| Tab | Content |
|---|---|
| Overview | Length distribution, nucleotide frequencies, GC content, dataset summary |
| Motif Logos | Interactive sequence logo per motif; 6-logo grid |
| Frequency Analysis | Ranked bar chart, volcano plot, full searchable table |
| Enrichment Heatmaps | Positional heatmap, density curves |
| Secondary Structure | Structure statistics, motif enrichment by structure type, structure-motif correlation |
| Downloads | CSV/PNG/PDF downloads |

---

## Pipeline Details

### Step 1 — Preprocessing (`preprocessing.R`)

```
read_fasta_file()           → RNAStringSet
validate_rna_sequences()    → remove invalid characters
clean_rna_sequences()       → uppercase, T→U, filter by length
remove_duplicate_sequences()→ deduplicate by sequence content
calculate_sequence_stats()  → per-sequence and global metrics
```

**Outputs:** sequence length, A/U/G/C frequencies, GC content, dataset summary.

### Step 2 — Motif Discovery (`motif_detection.R`)

```
count_kmers_dataset()          → total & sequence counts per k-mer
calculate_enrichment_scores()  → log2(observed/expected) vs uniform model
calculate_pvalues()            → binomial test + Bonferroni correction
compute_positional_distribution() → relative positions across sequences
build_pwm_from_kmers()         → Position Weight Matrix (universalmotif)
rank_motifs()                  → combined frequency + enrichment + significance rank
```

**Background model:** Uniform (p = 0.25 per nucleotide).

### Step 3 — Secondary Structure Prediction (`structure_prediction.R`)

```
predict_secondary_structure()    → RNAstructure MFE folding
predict_structures_dataset()     → bulk structure prediction for all sequences
parse_structure_elements()       → count stems, loops, bulges, hairpins
map_position_to_structure()      → assign structure type to each position
structure_statistics()           → aggregate structure metrics
```

**Details:** Uses the RNAstructure R package to compute minimum free energy (MFE) secondary structures
in dot-bracket notation. Each position is classified as stem (paired), loop, bulge, or hairpin.

### Step 4 — Motif-Structure Analysis (`motif_structure_analysis.R`)

```
find_motif_positions_in_sequence()      → locate motif occurrences
classify_motif_positions()              → assign structure context to each occurrence
map_motifs_to_structures_dataset()      → full dataset mapping
calculate_motif_structure_enrichment()  → Fisher's exact test, log-odds ratios
get_motif_structure_distribution()      → count motifs per structure type
compute_structure_correlation_metrics() → correlation between motif presence & structure
```

**Outputs:** Enrichment p-values, distributions, correlations showing which structural elements
preferentially contain specific motifs.

### Step 5 — Visualization (`visualization.R`)

All plots use a consistent dark theme (`THEME_RNA()`) and the RNA color scheme:
- A = `#E74C3C` (red) · U = `#3498DB` (blue) · G = `#2ECC71` (green) · C = `#F39C12` (orange)

---

## Output Files

| File | Description |
|---|---|
| `motif_table.csv` | Full ranked motif table with all statistics |
| `sequence_stats.csv` | Per-sequence statistics |
| `positional_distribution.csv` | Motif positions across sequences |
| `structure_predictions.csv` | Predicted secondary structures (dot-bracket notation) |
| `motif_structure_mapping.csv` | Motif occurrences mapped to structure elements |
| `structure_enrichment.csv` | Enrichment scores and p-values for motifs in each structure type |
| `motif_structure_distribution.csv` | Motif frequency by structure type |
| `structure_correlation.csv` | Correlation between motifs and structural elements |
| `freq_bar_k{k}.png` | Frequency bar chart |
| `volcano_k{k}.png` | Enrichment volcano plot |
| `logos_grid_k{k}.png` | Sequence logos grid |
| `rna_motif_report.pdf` | Full PDF report |

---

## Deployment

### Local (development)

```r
shiny::runApp("shiny_app/", port = 3838)
```

### shinyapps.io

```r
library(rsconnect)
rsconnect::deployApp("shiny_app/",
                     appName = "rna-motif-explorer")
```

### Docker

```bash
# Build
docker build -t rna-motif-explorer .

# Run
docker run -p 3838:3838 rna-motif-explorer

# Access at http://localhost:3838
```

### Shiny Server (Linux)

```bash
# Copy shiny_app/ to Shiny Server directory
sudo cp -r shiny_app/ /srv/shiny-server/rna-motif-explorer/
sudo cp -r scripts/   /srv/shiny-server/rna-motif-explorer/

# Restart Shiny Server
sudo systemctl restart shiny-server
```

---

## API Integration

The `api/plumber_api.R` exposes the full pipeline as a REST API:

```r
# Launch API server
library(plumber)
plumber::plumb("api/plumber_api.R")$run(port = 8080)
```

### Endpoints

| Method | Endpoint | Description |
|---|---|---|
| GET | `/health` | Health check |
| GET | `/config` | Default parameters |
| POST | `/preprocess` | Upload + preprocess FASTA |
| POST | `/discover` | Run motif discovery on sequences |
| POST | `/top_kmers` | Lightweight top k-mer endpoint |
| POST | `/structure/predict` | Predict secondary structure for sequences |
| POST | `/structure/motif_enrichment` | Get motif enrichment across structure types |
| POST | `/structure/correlation` | Compute motif-structure correlation metrics |

### Example (curl)

```bash
curl -X GET http://localhost:8080/health

curl -X POST http://localhost:8080/top_kmers \
  -H "Content-Type: application/json" \
  -d '{"sequences": ["AUGCGUACGAUGCUA", "GCUAGCUAGCUAGCU"], "k": 5, "top_n": 5}'
```

---

## Future Extensibility

### Planned Enhancements

| Feature | Implementation Path |
|---|---|
| React frontend | Connect to Plumber API; replace Shiny UI |
| Cloud deployment | AWS/GCP with Docker + load balancer |
| Multi-sample comparison | Add sample grouping + differential enrichment |
| MEME/JASPAR integration | Import known RNA motif databases |
| ~~Structure prediction~~ | ✅ **DONE** — RNAstructure integration with enrichment analysis |
| Alignment-based logos | Smith-Waterman local alignment before PWM |
| GPU k-mer counting | Rewrite `count_kmers_dataset` with data.table or C++ via Rcpp |
| Authentication | Add shinymanager for user login |
| Job queue | Add future/promises for async large jobs |

### Adding a New Analysis Module

1. Create `scripts/my_analysis.R` with documented functions.
2. Add `source("../scripts/my_analysis.R")` in `server.R`.
3. Add UI elements in `ui.R`.
4. Expose via new Plumber endpoint in `api/plumber_api.R`.

---

## Troubleshooting

**"No k-mers pass the current filters"**
→ Lower `Min Sequence Frequency` or `Min Occurrence Count` in the sidebar.

**Logo not rendering**
→ Ensure `ggseqlogo` is installed. The fallback shows a plain ggplot.

**PDF generation fails**
→ Install TinyTeX: `tinytex::install_tinytex()`. The CSV reports always work.

**Very slow on large datasets (>10,000 sequences)**
→ Increase `Min Sequence Frequency` to filter rare k-mers faster. Consider
using `data.table` for k-mer counting (planned upgrade).

**DNA FASTA input needs RNA conversion**
→ The cleaner in `preprocessing.R` preserves uploaded bases, then auto-converts T→U.

---

## License

MIT License — see LICENSE file.
