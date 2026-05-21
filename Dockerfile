# =============================================================================
# RNA Motif Discovery Explorer
# File: Dockerfile
# Purpose: Containerize the Shiny app for cloud deployment
# Build:  docker build -t rna-motif-explorer .
# Run:    docker run -p 3838:3838 rna-motif-explorer
# =============================================================================

FROM rocker/shiny:4.3.2

LABEL maintainer="RNA Motif Explorer"
LABEL description="RNA Motif Discovery Explorer - Shiny App"

# System dependencies for Bioconductor packages
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    pandoc \
    texlive-latex-base \
    texlive-fonts-recommended \
    texlive-latex-extra \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "\
  install.packages(c( \
    'BiocManager', 'shiny', 'shinydashboard', 'shinyjs', \
    'DT', 'plotly', 'ggplot2', 'ggseqlogo', \
    'dplyr', 'tidyr', 'tibble', 'stringr', \
    'writexl', 'rmarkdown', 'knitr', 'gridExtra', \
    'scales', 'RColorBrewer', 'plumber' \
  ), repos = 'https://cloud.r-project.org', quiet = TRUE); \
  BiocManager::install(c('Biostrings', 'universalmotif'), \
    ask = FALSE, quiet = TRUE)"

# Copy project files
COPY . /srv/shiny-server/rna-motif-explorer/

# Copy Shiny server config
RUN echo '\
server { \n\
  listen 3838; \n\
  location / { \n\
    site_dir /srv/shiny-server/rna-motif-explorer/shiny_app; \n\
    log_dir /var/log/shiny-server; \n\
    directory_index on; \n\
  } \n\
}' > /etc/shiny-server/shiny-server.conf

# Set working directory to app
WORKDIR /srv/shiny-server/rna-motif-explorer

EXPOSE 3838

CMD ["/usr/bin/shiny-server"]
