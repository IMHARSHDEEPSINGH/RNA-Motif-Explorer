# =============================================================================
# RNA Motif Discovery Explorer
# File: shiny_app/app.R
# Purpose: Application entry point — sources ui.R and server.R and launches.
# Run with: shiny::runApp("shiny_app/")  from the project root.
# =============================================================================

library(shiny)

# Source UI and server from the same directory
source("ui.R")
source("server.R")

shinyApp(ui = ui, server = server)
