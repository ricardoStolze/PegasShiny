
cat("\014")

library(shiny)
library(vcfR)
library(pegas)
library(pheatmap)
library(dplyr)
library(igraph)
library(shinycssloaders)
library(shinyBS)
library(graphics)

source('appUI.R')
source('appServer.R')

# Run the application
shinyApp(ui = ui, server = server)
