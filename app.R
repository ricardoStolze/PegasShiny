

library(shiny)
library(vcfR)
library(pegas)
library(pheatmap)
library(dplyr)
library(igraph)
library(shinycssloaders)
library(shinyBS)
library(graphics)
library(shinyFeedback)
library(tools)

source('appUI.R')
source('appServer.R')

# Run the application
shinyApp(ui = ui, server = server)
