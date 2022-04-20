library(shiny)
library(tidyverse)
library(viridis)
library(shinydashboard)
library(dashboardthemes)



ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", h3("Input a csv file"),
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      "File may have one or more gene sets. Each column should have the gene set name in the first row, with one gene per row thereafter. Example file input below:",
      tableOutput("example"),
      selectInput("select", h3("Sequence type"), choices = list("Full gene body" = 1,
                                                                "Coding sequence only" = 2), selected = 1),
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", plotOutput("heatmap_gene", width = "1500px", height = "1500px"),textOutput("text1")),
        tabPanel("Plot", plotOutput("contents", width = "1500px", height = "1500px"),textOutput("text2")), 
        tabPanel("Results table", tableOutput("contents2")),
        tabPanel("Download results", downloadButton("download1"))
      )
    )))

