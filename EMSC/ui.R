library(shiny)
library(tidyverse)
library(viridis)
library(shinydashboard)
library(dashboardthemes)



ui <- fluidPage(
  h3("Environmental Mutation Calculator"),
  "Calculate mutation enrichment with user-defined gene sets versus randomly sampled genes using exact bionomial tests. Substitution mutations called from WGS of iPSC clones dosed with 12 classes of mutagenic carcinogens. Data from Kucab et al., 2019. PMID: 30982602. See Baker et al., 2021 (https://doi.org/10.1101/2021.12.17.473207) for more details.",
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
      selectInput("select", h3("Sequence type"), choices = list("Entire gene" = 1,
                                                                "Coding sequence only" = 2), selected = 1),
      "Prepared and maintained by Brennan Baker and Brandon Pearson, Columbia University. Contact: blp2125@cumc.columbia.edu",
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", plotOutput("heatmap_gene", width = "1500px", height = "1500px"),textOutput("text1")),
        tabPanel("Valid genes",textOutput("text3"), tableOutput("contents3")),
        tabPanel("Results table", tableOutput("contents2")),
        tabPanel("Download results", downloadButton("download1"))
      )
    )))

