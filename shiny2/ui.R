library(shiny)
library(gapminder)
library(plotly)
library(ggplot2)
library(readr)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(randomcoloR)
library(tsne)

# Define UI for application that draws a scatter plot
shinyUI(fluidPage(
          sidebarPanel(
              useShinyjs(),
              width = 2,
              style = "height: 750px;",
              
              div(style="display: inline-block;vertical-align:top; width: 100px;",
                selectInput("gprotein", "G-prots/B-arrs",
                            choices = c('GoB', 'GNAI2', 'GNA14', 'GNA12', 'GoA', 'Barr2-GRK2', 'GNAI3', 'GNA15', 'Barr2', 'Barr1-GRK2', 'GNAQ', 'GNAO1', 'GNAI1', 'GNAS', 'GNAZ', 'GNA11', 'GNA13', 'GNAL'),
                            selected = "GNAI3")
                ),
              
              div(style="display: inline-block;vertical-align:top; width: 100px;",
                  selectInput("pca_type", "PCA type", choices = c('All'), selected = "All")
              ),
              
              div(style="vertical-align:top; width: 200px;",
                radioButtons("option", "Default layer to display:", c("GPCRome" = "gpcrome", "Best KM" = "bestKM", "Best ML" = "bestML", "None" = "none")),
              ),
              
              div(style="display: inline-block; width: 100px;",
                  uiOutput("layer")
              ),

              div(style="display: inline-block; width: 100px;",
                  selectInput("assay", "Color by",
                              choices = c('Shedding', 'IUPHAR', 'ebBRET', 'STRING', 'Family', 'Class'),
                              selected = "Shedding")
              ),
              
              div(style="display: inline-block; width: 100px;",
                div(style="display: inline-block; padding-left: 3%;", strong('Taste recep?', style="display: inline-block; padding-bottom: 5px; vertical-align:top; width: 100px;")),
                switchInput(inputId = "taste", label="Show", value = TRUE, onLabel = "Yes", offLabel = "No")
              ),
              
              div(style="display: inline-block; width: 10px;"),
              div(style="display: inline-block; width: 100px;",
                div(style="display: inline-block; padding-left: 3%;", strong('Scaled?', style="display: inline-block; padding-bottom: 5px; vertical-align:top; width: 75px;")),
                switchInput(inputId = "scaled", label="Scaled", value = FALSE, onLabel = "Yes", offLabel = "No")
              ),
              
              div(style="vertical-align:top; width: 150px;",
                uiOutput("numClusters")
                #numericInput('numClusters', 'Cluster count', 3, min = 1, max = 9)
                #selectInput("numClusters", "Cluster count",
                 #           choices = c(1:10),
                  #          selected = 3)
              ),
              
              hr(style = "border: 1px solid;"),
              
              div(style="display: inline-block; padding-left: 3%;", strong('Color all GPCR families?', style="display: inline-block; padding-bottom: 5px; vertical-align:top;")),
              div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
              switchInput(inputId = "all", label="All", value = TRUE, onLabel = "Yes", offLabel = "No"),
              
              div(style="display: inline-block; padding-left: 3%;", strong('OR', style="display: inline-block; padding-bottom: 10px; vertical-align:top;")),
              selectizeInput("family", "Select one or more families:",
                             choices = NULL, multiple = TRUE),
          ),
        mainPanel(
          tabsetPanel(
            id = "tabs",
            tabPanel(title = "PCA",
                      column(10,
                         align = "center",
                         plotlyOutput(outputId = "scatter", width = '100%', height = '100%'),
                         DT::dataTableOutput("table", width = '75%')
                      )
                    ),
            tabPanel(title = "t-SNE",
                     column(10,
                            align = "center",
                            plotlyOutput(outputId = "scatter2", width = '100%', height = '100%'),
                            DT::dataTableOutput("table2", width = '75%')
                     )
            ),
            tabPanel(title = "UMAP",
                     column(10,
                            align = "center",
                            plotlyOutput(outputId = "scatter3", width = '100%', height = '100%'),
                            DT::dataTableOutput("table3", width = '75%')
                     )
            )
  
          )
        )
)
)
