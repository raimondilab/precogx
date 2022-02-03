library(shiny)
library(gapminder)
library(plotly)
library(ggplot2)
library(readr)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(randomcoloR)

# Define UI for application that draws a scatter plot
shinyUI(navbarPage("PCA",
          tabPanel("Scatter plot",
            fluidPage(
              tags$style(
                HTML('
                     #structure_panel {
                     display: flex;
                     align-items: center;
                     justify-content: center;
                     top: 50%;
                     }
                     ')
                  ),
              fluidRow(
                column(2,
                      align="center",
                      wellPanel(
                                useShinyjs(),
                                style = "height: 750px;",
                                
                                div(style="display: inline-block; vertical-align:top;",
                                  selectInput("gprotein", "G-proteins/B-arrs",
                                              choices = c('GoB', 'GNAI2', 'GNA14', 'GNA12', 'GoA', 'Barr2-GRK2', 'GNAI3', 'GNA15', 'Barr2', 'Barr1-GRK2', 'GNAQ', 'GNAO1', 'GNAI1', 'GNAS', 'GNAZ', 'GNA11', 'GNA13', 'GNAL'),
                                              selected = "GNAI3")
                                  ),
                                
                                div(style="vertical-align:top; width: 150px;",
                                    selectInput("pca_type", "PCA type",
                                                choices = c('GPCRome', 'Best'),
                                                selected = "GPCRome")
                                ),
                                
                                div(style="display: inline-block; padding-left: 3%;", strong('Taste receptors?', style="display: inline-block; padding-bottom: 5px; vertical-align:top;")),
                                switchInput(inputId = "taste", label="Show", value = TRUE, onLabel = "Yes", offLabel = "No"),
   
                                div(style="vertical-align:top; width: 150px;",
                                    selectInput("assay", "Color by",
                                                choices = c('Shedding', 'IUPHAR', 'ebBRET', 'STRING', 'Family', 'Class'),
                                                selected = "Shedding")
                                ),
                                
                                div(style="display: inline-block; padding-left: 3%;", strong('Scaled data?', style="display: inline-block; padding-bottom: 5px; vertical-align:top;")),
                                switchInput(inputId = "scaled", label="Scaled", value = TRUE, onLabel = "Yes", offLabel = "No"),
                                
                                div(style="vertical-align:top; width: 150px;",
                                  numericInput('numClusters', 'Cluster count', 3, min = 1, max = 9)
                                ),
                                
                                hr(style = "border: 1px solid;"),
                                
                                div(style="display: inline-block; padding-left: 3%;", strong('Color all GPCR families?', style="display: inline-block; padding-bottom: 5px; vertical-align:top;")),
                                div(style="display: inline-block; padding-left: 3%;", actionButton("inputhelp", "", icon = icon("question-circle"))),
                                switchInput(inputId = "all", label="All", value = TRUE, onLabel = "Yes", offLabel = "No"),
                                
                                div(style="display: inline-block; padding-left: 3%;", strong('OR', style="display: inline-block; padding-bottom: 10px; vertical-align:top;")),
                                selectizeInput("family", "Select one or more families:",
                                               choices = NULL, multiple = TRUE),
                            )
                      ),
                column(10,
                       align = "center",
                       plotlyOutput(outputId = "scatter", width = '100%', height = '100%'),
                       DT::dataTableOutput("table", width = '75%')
                )
              )
              
            )
            )
  )
)
