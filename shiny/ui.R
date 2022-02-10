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
                                
                                div(style="display: inline-block;vertical-align:top; width: 100px;",
                                  selectInput("gprotein", "G-prots/B-arrs",
                                              choices = c('GoB', 'GNAI2', 'GNA14', 'GNA12', 'GoA', 'Barr2-GRK2', 'GNAI3', 'GNA15', 'Barr2', 'Barr1-GRK2', 'GNAQ', 'GNAO1', 'GNAI1', 'GNAS', 'GNAZ', 'GNA11', 'GNA13', 'GNAL'),
                                              selected = "GNAI3")
                                  ),
                                
                                div(style="display: inline-block;vertical-align:top; width: 100px;",
                                    selectInput("pca_type", "PCA type",
                                                choices = c('All'),
                                                selected = "All")
                                ),
                                
                                div(style="vertical-align:top; width: 200px;",
                                  radioButtons("option", "Default layer to display:",
                                               c("GPCRome" = "gpcrome",
                                                 "Best" = "best",
                                                 "None" = "none")),
                                ),
                                
                                div(style="display: inline-block; width: 100px;",
                                    uiOutput("layer")
                                    #selectInput("layer", "Layer",
                                    #            choices = c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33'),
                                    #            selected = "33")
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
                                  switchInput(inputId = "scaled", label="Scaled", value = TRUE, onLabel = "Yes", offLabel = "No")
                                ),
                                
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
