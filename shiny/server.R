library(shiny)
library(plotly)
library(ggplot2)
library(readr)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)
library(randomcoloR)

plot <- function(gprotein, assay, family, all, pca_type, taste) {
  if (pca_type == "GPCRome" && taste == FALSE){
    gproteinPath <- paste("33layer_PCA_without_taste/", gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "GPCRome" && taste == TRUE){
    gproteinPath <- paste("33layer_PCA/", gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "Best" && taste == TRUE){
    gproteinPath <- paste("best_pca/", gprotein, "_", pca_type, ".tsv", sep="")
  }
  else {
    gproteinPath <- paste("best_pca_without_taste/", gprotein, "_", pca_type, ".tsv", sep="")
  }
  
  data <- read_tsv(gproteinPath)
  print (family)
  if (all == TRUE && assay %in% c('Shedding', 'ebBRET', 'IUPHAR', 'STRING')) {
    p <- data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color=assay)) +
      geom_point()+
      scale_color_manual(values=c('Coupled' = "forestgreen", 'Not-coupled' = "darkred", "-" = "grey")) 
  }
  else if (is.null(family) == FALSE && assay %in% c('Family')) {
    new_data <- mutate(data, Selected = ifelse(`Family` %in% family, `Family`, "Others"))
    print (levels(new_data$Selected))
    
    families <- unique(new_data$Selected)
    families <- families[families != "Others"]
    n <- length(families)
    #palette <- distinctColorPalette(n)
    palette <- randomColor(n, luminosity="bright")
    
    families <- c(families, 'Others')
    palette <- c(palette, 'grey')
    names(palette) <- families
    
    p <- new_data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color="Selected")) +
      geom_point()+
      scale_color_manual(name=new_data$Selected, values=palette) 
  }
  else if (all == TRUE && assay %in% c('Class')) {
    print ('class')
    p <- data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color=assay)) +
      geom_point() +
      scale_color_manual(values=c('classA' = "forestgreen", 'classB' = "darkred", "classC" = "chocolate", "Frizzeled" = "gold", "Taste" = "pink","Other" = "grey")) 
  }
  else {
    p <- data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color=assay)) +
      geom_point()
  }
  return(ggplotly(p, height = 900))
}

# Define server logic required to draw a scatter plot
shinyServer(
  function(input, output, session) {
    
    gproteinPath <- paste("33layer_PCA/", "GNAS", "_", "GPCRome", ".tsv", sep="")
    data <- read_tsv(gproteinPath)
    updateSelectizeInput(session, 'family', choices = levels(factor(data$`Family`)), server = TRUE)
    
    observeEvent(input$assay, {
      print (input$assay)
      if (input$assay == "Family"){
        shinyjs::enable("all")
      }
      else {
        shinyjs::disable("all")
      }
    })
    
    observeEvent(input$all, {
      if (input$all == TRUE){
        shinyjs::disable("family")
      }
      else {
        shinyjs::enable("family")
      }
    })
    observeEvent(input$assay, {
    output$scatter <- renderPlotly({plot(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$taste)})
    })
    
    }
)
