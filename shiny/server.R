library(shiny)
library(gapminder)
library(plotly)
library(ggplot2)
library(readr)
library(shinyjs)
library(dplyr)
library(shinythemes)
library(shinyWidgets)
library(randomcoloR)
library(tidyverse)

RV <- reactiveValues(data = data.frame())

table2 <- function(gprotein, assay, family, all, pca_type, taste, numClusters, scaled) {
  if (scaled == TRUE) {
    scaled = "_scaled"
  }
  else {
    scaled = ""
  }
  if (pca_type == "GPCRome" && taste == FALSE){
    gproteinPath <- paste("33layer_PCA_without_taste", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "GPCRome" && taste == TRUE){
    gproteinPath <- paste("33layer_PCA", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "Best" && taste == TRUE){
    gproteinPath <- paste("best_PCA", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else {
    gproteinPath <- paste("best_PCA_without_taste", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  
  data <- read_tsv(gproteinPath)
  selectedData <- reactive({
    #iris[, c(input$xcol, input$ycol)]
    data[, c('PC1', 'PC2')]
  })
  
  clusters <- reactive({
    kmeans(selectedData(), numClusters)
  })
  
  clust <- clusters()$cluster
  data <- mutate(data, clust=clust)
  var <- unique(data[[assay]])
  print (var)
  df <- data_frame()
  for(i in 1:length(var)) {       # for-loop over columns
    clas <- var[i]
    newData <- data %>%
      filter((!!sym(assay)) %in% var[i]) %>%
      group_by(clust) %>%
      summarise(n = n()) %>%
      mutate(Freq = n/sum(n)) %>%
      mutate(clas = var[i])
    names(newData)[names(newData) == "Freq"] <- var[i]
    if (dim(df) == 0){
      df <- newData[c('clust', var[i])]
    }
    else {
      df <- full_join(df, newData[c('clust', var[i])], by="clust", quiet=TRUE)
    }
    
  }
  df <- df[order(df$clust),]
  print (df)
  return(df)
}

plot <- function(gprotein, assay, family, all, pca_type, taste, numClusters, scaled) {
  print (taste)
  if (scaled == TRUE) {
    scaled = "_scaled"
  }
  else {
    scaled = ""
  }
  if (pca_type == "GPCRome" && taste == FALSE){
    gproteinPath <- paste("33layer_PCA_without_taste", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "GPCRome" && taste == TRUE){
    gproteinPath <- paste("33layer_PCA", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "Best" && taste == TRUE){
    gproteinPath <- paste("best_PCA", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else {
    gproteinPath <- paste("best_PCA_without_taste", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  
  data <- read_tsv(gproteinPath)
  selectedData <- reactive({
    #iris[, c(input$xcol, input$ycol)]
    data[, c('PC1', 'PC2')]
  })
  
  clusters <- reactive({
    kmeans(selectedData(), numClusters)
  })
  #print (clusters()$cluster)
  clust <- clusters()$cluster
  data <- mutate(data, clust=clust)
  ###############################################################################
  var <- unique(data[[assay]])
  print (var)
  df <- data_frame()
  for(i in 1:length(var)) {       # for-loop over columns
    clas <- var[i]
    newData <- data %>%
      filter((!!sym(assay)) %in% var[i]) %>%
      group_by(clust) %>%
      summarise(n = n()) %>%
      mutate(Freq = round(n/sum(n), digits=2)) %>%
      mutate(clas = var[i])
    names(newData)[names(newData) == "Freq"] <- var[i]
    if (dim(df) == 0){
      df <- newData[c('clust', var[i])]
    }
    else {
      df <- full_join(df, newData[c('clust', var[i])], by="clust", quiet=TRUE)
    }
    
  }
  df <- df[order(df$clust),]
  ###############################################################################
  #print (data)
  print (family)
  if (all == TRUE && assay %in% c('Shedding', 'ebBRET', 'IUPHAR', 'STRING')) {
    p <- data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color=assay, size=1)) +
      geom_point(aes(shape = factor(clust), size=1))+
      scale_color_manual(values=c('Coupled' = "forestgreen", 'Not-coupled' = "darkred", "-" = "grey"))+
      scale_shape_manual(values=c("1" = 15, "2" = 16, "3" = 17, "4" = 18, "5" = 19, "6" = 20, "7" = 21, "8" = 22, "9" = 23)) 
  }
  else if (is.null(family) == FALSE && assay %in% c('Family')) {
    new_data <- mutate(data, Selected = ifelse(`Family` %in% family, `Family`, "Others"))
    #print (levels(new_data$Selected))
    
    families <- unique(new_data$Selected)
    families <- families[families != "Others"]
    n <- length(families)
    #palette <- distinctColorPalette(n)
    palette <- randomColor(n, luminosity="bright")
    
    families <- c(families, 'Others')
    palette <- c(palette, 'grey')
    names(palette) <- families
    #print (palette)
    
    p <- new_data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color="Selected")) +
      #geom_point()+
      geom_point(aes(shape = factor(clust), size=1))+
      scale_color_manual(values=palette)+
      scale_shape_manual(values=c("1" = 15, "2" = 16, "3" = 17, "4" = 18, "5" = 19, "6" = 20, "7" = 21, "8" = 22, "9" = 23)) 
      #scale_color_manual(name=new_data$Selected, values=palette) 
  }
  else if (all == TRUE && assay %in% c('Class')) {
    #print ('class')
    p <- data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color=assay)) +
      #geom_point() +
      geom_point(aes(shape = factor(clust), size=1))+
      scale_color_manual(values=c('classA' = "forestgreen", 'classB' = "darkred", "classC" = "chocolate", "Frizzeled" = "gold", "Taste" = "pink","Other" = "grey"))+
      scale_shape_manual(values=c("1" = 15, "2" = 16, "3" = 17, "4" = 18, "5" = 19, "6" = 20, "7" = 21, "8" = 22, "9" = 23)) 
  }
  else {
    p <- data %>%
      ggplot(aes_string(text="Gene", "PC1", "PC2", color=assay)) +
      #geom_point()
      geom_point(aes(shape = factor(clust), size=1))+
      scale_shape_manual(values=c("1" = 15, "2" = 16, "3" = 17, "4" = 18, "5" = 19, "6" = 20, "7" = 21, "8" = 22, "9" = 23)) 
  }
  #print (taste)
  List <- list('df'= df, 'plot' = ggplotly(p, height = 600))
  return(List)
  #return(ggplotly(p, height = 500))
}

# Define server logic required to draw a scatter plot
shinyServer(
  function(input, output, session) {
    observeEvent(input$inputhelp, {
      showModal(modalDialog(
        title = "Help",
        HTML("This will work only when you select Color by as <kbd>Family</kbd>"
        )
      ))
    })
    
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
    observeEvent(c(input$assay, input$gprotein, input$pca_type, input$taste, input$scaled, input$scaled, input$numClusters), {
    #output$scatter <- renderPlotly({plot(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$taste, input$numClusters, input$scaled)})
    X <- plot(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$taste, input$numClusters, input$scaled)
    output$scatter <- renderPlotly({X$plot})
    RV$data <- X$df
    })
    
    #observeEvent(c(input$assay, input$gprotein, input$pca_type, input$taste, input$scaled, input$scaled, input$numClusters), {RV$data <- table(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$taste, input$numClusters, input$scaled)})
    output$table <- DT::renderDataTable({
      DT::datatable(RV$data,
                    filter = "top",
                    height = 1,
                    class = 'cell-border strip hover',
                    extensions = list("ColReorder" = NULL,
                                      "Buttons" = NULL,
                                      "FixedColumns" = list(leftColumns=1)),
                    options = list(
                      dom = 'lBRrftpi',
                      autoWidth=TRUE,
                      pageLength = -1,
                      lengthMenu = list(c(3, 5, 10, 15, 50, -1), c('3', '5', '10', '15','50', 'All')),
                      ColReorder = TRUE,
                      buttons =
                        list(
                          'copy',
                          'print',
                          list(
                            extend = 'collection',
                            buttons = c('csv', 'excel', 'pdf'),
                            text = 'Download'
                          )
                        )
                    )
      )
    }, server = TRUE)
    
    }
)
