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

plot <- function(gprotein, assay, family, all, pca_type, layer, taste, numClusters, scaled, dimr) {
  if (is.null(numClusters) == TRUE){
    numClusters = 1
  }
  else if (grepl("Best", numClusters) == TRUE) {
    numClusters = strsplit(numClusters, split=" ")
    print ('Here')
    numClusters <- numClusters[[1]][1]
    print (numClusters)
  }
  print ('We are here')
  print (numClusters)
  #layer <- '33'
  if (is.null(layer) == TRUE){
    layer = "33"
  }
  else if (grepl("Best KM", layer) == TRUE) {
    layer = strsplit(layer, split=" ")
    print ('Here')
    layer <- layer[[1]][1]
    print (layer)
  }
  else if (grepl("Best ML", layer) == TRUE) {
    layer = strsplit(layer, split=" ")
    print ('Here')
    layer <- layer[[1]][1]
    print (layer)
  }
  print (taste)
  
  if (dimr == "PCA"){
    if (scaled == TRUE) {
      scaled = "_scaled"
    }
    else {
      scaled = ""
    }
  }
  else{
    scaled = ""
  }
  if (pca_type == "GPCRome" && taste == FALSE){
    gproteinPath <- paste("33layer_PCA_without_taste", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "GPCRome" && taste == TRUE){
    gproteinPath <- paste("33layer_PCA", scaled, '/', gprotein, "_", pca_type, ".tsv", sep="")
  }
  else if (pca_type == "All" && taste == FALSE && dimr == "PCA"){
    gproteinPath <- paste("all_layer_without_taste", scaled, '/', gprotein, "_", pca_type, "_", layer, ".tsv", sep="")
  }
  else if (pca_type == "All" && taste == TRUE  && dimr == "PCA"){
    gproteinPath <- paste("all_layers", scaled, '/', gprotein, "_", pca_type, "_", layer , ".tsv", sep="")
  }
  else if (pca_type == "All" && taste == TRUE  && dimr == "tSNE"){
    gproteinPath <- paste("embedding_all", scaled, '/', gprotein, "_", pca_type, "_", layer, '_tSNE', ".tsv", sep="")
    print ('tSNE')
  }
  else if (pca_type == "All" && taste == FALSE  && dimr == "tSNE"){
    gproteinPath <- paste("embedding_without_taste", scaled, '/', gprotein, "_", pca_type, "_", layer, '_tSNE', ".tsv", sep="")
    print ('tSNE')
  }
  else if (pca_type == "All" && taste == TRUE  && dimr == "UMAP"){
    gproteinPath <- paste("embedding_all", scaled, '/', gprotein, "_", pca_type, "_", layer, '_UMAP', ".tsv", sep="")
    print ('tSNE')
  }
  else if (pca_type == "All" && taste == FALSE  && dimr == "UMAP"){
    gproteinPath <- paste("embedding_without_taste", scaled, '/', gprotein, "_", pca_type, "_", layer, '_UMAP', ".tsv", sep="")
    print ('tSNE')
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
  if (TRUE) {
    var <- unique(data[[assay]])
    #var <- var[!is.na(var)]
    #var %>% replace(is.na(.), "Unknown")
    var[is.na(var)] <- "NA"
    #print (var)
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
  }
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
    
    observeEvent(input$gprotein, {
      ## Layers ML
      gproteinPath <- paste("layers.txt", sep="")
      data <- read_tsv(gproteinPath)
      subData <- subset(data, Gprotein == input$gprotein)
      print(subData$best_layer)
      bestLayer <- subData$best_layer
      print ('best layer')
      layers <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33')
      for (layer in layers) {
        if (layer == bestLayer){
          layers[as.numeric(bestLayer)+1] <- paste(bestLayer, ' (Best ML)')
        }
      }
      
      print (input$layer)
      print (input$option)
      if (input$option == 'gpcrome'){
        defaultLayer = layers[-1]
      }
      else if (input$option == 'best') {
        defaultLayer <- layers[as.numeric(bestLayer)+1]
      }
      else {
        if (is.null(input$layer) == TRUE){
          defaultLayer <- '0'
        }
        else {
          if (grepl("Best ML", input$layer) == TRUE) {
            defaultLayer = strsplit(layer, split=" ")
            defaultlayer <- defaultLayer[[1]][1]
          }
          else {
            defaultLayer <- input$layer
          }
        }
      }
      
      #output$layer <- renderUI({
      #  selectInput("layer", "Layer",
      #              choices = layers,
      #              selected = defaultLayer
      #              #selected = layers[as.numeric(bestLayer)+1]
      #  )
      #})
      
      print ('clusters')
      ## Clusters
      gproteinPath <- paste("optimalClustering.tsv", sep="")
      data <- read_tsv(gproteinPath)
      print (data)
      #subData <- filter(data, Gprotein %in% c(input$gprotein))
      subData <- data
      subData <- subData[order(subData$WSS),]
      optimalCluster <- (subData$Best)[[1]]
      clusters <- c(1:10)
      for (cluster in clusters) {
        if (cluster == optimalCluster){
          clusters[as.numeric(optimalCluster)] <- paste(optimalCluster, ' (Best)')
        }
      }
      
      print (optimalCluster)
      print (clusters)
      
      output$numClusters <- renderUI({
        selectInput("numClusters", "Cluster count",
                    choices = clusters,
                    selected = paste(optimalCluster, ' (Best KM)')
                    #selected = layers[as.numeric(bestLayer)+1]
        )
      })
      
      ## Layers
      gproteinPath <- paste("optimalClustering.tsv", sep="")
      data <- read_tsv(gproteinPath)
      print (data)
      #subData <- filter(data, Gprotein %in% c(input$gprotein))
      subData <- data
      subData <- subData[order(subData$WSS),]
      optimalLayer <- (subData$Layer)[[1]]
      #layers <- c(0:33)
      for (layer in layers) {
        if (layer == optimalLayer){
          layers[as.numeric(optimalLayer)] <- paste(optimalLayer, ' (Best KM)')
        }
      }
      
      print (optimalLayer)
      print (layers)
      
      if (input$option == 'gpcrome'){
        defaultLayer = layers[-1]
      }
      else if (input$option == 'bestKM') {
        defaultLayer <- layers[as.numeric(optimalLayer)+1]
      }
      else {
        if (is.null(input$layer) == TRUE){
          defaultLayer <- '0'
        }
        else {
          if (grepl("Best KM)", input$layer) == TRUE) {
            defaultLayer = strsplit(layer, split=" ")
            defaultlayer <- defaultLayer[[1]][1]
          }
          else {
            defaultLayer <- input$layer
          }
        }
      }
      
      output$layer <- renderUI({
        selectInput("layer", "Layers",
                    choices = layers,
                    selected = defaultLayer
                    #selected = paste(optimalLayer, ' (Best)')
                    #selected = layers[as.numeric(bestLayer)+1]
        )
      })
      
    })
      
    observeEvent(input$inputhelp, {
      showModal(modalDialog(
        title = "Help",
        HTML("This will work only when you select Color by as <kbd>Family</kbd>"
        )
      ))
    })
    
    print ('We are here')
    gproteinPath <- paste("all_layers/", "GNAS", "_", "All_0", ".tsv", sep="")
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
    observeEvent(c(input$assay, input$gprotein, input$pca_type, input$layer, input$taste, input$scaled, input$scaled, input$numClusters), {
      #output$scatter <- renderPlotly({plot(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$taste, input$numClusters, input$scaled)})
      X <- plot(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$layer, input$taste, input$numClusters, input$scaled, 'PCA')
      output$scatter <- renderPlotly({X$plot})
      RV$data <- X$df
      
      Y <- plot(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$layer, input$taste, input$numClusters, input$scaled, 'tSNE')
      output$scatter2 <- renderPlotly({Y$plot})
      RV$data2 <- Y$df
      
      Z <- plot(input$gprotein, input$assay, input$family, input$all, input$pca_type, input$layer, input$taste, input$numClusters, input$scaled, 'UMAP')
      output$scatter3 <- renderPlotly({Z$plot})
      RV$data3 <- Z$df
    })
    
    #output$scatter2 <- renderPlotly(plottSNE(input$gprotein, inmput$assay))
    
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
    
    output$table2 <- DT::renderDataTable({
      DT::datatable(RV$data2,
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
    
    output$table3 <- DT::renderDataTable({
      DT::datatable(RV$data3,
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
