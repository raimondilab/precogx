library(umap)
library(tidyverse)

pathToFolder <- "pca_shed_scaled"
#pathToFolder <- "embedding_all"
#pathToFolder <- "embedding_without_taste"
#files <- list.files(path=pathToFolder, pattern="*_embedding.tsv", full.names=TRUE, recursive=FALSE)
for (num in 0:33){
  x <- paste(pathToFolder, '/', num, '_embedding.tsv', sep='')
  print (x)
  data <- read_tsv(x)
  umap <- umap(data, n_components = 2, random_state = 15)
  umap <- umap[["layout"]]
  print (umap)
  #umap <- data.frame(umap)
  #print(umap)
  name <- strsplit(x, split="/")
  name <- name[[1]][2]
  name <- strsplit(name, split="_")
  name <- name[[1]][1]
  print (name)
  
  write.table(umap, file=paste(pathToFolder, '/', name, '_UMAP.tsv', sep=''), quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE)
}
