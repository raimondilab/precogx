library(tsne)
library(tidyverse)

pathToFolder <- "embedding_all"
#pathToFolder <- "embedding_without_taste"
#files <- list.files(path=pathToFolder, pattern="*_embedding.tsv", full.names=TRUE, recursive=FALSE)
for (num in 0:33){
  x <- paste(pathToFolder, '/', num, '_embedding.tsv', sep='')
  print (x)
  data <- read_tsv(x)
  tsne <- tsne(data)
  tsne <- data.frame(tsne)
  #print(tsne)
  name <- strsplit(x, split="/")
  name <- name[[1]][2]
  name <- strsplit(name, split="_")
  name <- name[[1]][1]
  print (name)
  
  write.table(tsne, file=paste(pathToFolder, '/', name, '_tSNE.tsv', sep=''), quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE)
}