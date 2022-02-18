library('mclust')
library('tidyverse')
data(diabetes)
print (diabetes)
class <- diabetes$class
table(class)
## class
## Chemical   Normal    Overt 
##       36       76       33
X <- diabetes[,-1]
head(X)
##   glucose insulin sspg
## 1      80     356  124
## 2      97     289  117
## 3     105     319  143
## 4      90     356  199
## 5      90     323  240
## 6      86     381  157
clPairs(X, class)

d_clust <- Mclust(as.matrix(X), G=1:20)
dim(d_clust$z)
m.best <- dim(d_clust$z)[2]
cat("model-based optimal number of clusters:", m.best, "\n")
# 4 clusters
plot(d_clust)


scaled <- '_scaled'
pca_type <- 'All'
values <- data.frame('Gprotein' = character(), 'Layer' = integer(), 'Best' = integer(), 'WSS' = double())
#gproteins <- c('GNAI3', 'GNAI2', 'GNA14', 'GNA12', 'GoA', 'Barr2-GRK2', 'GoA', 'GNA15', 'Barr2', 'Barr1-GRK2', 'GNAQ', 'GNAO1', 'GNAI1', 'GNAS', 'GNAZ', 'GNA11', 'GNA13', 'GNAL')
#gproteins <- c('GNAI3', 'GNAI2', 'GNA14')
#gproteins <- c('GNAS', 'GNAZ', 'GNA11', 'GNA13', 'GNAL')
gproteins <- c('GNA12')
for (gprotein in gproteins) {
  for (layer in 0:33) {
    gproteinPath <- paste("all_layers", scaled, '/', gprotein, "_", pca_type, "_", layer , ".tsv", sep="")
    
    print (gproteinPath)
    data <- read_tsv(gproteinPath)
    data <- data[c('PC1', 'PC2')]
    
    d_clust <- Mclust(as.matrix(data), G=1:10)
    dim(d_clust$z)
    optimalClusters <- dim(d_clust$z)[2]
    
    optimalClusters
    #k.max <- 5
    
    wss <- sapply(optimalClusters:optimalClusters, function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
    #wss <- sapply(1:optimalClusters, function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
    #print (wss)
    #print (optimalClusters)
    values[nrow(values) + 1,] = c(gprotein, layer, optimalClusters, wss)
    
    #plot(1:optimalClusters, wss,
    #     type="b", pch = 19, frame = FALSE, 
    #     xlab="Number of clusters K",
    #     ylab="Total within-clusters sum of squares")
  }
}
values
values$WSS
plot(values$Layer, values$WSS,
          type="b", pch = 19, frame = FALSE, 
          xlab="Number of Layers",
          ylab="WSS scores")

write.table(values, file='optimalClustering.tsv', quote = FALSE, row.names=FALSE, sep='\t')

