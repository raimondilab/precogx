import numpy as np
from sklearn.decomposition import PCA
import pandas as pd

pathToFolder = "embedding_all"
#pathToFolder <- "embedding_without_taste"
#files <- list.files(path=pathToFolder, pattern="*_embedding.tsv", full.names=TRUE, recursive=FALSE)
for num in range(0, 34):
  data = pd.read_csv(pathToFolder+'/'+str(num)+'_embedding.tsv', sep='\t')
  pca = PCA(n_components=0.95, svd_solver="full")
  newData = pca.fit_transform(data)
  np.savetxt(pathToFolder+'/'+str(num)+'_PCA.tsv', newData, delimiter='\t')


  #write.table(tsne, file=paste(pathToFolder, '/', name, '_tSNE.tsv', sep=''), quote=FALSE, sep='\t', col.names = FALSE, row.names=FALSE)
