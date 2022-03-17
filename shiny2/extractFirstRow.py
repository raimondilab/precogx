import os, sys
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import normalized_mutual_info_score
import numpy as np

## For Shedding, ebBRET, IUPHAR and STRING
STRING = ['Barr1-GRK2', 'Barr2', 'Barr2-GRK2']
ebBRET = ['GNAS', 'GNAI1', 'GNAI2', 'GoA', 'GoB', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA11', 'GNA14', 'GNA15', 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2']
Shedding = ['GNAS', 'GNAL', 'GNAI1', 'GNAI3', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15']
Shedding.sort()
ebBRET.sort()
IUPHAR = list(set(ebBRET+Shedding))
IUPHAR.remove('Barr2')
IUPHAR.remove('Barr2-GRK2')
IUPHAR.remove('Barr1-GRK2')
IUPHAR.sort()

#for dimr in ['UMAP', 'tSNE', 'PCA']:
for dimr in ['PCA']:
    df = pd.read_csv('systematic'+dimr+'.tsv', sep='\t')
    for embedding in ['shed', 'ebret']:
        dfEmbed = df[df['Embedding'] == embedding]
        #print (dfEmbed)

        for assay in ['STRING', 'IUPHAR', 'ebBRET', 'Shedding']:
            if assay == 'Shedding':
                gproteins = Shedding
            elif assay == 'ebBRET':
                gproteins = ebBRET
            elif assay == 'IUPHAR':
                gproteins = IUPHAR
            elif assay == 'STRING':
                gproteins = STRING

            dfAssay = dfEmbed[dfEmbed['Assay'] == assay]
            #print (dfAssay)

            for gprotein in gproteins:
                dfGprotein = dfAssay[dfAssay['G-protein'] == gprotein]
                firstRow = dfGprotein.to_numpy()[0].tolist()
                firstRow = [str(value) for value in firstRow]
                print ('\t'.join(firstRow))
