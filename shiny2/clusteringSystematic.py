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
    data = []
    for embedding in ['shed', 'ebret']:
        for assay in ['STRING', 'IUPHAR', 'ebBRET', 'Shedding']:

            if assay == 'Shedding':
                gproteins = Shedding
            elif assay == 'ebBRET':
                gproteins = ebBRET
            elif assay == 'IUPHAR':
                gproteins = IUPHAR
            elif assay == 'STRING':
                gproteins = STRING

            #gproteins = ['GNAS', 'GNA12']

            for gprotein in gproteins:
                for layer in range(0,34):
                    if dimr == 'PCA':
                        df = pd.read_csv('pca_'+str(embedding)+'_scaled/'+gprotein+'_All_'+str(layer)+'_'+dimr+'.tsv', sep='\t')
                    else:
                        df = pd.read_csv('pca_'+str(embedding)+'_scaled/'+gprotein+'_All_'+str(layer)+'_'+dimr+'.tsv', sep='\t')
                    df[assay].replace({"Coupled": 1, "Not-coupled": 0}, inplace=True)
                    df = df[df[assay] != '-']
                    #print(df[['PC1', 'PC2']].to_numpy())
                    trueLabels = df[assay].to_numpy()
                    #print (trueLabels)
                    extract = df.columns[df.columns.str.contains('PC*')]
                    header = df.columns.tolist()
                    extract = []
                    for i in range(0, len(header)):
                        head = header[i]
                        if head[0:2] == 'PC':
                            extract.append(head)
                    #print (extract)
                    #X = df[['PC1', 'PC2']].to_numpy()
                    X = df[extract].to_numpy()
                    #print (X)
                    kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
                    #print (kmeans.labels_)
                    predLabels = kmeans.labels_
                    score = normalized_mutual_info_score(trueLabels, predLabels)
                    print (dimr + '\t' + embedding + '\t' + gprotein + '\t' + str(layer) + '\t' + assay + '\t' + str(round(score,3)))

                    row = []
                    row.append(dimr)
                    row.append(embedding)
                    row.append(assay)
                    row.append(gprotein)
                    row.append(layer)
                    row.append(round(score,3))
                    data.append(row)
            #sys.exit()

    df = pd.DataFrame(data, columns=['DimR', 'Embedding', 'Assay', 'G-protein', 'Layer', 'Score'])
    df = df.sort_values(by=['Embedding', 'Assay', 'G-protein', 'Score'], ascending=[False, False, False, False])
    df.to_csv('systematic'+dimr+'.tsv', sep='\t', index=False)


## For classess
#for dimr in ['UMAP', 'tSNE', 'PCA']:
for dimr in ['PCA']:
    data = []
    for embedding in ['shed', 'ebret']:
        assay = "Class"
        for layer in range(0,34):
            if dimr == 'PCA':
                df = pd.read_csv('pca_'+str(embedding)+'_scaled/GNAS_All_'+str(layer)+'_'+dimr+'.tsv', sep='\t')
            else:
                df = pd.read_csv('pca_'+str(embedding)+'_scaled/GNAS_All_'+str(layer)+'_'+dimr+'.tsv', sep='\t')
            df[assay].replace({"classA":0, "classB":1, "classC":2, "Frizzeled":3, "Taste":4}, inplace=True)
            df = df[df[assay] != 'Other']
            #print(df[['PC1', 'PC2']].to_numpy())
            trueLabels = df[assay].to_numpy()
            #print (trueLabels)
            extract = df.columns[df.columns.str.contains('PC*')]
            header = df.columns.tolist()
            extract = []
            for i in range(0, len(header)):
                head = header[i]
                if head[0:2] == 'PC':
                    extract.append(head)
            #print (extract)
            #X = df[['PC1', 'PC2']].to_numpy()
            X = df[extract].to_numpy()
            #print (X)
            kmeans = KMeans(n_clusters=5, random_state=0).fit(X)
            #print (kmeans.labels_)
            predLabels = kmeans.labels_
            score = normalized_mutual_info_score(trueLabels, predLabels)
            print (str(layer) + '\t' + str(round(score,3)))

            row = []
            row.append(dimr)
            row.append(embedding)
            row.append(layer)
            row.append(round(score,3))
            data.append(row)

    df = pd.DataFrame(data, columns=['DimR', 'Embedding', 'Layer', 'Score'])
    df = df.sort_values(by=['Embedding', 'Score'], ascending=[False, False])
    df.to_csv('systematicClass'+dimr+'.tsv', sep='\t', index=False)


## For families
'''
df = pd.read_csv('embedding_all/GNAS_All_0_tSNE.tsv', sep='\t')
df = df.dropna(subset=['Family'])
families = list(set(df['Family'].to_numpy()))

dfAllFamilies = pd.DataFrame(columns=['DimR', 'Family', 'Layer', 'Score'])
for family in families:
    for dimr in ['UMAP', 'tSNE', 'PCA']:
        data = []
        for layer in range(0,34):
            if dimr == 'PCA':
                df = pd.read_csv('all_layers_scaled/GNAS_All_'+str(layer)+'.tsv', sep='\t')
            else:
                df = pd.read_csv('embedding_all/GNAS_All_'+str(layer)+'_'+dimr+'.tsv', sep='\t')

            #df = df[df['Family'] != NaN]
            df = df.dropna(subset=['Family'])
            #print (df.to_numpy())
            df['FamilyValue'] = np.where(df['Family']== family, 1, 0)
            #print (df[''])
            #sys.exit()
            assay = "FamilyValue"
            #df[assay].replace({"classA":0, "classB":1, "classC":2, "Frizzeled":3, "Taste":4}, inplace=True)
            #df = df[df[assay] != 'Other']
            #print(df[['PC1', 'PC2']].to_numpy())
            trueLabels = df[assay].to_numpy()
            #print (trueLabels)
            X = df[['PC1', 'PC2']].to_numpy()
            #print (X)
            kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
            #print (kmeans.labels_)
            predLabels = kmeans.labels_
            score = normalized_mutual_info_score(trueLabels, predLabels)
            print (str(family) + '\t' + str(layer) + '\t' + str(round(score,3)))

            row = []
            row.append(dimr)
            row.append(family)
            row.append(layer)
            row.append(round(score,3))
            data.append(row)

        df = pd.DataFrame(data, columns=['DimR', 'Family', 'Layer', 'Score'])
        df = df.sort_values(by=['Score'], ascending=[False])
        df.to_csv('systematicFamily/'+family.replace(' ', '_').replace('>', '_').replace('<', '_').replace('/', '_')+'_'+dimr+'.tsv', sep='\t', index=False)

        dfAllFamilies = pd.concat([dfAllFamilies, df])
        #sys.exit()

dfAllFamilies.to_csv('systematicFamily/AllFamilies.tsv', sep='\t', index=False)
'''
