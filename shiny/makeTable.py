import os, sys, json
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from Bio import SearchIO
from Bio.Blast import NCBIXML
from matplotlib import cm
import pandas as pd

class GPCR:
    def __init__(self, gene, acc):
        self.gene = gene
        self.acc = acc
        self.assay = 'Shedding'
        self.iuphar = '-'
        self.shedding = '-'
        self.ebbret = '-'
        self.string = '-'
        self.family = 'NA'
        self.cls = 'NA'
        self.x = None
        self.y = None

path = os.getcwd()

def extract_pca(gprotein, pca_type, taste):
    if pca_type == 'Best' and taste == 1:
        Xs_train_pca = np.load(path+'/../static/best_PCA/'+gprotein+'.npy', allow_pickle=True)
    elif pca_type == 'GPCRome' and taste == 1:
        #Xs_train_pca = np.load(path+'/static/33layer_PCA/33layer.npy', allow_pickle=True)
        Xs_train_pca = np.load(path+'/../static/best_PCA/GNAZ.npy', allow_pickle=True)
    elif pca_type == 'Best' and taste == 0:
        Xs_train_pca = np.load('best_pca_without_taste/'+gprotein+'.npy', allow_pickle=True)
    elif pca_type == 'GPCRome' and taste == 0:
        #Xs_train_pca = np.load(path+'/static/33layer_PCA/33layer.npy', allow_pickle=True)
        Xs_train_pca = np.load('33layer_PCA_without_taste/33layer.npy', allow_pickle=True)
    #score_coupling, score_uncoupling, Xs_train_pca_coupling, Xs_train_pca_uncoupling, Xs_train_pca_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = filter_gpcr_list(Xs_train_pca, assay, gprotein)
    filter_gpcr_list(Xs_train_pca, gprotein, pca_type, taste)
    '''
    #print ('train', Xs_train_pca_coupling)
    score_coupling = score_coupling.tolist()
    score_uncoupling = score_uncoupling.tolist()
    x_train_coupling = Xs_train_pca_coupling[:,0].tolist()
    x_train_uncoupling = Xs_train_pca_uncoupling[:,0].tolist()
    x_train_grey = Xs_train_pca_grey[:,0].tolist()
    y_train_coupling = Xs_train_pca_coupling[:,1].tolist()
    y_train_uncoupling = Xs_train_pca_uncoupling[:,1].tolist()
    y_train_grey = Xs_train_pca_grey[:,1].tolist()
    print (len(x_train_coupling))
    print (len(y_train_coupling))
    return score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey
    '''

def filter_gpcr_list(X, gprotein, pca_type, taste):
    gpcr_list = [];
    dic = {};
    if taste == 1:
        for line in open(path+'/../static/gpcr_list_new_GN.txt', 'r'):
            gene = line.replace('\n', '').split('\t')[1]
            acc = line.replace('\n', '').split('\t')[0]
            gpcr_list.append(gene)
            #gpcr_list.append(line.replace('\n', '').split('\t')[1])
            if gene not in dic:
                dic[gene] = GPCR(gene, acc)
    else:
        print (path)
        for line in open(path+'/../static/gpcr_list_new_without_TASTE_GN.txt', 'r'):
            gene = line.replace('\n', '').split('\t')[1]
            acc = line.replace('\n', '').split('\t')[0]
            gpcr_list.append(gene)
            #gpcr_list.append(line.replace('\n', '').split('\t')[1])
            if gene not in dic:
                dic[gene] = GPCR(gene, acc)

    coupledGenes = []
    notcoupledGenes = []
    #print (assay)
    #assay = 'ebBRET'

    ## Shedding
    num = -1
    for line in open(path+'/../data/shedding.tsv', 'r'):
        if line[0] != '#' and num!= -1:
            gene = line.split('\t')[0]
            acc = line.split('\t')[1]
            score = float(line.split('\t')[num+2])

            if gene in dic:
                if score >= -1.0:
                    coupledGenes.append(gene+'|'+acc)
                    dic[gene].shedding = 'Coupled'
                else:
                    dic[gene].shedding = 'Not-coupled'

            print (gene, score, dic[gene].shedding)

        else:
            flag = 0
            header = line.replace('\n', '').split('\t')[2:]
            for num, gprot in enumerate(header):
                if gprot == gprotein:
                    flag = 1
                    break
            if flag == 0:
                num = -1

    ## ebBRET
    num = -1
    for line in open(path+'/../data/ebbret.tsv', 'r', encoding="utf-8"):
        if line[0] != '#':
            gene = line.split('\t')[0]
            acc = line.split('\t')[1]
            score = float(line.split('\t')[num+2])

            if gene in dic:
                if score > 0:
                    coupledGenes.append(gene+'|'+acc)
                    dic[gene].ebbret = 'Coupled'
                else:
                    notcoupledGenes.append(gene+'|'+acc)
                    dic[gene].ebbret = 'Not-coupled'
        else:
            header = line.replace('\n', '').split('\t')[2:]
            for num, gprot in enumerate(header):
                if gprot == gprotein:
                    #print (gprot)
                    break

    ## IUPHAR
    iuphar_map = {
                  'GNAS': 'Gs', 'GNAL': 'Gs',
                  'GNAI1': 'Gi/Go', 'GNAI2': 'Gi/Go', 'GNAI3': 'Gi/Go', 'GNAO1': 'Gi/Go', 'GNAZ': 'Gi/Go', 'GoA': 'Gi/Go', 'GoB': 'Gi/Go',
                  'GNA12': 'G12/G13', 'GNA13': 'G12/G13',
                  'GNAQ': 'Gq/G11', 'GNA11': 'Gq/G11', 'GNA14': 'Gq/G11', 'GNA15': 'Gq/G11'
                  }
    if gprotein in iuphar_map:
        gprotein_fam = iuphar_map[gprotein]
        for line in open(path+'/../data/iuphar.tsv', 'r'):
            if line[0] != '#' and line.split('\t')[1] != '':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                if gene in dic:
                    if gprotein_fam in line:
                        coupledGenes.append(gene+'|'+acc)
                        dic[gene].iuphar = 'Coupled'
                    else:
                        notcoupledGenes.append(gene+'|'+acc)
                        dic[gene].iuphar = 'Not-coupled'

    ## STRING
    string_map = {
                  'Barr1-GRK2': 'ARRB1',
                  'Barr2': 'ARRB2',
                  'Barr2-GRK2': 'ARRB2'
                  }
    if gprotein in string_map:
        barr = string_map[gprotein]
        num = -1
        for line in open(path+'/../data/string.tsv', 'r', encoding="utf-8"):
            gene = line.split('\t')[0]
            acc = line.split('\t')[1]
            if gene in dic:
                if barr in line:
                    coupledGenes.append(gene+'|'+acc)
                    dic[gene].string = 'Coupled'
                else:
                    notcoupledGenes.append(gene+'|'+acc)
                    dic[gene].string = 'Not-coupled'
    #print (genes_to_consider_coupling)
    #print (genes_to_consider_uncoupling)

    for line in open(path+'/../static/predictor/data_precog2/IUPHAR_couplings.tsv', 'rt'):
        gene = line.split('\t')[1]
        family = line.split('\t')[3]
        if gene in dic:
            dic[gene].family = family

    for line in open(path+'/../data/classification.txt', 'rt'):
        if 'Uniprot_acc' not in line and line[0] != '\n':
            gene = line.split('\t')[1]
            #print (line.split('\t'))
            cls = line.split('\t')[3].replace('\n', '')
            if gene in dic:
                dic[gene].cls = cls

    data = []
    for gene, vector in zip(gpcr_list, X):
        dic[gene].x = vector[0]
        dic[gene].y = vector[1]
        row = []
        row.append(gene)
        row.append(dic[gene].iuphar)
        row.append(dic[gene].shedding)
        row.append(dic[gene].ebbret)
        row.append(dic[gene].string)
        row.append(dic[gene].family)
        row.append(dic[gene].cls)
        row.append(dic[gene].x)
        row.append(dic[gene].y)
        data.append(row)

    #print (data)
    df = pd.DataFrame(data, columns = ['Gene', 'IUPHAR', 'Shedding', 'ebBRET', 'STRING', 'Family', 'Class', 'PC1', 'PC2'])
    print (df)
    df.to_csv(gprotein+'_'+pca_type+'.tsv', sep='\t', index=False)

ebBRET = ['GNAS', 'GNAI1', 'GNAI2', 'GoA', 'GoB', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA11', 'GNA14', 'GNA15', 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2'];
shedding = ['GNAS', 'GNAL', 'GNAI1', 'GNAI3', 'GNAO1', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15'];
interactors = list(set(ebBRET + shedding))

for interactor in interactors:
    print (interactor)
    #extract_pca(interactor, 'GPCRome')
    if interactor != 'GNAO1':
        extract_pca(interactor, 'GPCRome', 0)

print (interactors)
