#!/home/marin/.conda/envs/transformer/bin/python
# coding: utf-8

import argparse
import gzip
import os, sys
import random
from collections import Counter
from tqdm import tqdm
import torch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import esm
import scipy
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
from sklearn.preprocessing import MinMaxScaler
import numpy as np
from joblib import dump, load

def main(d, uniq_id, gprotein, input_fasta,input_embedding,model,feature_type,embedding):
    train_set = model.split('/')[-1].split('_')[3]
    EMB_LAYER = int(model.split('/')[-1].split('_')[2])
    num_pca = float(model.split('/')[-1].split('_')[1])
    if num_pca > 1.0:
        num_pca = int(num_pca)

    ## TRAIN SET PATH
    if train_set == 'shed':
        TRAIN_FASTA_PATH = "/data/Users/marin/transformers/fasta/" + gprotein + "-shed.fa" # OR -shed.fa or IUPHAR.fa or -ebret.fa
    else:
        TRAIN_FASTA_PATH = "/data/Users/marin/transformers/fasta/" + gprotein + "-ebret1.fa"

    if train_set == 'shed':
        TRAIN_EMB_PATH = "/data/Users/marin/transformers/fasta/" + gprotein + "-shed_"+embedding+"_all/" # Path to embedings if you want the esm1 remove b
    else:
        TRAIN_EMB_PATH = "/data/Users/marin/transformers/fasta/" + gprotein + "-ebret_"+embedding+"_all/" # Path to embedings if you want the esm1 remove b

    ## TEST SET PATH
    g_fam = {'GNAS': 'Gs', 'GNAL': 'Gs', 'Gs': 'Gs', 'Gi1':'Gio', 'Gi2': 'Gio', 'Gz': 'Gio', 'GoA': 'Gio', 'GoB': 'Gio', 'GNAI1': 'Gio', 'GNAI3': 'Gio','GNAZ': 'Gio', 'GNAO1': 'Gio','G12':'G1213', 'G13': 'G1213', 'GNA12': 'G1213', 'GNA13': 'G1213','Gq':'Gq11', 'G11':'Gq11', 'G14':'Gq11','G15':'Gq11', 'GNAQ': 'Gq11', 'GNA14': 'Gq11','GNA15': 'Gq11'}

    TEST_FASTA_PATH = input_fasta
    TEST_EMB_PATH = input_embedding

    ## TEST SET
    Xtest = []
    gpcr_test = []
    for header, _seq in esm.data.read_fasta(TEST_FASTA_PATH):
        if header.split('>')[1]:
            gpcr_test.append(header.split('>')[1])
            fn = TEST_EMB_PATH+header[1:]+'.pt'
            embs = torch.load(fn)
            Xtest.append(embs['mean_representations'][EMB_LAYER])
            #if len(Xtest) == 50:
             #   break
    Xtest = torch.stack(Xtest, dim=0).numpy()

    ###Xs_train = Xs
    Xs_test = Xtest
    ###ys_train = np.array(ys)

    scaler = load('static/predictor/scaler/all/scaler_'+gprotein+'_'+str(num_pca)+'_'+str(EMB_LAYER)+'_'+str(train_set)+'_'+embedding)
    pca = load('static/predictor/pca/all/pca_'+gprotein+'_'+str(num_pca)+'_'+str(EMB_LAYER)+'_'+str(train_set)+'_'+embedding)

    ###pca = PCA(num_pca) #204 iuphar,118 shedding ,67 ebret,bret new 78
    ###pca.fit(Xs_train)
    ###Xs_train_pca = pca.transform(Xs_train)
    Xs_test_pca = pca.transform(Xs_test)

    ## Load handcrafted features
    if feature_type > 0:
        if train_set == 'shed':
            if feature_type == 1:
                path = 'feature_files/shedding_features/seq/'
            else:
                path = 'feature_files/shedding_features/struc/'
        else:
            if feature_type == 1:
                path = 'feature_files/ebbret_features/seq/'
            else:
                path = 'feature_files/ebbret_features/struc/'

        class handcrafted_features:
            def __init__(self, name, values):
                self.name = name
                self.values = values

        features = {}
        for files in [path+gprotein+'_train.txt', 'output/'+uniq_id+'/'+train_set+'/seq_features/'+gprotein+'.txt']:
            for line in open(files, 'r'):
                if line.split('\t')[0] != 'GPCR':
                    name = line.split('\t')[0]
                    if 'train' in files:
                        values = np.array(line.split('\t')[1:-1]).astype(float)
                    else:
                        values = np.array(line.split('\t')[1:]).astype(float)
                    features[name] = handcrafted_features(name, values)
                    #if len(values) < 24:
                     #   print (name, train_set)
        '''
        X_seq = []
        X_pca = []
        for name, pca_values in zip(gpcr_train, Xs_train_pca):
            if name in features:
                X_seq.append(features[name].values)
                X_pca.append(pca_values)

        scaler = MinMaxScaler()

        Xs_train_pca = np.concatenate((np.array(X_pca), np.array(X_seq)), axis=1)
        #
        scaler.fit(Xs_train_pca)
        Xs_train_pca = scaler.transform(Xs_train_pca)
        #
        '''

        X_seq = []
        X_pca = []
        gpcr_test_new = []
        for name, pca_values in zip(gpcr_test, Xs_test_pca):
            if name in features:
                gpcr_test_new.append(name)
                X_seq.append(features[name].values)
                X_pca.append(pca_values)

        gpcr_test = gpcr_test_new
        #X_seq = scaler.transform(X_seq)
        Xs_test_pca = np.concatenate((X_pca, X_seq), axis=1)
        #
        Xs_test_pca = scaler.transform(Xs_test_pca)
        #
    else:
        ###scaler = MinMaxScaler()
        #
        ###scaler.fit(Xs_train_pca)
        ###Xs_train_pca = scaler.transform(Xs_train_pca)
        Xs_test_pca = scaler.transform(Xs_test_pca)
        #

    ###dump(scaler, 'scaler/scaler_'+gprotein+'_'+str(num_pca)+'_'+str(EMB_LAYER)+'_'+str(train_set)+'_'+embedding)
    ###dump(pca, 'pca/pca_'+gprotein+'_'+str(num_pca)+'_'+str(EMB_LAYER)+'_'+str(train_set)+'_'+embedding)
    Xs_test_pca_copy = Xs_test_pca[:]

    model_name = load(model)
    try:
        ys_pred = model_name.predict_proba(Xs_test_pca)
    except:
        ys_pred = model_name.predict(Xs_test_pca)

    for name, y in zip(gpcr_test, ys_pred):
        if name not in d:
            d[name] = {}
        d[name][gprotein] = y[1]

    return (Xs_test_pca_copy)
