#!/usr/bin/env python3
# coding: utf-8

import argparse, requests, urllib
import gzip
import os, sys, json
import random
from collections import Counter
from tqdm import tqdm
import torch
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import esm
import shutil
import scipy
from sklearn.model_selection import GridSearchCV, train_test_split
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold
from sklearn.neighbors import KNeighborsClassifier, KNeighborsRegressor
from sklearn.svm import SVC, SVR
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression, SGDRegressor
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.metrics import f1_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import matthews_corrcoef
from sklearn.preprocessing import MinMaxScaler
from xgboost import XGBClassifier, plot_importance
from xgboost import plot_tree, to_graphviz
import numpy as np
from joblib import dump, load
import time, random, string
from Bio import SeqIO
import predict
import extract

class GPCR:
    def __init__(self, name):
        self.name = name
        self.seq = ''
        self.var = []

def main(input, input_file, assay):
    while True:
        uniq_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k = 5))
        if os.path.exists('output/'+uniq_id) == False:
            os.system('mkdir static/predictor/output/'+uniq_id)
            os.system('mkdir static/predictor/output/'+uniq_id+'/shed/')
            os.system('mkdir static/predictor/output/'+uniq_id+'/shed/seq_features/')
            os.system('mkdir static/predictor/output/'+uniq_id+'/ebret/')
            os.system('mkdir static/predictor/output/'+uniq_id+'/ebret/seq_features/')
            os.system('mkdir static/predictor/output/'+uniq_id+'/embed/')
            input_embedding = 'static/predictor/output/'+uniq_id+'/embed/'
            break

    print ('Your output will be stored at: static/predictor/output/'+uniq_id)
    if input_file == None:
        input = formatInput(input)
    else:
        input = formatInputFile(input_file)

    open('static/predictor/output/'+uniq_id+'/input.fasta', 'w').write(input)
    input_file = 'static/predictor/output/'+uniq_id+'/input.fasta'

    '''
    record = SeqIO.read(input_file, "fasta")

    if len(record)> 1024:
        print('Sequence too long to create embeddings-max 1024 aa')
        sys.exit()
    '''
    ###### Input models ######

    if assay == 'all':
        #path_to_model = '/data/Users/marin/transformers/runs/final_models/best_all/'
        path_to_model = 'static/predictor/best_all/'
        data = [['GNAS_0.95_28_shed_esm1b.joblib', 0],['GNAL_0.95_33_shed_esm1b.joblib', 0],['GNAI1_0.95_31_ebret_esm1b.joblib', 0],
            ['GNAI2_0.95_20_ebret_esm1b.joblib', 0],['GNAI3_0.95_29_shed_esm1b.joblib',0],
            ['GoA_0.95_20_ebret_esm1b.joblib', 0],['GoB_0.95_33_ebret_esm1b.joblib', 0],['GNAZ_0.95_32_ebret_esm1b.joblib', 0],['GNA11_0.95_25_ebret_esm1b.joblib', 0],
            ['GNA14_0.95_18_shed_esm1b.joblib', 0],['GNA15_0.95_29_shed_esm1b.joblib', 0],['GNAQ_0.95_31_shed_esm1b.joblib', 0],['GNA12_0.95_18_shed_esm1b.joblib',0],
            ['GNA13_0.95_18_shed_esm1b.joblib',0],['Barr1-GRK2_0.95_0_ebret_esm1b.joblib',0],
            ['Barr2_0.95_33_ebret_esm1b.joblib',0],['Barr2-GRK2_0.95_0_ebret_esm1b.joblib',0]]

    # run hmmsearch
    os.system('hmmsearch data/7tm_1.hmm '+input_file+' > static/predictor/output/'+uniq_id+'/temp_hmm_file.txt')

     # create embeddings
    for row in data:
        if "esm1b"  in  row[0].split('.')[1].split('_')[-1]:
            model_location= "esm1b_t33_650M_UR50S"

    md=[]
    for row in data:
        md=[int(row[0].split('.')[1].split('_')[1]) for row in data]
        repr_layer=list(set(md))
        repr_layer.sort()

    extract.main(model_location,input_file,str(input_embedding),repr_layer)

    #print ('Generating hand-creafted features')
    #os.system('./precog.py '+input_file+' '+str(uniq_id))
    #os.system('./precog2.py '+input_file+' '+str(uniq_id))

    print ('Generating embeddings')
    d = {}; gproteins=[]
    for row in data:
        gprotein = row[0].split('_')[0]
        ###print (gprotein)
        gproteins.append(gprotein)
        model = path_to_model + row[0]
        feature_type = row[1]
        embedding = row[0].split('.')[1].split('_')[-1]

        Xs_test_pca_copy = predict.main(d, uniq_id, gprotein, input_file, input_embedding, model, int(feature_type), str(embedding))
        np.save('static/predictor/output/'+uniq_id+'/'+gprotein, Xs_test_pca_copy)

    name=str(input_file)
    #name1=name.split('/')[7].split('.')[0]
    l = '#PRECOGx\n'
    l += '#Input'
    for gprotein in gproteins:
        l += '\t' + str(gprotein)
    l += '\n'
    for name in d:
        l += name
        for gprotein in gproteins:
            l += '\t' + str(round(d[name][gprotein],3))
        l += '\n'
    #print (l)
    #shutil.rmtree('output/'+uniq_id+'/')
    open('static/predictor/output/'+uniq_id+'/out.tsv', 'w').write(l)
    print ('The output is saved at static/predictor/output/'+uniq_id)

    data = []
    dic = {}
    for line in open('static/predictor/output/'+uniq_id+'/out.tsv', 'r'):
        row = []
        if line[0] != '#':
            row = line.replace('\n', '').split('\t')
            data.append(row)

    dic = {'data': data}
    with open('static/predictor/output/'+uniq_id+'/out.json', 'w') as f:
        json.dump(dic, f)

    #shutil.rmtree('output/'+uniq_id+'/')
    return (uniq_id)

def formatInputFile(input_file):
    ## if input in FASTA format
    gpcrs = {}
    for line in open(input_file, 'r'):
        if line[0] == '>':
            given_name = line.split('>')[1].replace('\n', '').replace(' ','')
            if '/' in given_name:
                name = given_name.split('/')[0]
                variant = given_name.split('/')[1]
            else:
                name = given_name
                variant = 'wt'

            if name not in gpcrs:
                gpcrs[name] = GPCR(name)

            gpcrs[name].var.append(variant)
        else:
            gpcrs[name].seq += line.replace('\n', '')

    new_input = ''
    for name in gpcrs:
        for variant in gpcrs[name].var:
            if variant == 'WT':
                new_input += '>'+name+'_WT\n'
                new_input += gpcrs[name].seq + '\n'
            else:
                new_input += '>'+name+'_'+variant+'\n'
                position = int(variant[1:-1]) - 1
                newAA = variant[-1]
                new_input += gpcrs[name].seq[:position] + newAA + gpcrs[name].seq[position:] + '\n'
    print (new_input)
    #sys.exit()
    return (new_input)

def formatInput(input):
    ## if input in FASTA format
    gpcrs = {}
    if '>' in input:
        for line in input.split('\n'):
            if line[0] == '>':
                given_name = str(line.split('>')[1].replace('\n', '').replace(' ','').split()[0])
                if '/' not in given_name:
                    name = given_name
                    variant = 'WT'
                else:
                    name = given_name.split('/')[0]
                    variant = given_name.split('/')[1]
                if name not in gpcrs:
                    gpcrs[name] = GPCR(name)

                gpcrs[name].var.append(variant)
            else:
                gpcrs[name].seq += line.replace('\n', '').split()[0]
    else:
        ## if input is uniprot acc or gene name
        for line in input.split('\n'):
            given_name = str(line.replace('\n', '').replace(' ','').split()[0])
            if '/' not in given_name:
                name = given_name
                variant = 'WT'
            else:
                name = given_name.split('/')[0]
                variant = given_name.split('/')[1]

            if name not in gpcrs:
                #print (name, variant)
                gpcrs[name] = GPCR(name)
                gpcrs[name].seq = fetchSeq(name)

            gpcrs[name].var.append(variant)

    new_input = ''
    for name in gpcrs:
        if 'WT' not in gpcrs[name].var:
            gpcrs[name].var.append('WT')
        for variant in gpcrs[name].var:
            if gpcrs[name].seq != None:
                if variant == 'WT':
                    new_input += '>'+name+'_WT\n'
                    new_input += gpcrs[name].seq + '\n'
                else:
                    new_input += '>'+name+'_'+variant+'\n'
                    position = int(variant[1:-1]) - 1
                    newAA = variant[-1]
                    #print (position, newAA)
                    new_input += gpcrs[name].seq[:position] + newAA + gpcrs[name].seq[position+1:] + '\n'
    #print (new_input)
    #sys.exit()
    return (new_input)

def fetchSeq(name):
    #print ('here')
    response = requests.get('https://www.uniprot.org/uniprot/'+str(name))
    #print (response.status_code)
    if response.status_code == 200:
        response = urllib.request.urlopen('https://www.uniprot.org/uniprot/'+name+'.fasta')
        seq = ''
        for line in str(response.read().decode('utf-8')).split('\n'):
            if len(line.split())>0:
                if line[0] != '>':
                    seq += line.replace('\n', '')
    else:
        url = 'https://www.uniprot.org/uploadlists/'
        #print (name)
        params = {
        'from': 'GENENAME',
        'to': 'ACC',
        'format': 'tab',
        'taxon': '9606',
        'query': name
        }

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
           response = f.read()
        acc = None
        seq = None
        for line in response.decode('utf-8').split('\n'):
            if 'From' not in line:
                if len(line.split())>0:
                    acc = str(line.split('\t')[1].replace('\n', ''))
                    break
        #print (acc)
        if acc != None:
            response = requests.get('https://www.uniprot.org/uniprot/'+acc+'.fasta')
            if response.status_code == 200:
                response = urllib.request.urlopen('https://www.uniprot.org/uniprot/'+acc+'.fasta')
                seq = ''
                for line in response.read().decode('utf-8').split('\n'):
                    if len(line.split())>0:
                        if line[0] != '>':
                            seq += line.replace('\n', '')

    return seq

if __name__ == "__main__":
    ## Parser
    parser = argparse.ArgumentParser(description='Script to predict GPCR couplings using ESM and/or seq-based features', epilog='End of help')
    parser.add_argument('assay', help='Input what assay is used (ebret or shed)')
    parser.add_argument('--input_file', help='Input File (FASTA formatted)')
    parser.add_argument('--input', help='Input (FASTA formatted)')
    args = parser.parse_args()
    assay = args.assay
    input = args.input
    input_file = args.input_file
    uniq_id = main(input, input_file, assay)
