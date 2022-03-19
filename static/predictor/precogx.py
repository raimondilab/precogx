#!/usr/bin/env python3
# coding: utf-8

import argparse, requests, urllib
import gzip
import regex as re
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
#import precogxb_app.static.predictor.predict as predict
#import precogxb_app.static.predictor.extract as extract
#path = app.root_path
#sys.path.insert(1, path + '/static/predictor/')
#import extract
#import predict

class GPCR:
    def __init__(self, name):
        self.name = name
        self.seq = ''
        self.var = []
        self.shedding = {}
        self.iuphar = {}
        self.ebbret = {}

def main(numseqs, input, input_file, assay, path):
    #homeDir = os.getcwd()
    #homeDir = '/var/www/flask_apps/precogxb_app/'
    homeDir = path
    sys.path.insert(1, path + '/static/predictor/')
    import extract
    import predict
    #print ('hello', homeDir, os.getcwd(), os.listdir('.'), path)
    while True:
        uniq_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k = 5))
        if os.path.exists('output/'+uniq_id) == False:
            #os.system('mkdir static/predictor/output/'+uniq_id)
            os.system('mkdir ' + homeDir + '/static/predictor/output/'+uniq_id)
            os.system('mkdir ' + homeDir+ '/static/predictor/output/'+uniq_id+'/shed/')
            os.system('mkdir ' + homeDir + '/static/predictor/output/'+uniq_id+'/shed/seq_features/')
            os.system('mkdir ' + homeDir + '/static/predictor/output/'+uniq_id+'/ebret/')
            os.system('mkdir ' + homeDir +'/static/predictor/output/'+uniq_id+'/ebret/seq_features/')
            os.system('mkdir ' + homeDir + '/static/predictor/output/'+uniq_id+'/embed/')
            os.system('mkdir ' + homeDir + '/static/predictor/output/'+uniq_id+'/PCA/')
            os.system('mkdir ' + homeDir + '/static/predictor/output/'+uniq_id+'/attentions/')
            save_path = homeDir + '/static/predictor/output/'+uniq_id
            input_embedding = homeDir + '/static/predictor/output/' + uniq_id + '/embed/'
            '''
            os.system('mkdir ' + '/static/predictor/output/'+uniq_id)
            os.system('mkdir ' + '/static/predictor/output/'+uniq_id+'/shed/')
            os.system('mkdir ' + '/static/predictor/output/'+uniq_id+'/shed/seq_features/')
            os.system('mkdir ' + '/static/predictor/output/'+uniq_id+'/ebret/')
            os.system('mkdir ' + '/static/predictor/output/'+uniq_id+'/ebret/seq_features/')
            os.system('mkdir ' + '/static/predictor/output/'+uniq_id+'/embed/')
            os.system('mkdir ' + '/static/predictor/output/'+uniq_id+'/PCA/')
            #input_embedding = '/static/predictor/output/'+uniq_id+'/embed/'
            '''
            break

    print ('Your output will be stored at: static/predictor/output/'+uniq_id)
    if input_file == None:
        gpcrs, input = formatInput(numseqs, input)
    else:
        gpcrs, input = formatInputFile(numseqs, input_file)

    open(homeDir + '/static/predictor/output/'+uniq_id+'/input.fasta', 'w').write(input)
    input_file = homeDir + '/static/predictor/output/'+uniq_id+'/input.fasta'

    '''
    record = SeqIO.read(input_file, "fasta")

    if len(record)> 1024:
        print('Sequence too long to create embeddings-max 1024 aa')
        sys.exit()
    '''
    ###### Input models ######

    if assay == 'all':
        #path_to_model = '/data/Users/marin/transformers/runs/final_models/best_all/'
        path_to_model = homeDir + '/static/predictor/best_all/'
        data = [['GNAS_0.95_28_shed_esm1b.joblib', 0],['GNAL_0.95_33_shed_esm1b.joblib', 0],['GNAI1_0.95_31_ebret_esm1b.joblib', 0],
            ['GNAI2_0.95_20_ebret_esm1b.joblib', 0],['GNAI3_0.95_29_shed_esm1b.joblib',0],
            ['GoA_0.95_20_ebret_esm1b.joblib', 0],['GoB_0.95_33_ebret_esm1b.joblib', 0],['GNAZ_0.95_32_ebret_esm1b.joblib', 0],['GNA11_0.95_25_ebret_esm1b.joblib', 0],
            ['GNA14_0.95_18_shed_esm1b.joblib', 0],['GNA15_0.95_22_ebret_esm1b.joblib', 0],['GNAQ_0.95_31_shed_esm1b.joblib', 0],['GNA12_0.95_18_shed_esm1b.joblib',0],
            ['GNA13_0.95_18_shed_esm1b.joblib',0],['Barr1-GRK2_0.95_0_ebret_esm1b.joblib',0],
            ['Barr2_0.95_33_ebret_esm1b.joblib',0],['Barr2-GRK2_0.95_0_ebret_esm1b.joblib',0]]

    # run hmmsearch
    #os.system('hmmsearch data/7tm_1.hmm '+input_file+' > static/predictor/output/'+uniq_id+'/temp_hmm_file.txt')
    os.system('hmmsearch ' + homeDir + '/data/7tm_1.hmm '+input_file+' > ' + homeDir + '/static/predictor/output/'+uniq_id+'/temp_hmm_file.txt')
    #os.system('hmmsearch ' + homeDir + '/data/SCOP_7TM_348.hmm '+input_file+' > ' + homeDir + '/static/predictor/output/'+uniq_id+'/temp_hmm_file.txt')

    # BLAST GPCRDB
    os.system('blastp -query ' + input_file + ' -outfmt 5 -out ' + homeDir + '/static/predictor/output/' + uniq_id + '/GPCRDBblast.txt -db ' + homeDir + '/data/GPCRDB/blastdb/GPCRDB')

     # create embeddings
    for row in data:
        if "esm1b"  in  row[0].split('.')[1].split('_')[-1]:
            #model_location= "/var/www/flask_apps/precogxb_app/esm_pretrained/esm1b_t33_650M_UR50S.pt"
            #model_location= "/var/www/flask_apps/precogxb_app/esm_pretrained/esm1b_t33_650M_UR50S.pt"
            model_location = homeDir + "/esm_pretrained/esm1b_t33_650M_UR50S.pt"
            #print (model_location, 'model_location')
            if os.path.isfile(model_location) == False:
                #print ('not found')
                model_location= "esm1b_t33_650M_UR50S"


    md=[]
    for row in data:
        md=[int(row[0].split('.')[1].split('_')[1]) for row in data]
        repr_layer=list(set(md))
        repr_layer.sort()

    ##Run for all layers rather than just selected ones.
    ##This line overwrites the previous FOR loop
    repr_layer = [i for i in range(0, 34)]

    extract.main(model_location,input_file,save_path,repr_layer)

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

        Xs_test_pca_copy = predict.main(path, d, uniq_id, gprotein, input_file, input_embedding, model, int(feature_type), str(embedding))
        #np.save('static/predictor/output/'+uniq_id+'/PCA/'+gprotein, Xs_test_pca_copy)
        for name, row in zip(d, Xs_test_pca_copy):
            np.save(homeDir + '/static/predictor/output/'+uniq_id+'/PCA/'+gprotein+'_'+name, row)

    ###########################################################################
    ###########################################################################
    print ('Generating all layers for the unsupervised part')
    ## Save PCA for all layers
    ## TEST SET
    for layer in range(0, 34):
        TEST_FASTA_PATH = input_file
        TEST_EMB_PATH = input_embedding
        Xtest = []
        gpcr_test = []
        for header, _seq in esm.data.read_fasta(TEST_FASTA_PATH):
            if header.split('>')[1]:
                #print (header)
                gpcr_test.append(header.split('>')[1])
                fn = TEST_EMB_PATH+header[1:]+'.pt'
                embs = torch.load(fn)
                Xtest.append(embs['mean_representations'][layer])
                #if len(Xtest) == 50:
                 #   break
        Xtest = torch.stack(Xtest, dim=0).numpy()
        pca = load(homeDir + '/static/pca_all/pca_'+str(layer))
        Xs_test_pca = pca.transform(Xtest)
        for name, row in zip(d, Xs_test_pca):
            np.save(homeDir + '/static/predictor/output/'+uniq_id+'/PCA/'+str(layer)+'layer_'+name, row)
    print ('Done with generating all layers')
    ###########################################################################
    ###########################################################################

    for gpcr in gpcrs:
        OtherSources(gpcr, gpcrs, homeDir)

    name=str(input_file)
    #name1=name.split('/')[7].split('.')[0]
    l = '#PRECOGx\n'
    l += '#Input\tVAR'
    for gprotein in gproteins:
        l += '\t' + str(gprotein)
    l += '\n'
    for name in d:
        #print (name, 'here')
        l += '_'.join(name.split('_')[:-1]) + '\t' + name.split('_')[-1]
        for gprotein in gproteins:
            l += '\t' + str(round(d[name][gprotein],3))
        l += '\n'

        ## Put other information after WT
        if name.split('_')[-1] == 'WT':
            ## Add ebbret information
            gpcr = '_'.join(name.split('_')[:-1])
            l += gpcr + '\t' + 'IUPHAR'
            for gprotein in gproteins:
                if gprotein in gpcrs[gpcr].iuphar:
                    l += '\t' + gpcrs[gpcr].iuphar[gprotein]
                else:
                    l += '\t-'
            l += '\n'
            ## Add shedding information
            gpcr = '_'.join(name.split('_')[:-1])
            l += gpcr + '\t' + 'LogRAi'
            for gprotein in gproteins:
                if gprotein in gpcrs[gpcr].shedding:
                    l += '\t' + gpcrs[gpcr].shedding[gprotein]
                else:
                    l += '\t-'
            l += '\n'
            ## Add ebbret information
            gpcr = '_'.join(name.split('_')[:-1])
            l += gpcr + '\t' + 'Emax'
            for gprotein in gproteins:
                if gprotein in gpcrs[gpcr].ebbret:
                    l += '\t' + gpcrs[gpcr].ebbret[gprotein]
                else:
                    l += '\t-'
            l += '\n'

    #print (l)
    #shutil.rmtree('output/'+uniq_id+'/')
    open(homeDir + '/static/predictor/output/'+uniq_id+'/out.tsv', 'w').write(l)
    #print ('The output is saved at static/predictor/output/'+uniq_id)

    data = []
    dic = {}
    for line in open(homeDir + '/static/predictor/output/'+uniq_id+'/out.tsv', 'r'):
        row = []
        if line[0] != '#':
            row = line.replace('\n', '').split('\t')
            data.append(row)

    dic = {'data': data}
    with open(homeDir + '/static/predictor/output/'+uniq_id+'/out.json', 'w') as f:
        json.dump(dic, f)

    #shutil.rmtree('output/'+uniq_id+'/')
    return (uniq_id)

def OtherSources(gpcr_given, gpcrs, homeDir):
    for line in open(homeDir + '/data/shedding.tsv', 'r'):
        if line[0] != '#':
            gene_found = line.split('\t')[0]
            acc_found = line.split('\t')[1]

            pattern1 = re.compile("sp\\|.*\\|.*_.*")
            pattern2 = re.compile("tr\\|.*\\|.*_.*")

            if pattern1.match(gpcr_given) != None or pattern2.match(gpcr_given) != None:
                given_name = gpcr_given.split('|')[1]
            else:
                given_name = gpcr_given

            if gene_found == given_name or acc_found == given_name:
                values = line.replace('\n', '').split('\t')[2:]
                for value, gprotein in zip(values, gproteins):
                    gpcrs[gpcr_given].shedding[gprotein] = str(round(float(value), 2))
        else:
            gproteins = line.replace('\n', '').split('\t')[2:]

    for line in open(homeDir + '/data/ebbret.tsv', 'r', encoding="utf-8"):
        if line[0] != '#':
            gene_found = line.split('\t')[0]
            acc_found = line.split('\t')[1]

            pattern1 = re.compile("sp\\|.*\\|.*_.*")
            pattern2 = re.compile("tr\\|.*\\|.*_.*")

            if pattern1.match(gpcr_given) != None or pattern2.match(gpcr_given) != None:
                given_name = gpcr_given.split('|')[1]
            else:
                given_name = gpcr_given

            if gene_found == given_name or acc_found == given_name:
                values = line.replace('\n', '').split('\t')[2:]
                for value, gprotein in zip(values, gproteins):
                    gpcrs[gpcr_given].ebbret[gprotein] = value
        else:
            gproteins = line.replace('\n', '').split('\t')[2:]

    dic_gprot_family = {'Gs': ['GNAS', 'GNAL'],
                        'Gi/Go': ['GNAI1', 'GNAI2', 'GNAI3', 'GNAZ', 'GoA', 'GoB'],
                        'Gq/G11': ['GNAQ', 'GNA11', 'GNA14', 'GNA15'],
                        'G12/G13': ['GNA12', 'GNA13']
                        }
    for line in open(homeDir + '/data/iuphar.tsv', 'r'):
        if line[0] != '#':
            gene_found = line.split('\t')[0]
            acc_found = line.split('\t')[1]

            pattern1 = re.compile("sp\\|.*\\|.*_.*")
            pattern2 = re.compile("tr\\|.*\\|.*_.*")

            if pattern1.match(gpcr_given) != None or pattern2.match(gpcr_given) != None:
                given_name = gpcr_given.split('|')[1]
            else:
                given_name = gpcr_given

            if gene_found == given_name or acc_found == given_name:
                #print ('found in IUPHAR')
                pc_values = line.replace('\n', '').split('\t')[2]
                sc_values = line.replace('\n', '').split('\t')[3]
                for gprot_family in dic_gprot_family:
                    if gprot_family in pc_values:
                        for gprotein in dic_gprot_family[gprot_family]:
                            gpcrs[gpcr_given].iuphar[gprotein] = 'PC'
                    elif gprot_family in sc_values:
                        for gprotein in dic_gprot_family[gprot_family]:
                            gpcrs[gpcr_given].iuphar[gprotein] = 'SC'

def formatInputFileOld(numseqs, input_file):
    num = 0
    ## if input in FASTA format
    gpcrs = {}
    for line in open(input_file, 'r'):
        if num == numseqs:
            break
        if line[0] != '#':
            if line.split != []:
                if line[0] == '>':
                    print (line)
                    num += 1
                    given_name = line.replace('\n', '')

                    pattern1 = re.compile("\\>sp\\|.*\\|.*_.*")
                    pattern2 = re.compile("\\>tr\\|.*\\|.*_.*")

                    if pattern1.match(given_name) != None or pattern2.match(given_name):
                        given_name = given_name.split('>')[1].split(' ')[0]
                    else:
                        #given_name = str(line.split('>')[1].replace('\n', '').replace(' ','').split()[0])
                        given_name = str(line.split('>')[1])
                        given_name = ' '.join(given_name.split())


                    if '/' in given_name:
                        name = given_name.split('/')[0]
                        variant = given_name.split('/')[1]
                    else:
                        name = given_name
                        variant = 'WT'

                    if name not in gpcrs:
                        gpcrs[name] = GPCR(name)

                    gpcrs[name].var.append(variant)
                else:
                    print (line)
                    gpcrs[name].seq += line.replace('\n', '')

    print (gpcrs[name].seq)
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
    #print (new_input)
    #sys.exit()
    return (gpcrs, new_input)

def formatInputFile(numseqs, input_file):
    num = 0
    with open(input_file, 'r') as file:
        input = file.read()
    ## if input in FASTA format
    gpcrs = {}
    if '>' in input:
        for line in input.split('\n'):
            if num == numseqs:
                break
            if line.split() != []:
                if line[0] != '#':
                    if line[0] == '>':
                        num += 1
                        given_name = line.replace('\n', '')

                        pattern1 = re.compile("\\>sp\\|.*\\|.*_.*")
                        pattern2 = re.compile("\\>tr\\|.*\\|.*_.*")

                        if pattern1.match(given_name) != None or pattern2.match(given_name):
                            given_name = given_name.split('>')[1].split(' ')[0]

                        else:
                            #given_name = line.split('>')[1].replace('\n', '').replace(' ','').split()[0]
                            given_name = str(line.split('>')[1])
                            given_name = ' '.join(given_name.split())

                        if '/' not in given_name:
                            name = str(given_name)
                            variant = 'WT'
                        else:
                            name = '/'.join(given_name.split('/')[:-1])
                            variant = given_name.split('/')[-1]
                        if name not in gpcrs:
                            gpcrs[name] = GPCR(name)

                        gpcrs[name].var.append(str(variant))
                    else:
                        gpcrs[name].seq += line.replace('\n', '').split()[0]
    else:
        ## if input is uniprot acc or gene name
        for line in input.split('\n'):
            if num == numseqs:
                break
            if line.split() != []:
                if line[0] != '#':
                    num += 1
                    given_name = line.replace('\n', '')

                    pattern1 = re.compile("\\>sp\\|.*\\|.*_.*")
                    pattern2 = re.compile("\\>tr\\|.*\\|.*_.*")

                    if pattern1.match(given_name) != None or pattern2.match(given_name):
                        given_name = given_name.split(' ')[0]

                    else:
                        given_name = str(line.replace(' ','').split()[0])
                        #print ('this one', given_name)
                        #given_name = ' '.join(given_name.split())

                    if '/' not in given_name:
                        name = given_name
                        variant = 'WT'
                    else:
                        name = given_name.split('/')[0]
                        variant = given_name.split('/')[1]

                    if name not in gpcrs:
                        #print (name, variant, name.split())
                        gpcrs[name] = GPCR(name)
                        gpcrs[name].seq = fetchSeq(name)

                    gpcrs[name].var.append(variant)

    new_input = ''
    for gpcr in gpcrs:
        if 'WT' not in gpcrs[gpcr].var:
            gpcrs[gpcr].var.append('WT')
        for variant in gpcrs[gpcr].var:
            if gpcrs[gpcr].seq != None:
                if variant == 'WT':
                    new_input += '>'+str(gpcr)+'_WT\n'
                    new_input += str(gpcrs[gpcr].seq) + '\n'
                    #print (gpcr)
                else:
                    new_input += '>'+gpcr+'_'+variant+'\n'
                    position = int(variant[1:-1]) - 1
                    newAA = variant[-1]
                    #print (position, newAA)
                    new_input += gpcrs[gpcr].seq[:position] + newAA + gpcrs[gpcr].seq[position+1:] + '\n'

    #print (new_input)
    return (gpcrs, new_input)

def formatInput(numseqs, input):
    num = 0
    ## if input in FASTA format
    gpcrs = {}
    if '>' in input:
        for line in input.split('\n'):
            if num == numseqs:
                break
            if line.split() != []:
                if line[0] != '#':
                    if line[0] == '>':
                        num += 1
                        given_name = line.replace('\n', '')

                        pattern1 = re.compile("\\>sp\\|.*\\|.*_.*")
                        pattern2 = re.compile("\\>tr\\|.*\\|.*_.*")

                        if pattern1.match(given_name) != None or pattern2.match(given_name):
                            given_name = given_name.split('>')[1].split(' ')[0]

                        else:
                            #given_name = line.split('>')[1].replace('\n', '').replace(' ','').split()[0]
                            given_name = str(line.split('>')[1])
                            given_name = ' '.join(given_name.split())

                        if '/' not in given_name:
                            name = str(given_name)
                            variant = 'WT'
                        else:
                            name = '/'.join(given_name.split('/')[:-1])
                            variant = given_name.split('/')[-1]
                        if name not in gpcrs:
                            gpcrs[name] = GPCR(name)

                        gpcrs[name].var.append(str(variant))
                    else:
                        gpcrs[name].seq += line.replace('\n', '').split()[0]
    else:
        ## if input is uniprot acc or gene name
        for line in input.split('\n'):
            if num == numseqs:
                break
            if line.split() != []:
                if line[0] != '#':
                    num += 1
                    given_name = line.replace('\n', '')

                    pattern1 = re.compile("\\>sp\\|.*\\|.*_.*")
                    pattern2 = re.compile("\\>tr\\|.*\\|.*_.*")

                    if pattern1.match(given_name) != None or pattern2.match(given_name):
                        given_name = given_name.split(' ')[0]

                    else:
                        given_name = str(line.replace(' ','').split()[0])
                        #print ('this one', given_name)
                        #given_name = ' '.join(given_name.split())

                    if '/' not in given_name:
                        name = given_name
                        variant = 'WT'
                    else:
                        name = given_name.split('/')[0]
                        variant = given_name.split('/')[1]

                    if name not in gpcrs:
                        #print (name, variant, name.split())
                        gpcrs[name] = GPCR(name)
                        gpcrs[name].seq = fetchSeq(name)

                    gpcrs[name].var.append(variant)

    #print (name, gpcrs[name].seq)
    new_input = ''
    for gpcr in gpcrs:
        if 'WT' not in gpcrs[gpcr].var:
            gpcrs[gpcr].var.append('WT')
        for variant in gpcrs[gpcr].var:
            if gpcrs[gpcr].seq != None:
                if variant == 'WT':
                    new_input += '>'+str(gpcr)+'_WT\n'
                    new_input += str(gpcrs[gpcr].seq) + '\n'
                    #print (gpcr)
                else:
                    new_input += '>'+gpcr+'_'+variant+'\n'
                    position = int(variant[1:-1]) - 1
                    newAA = variant[-1]
                    #print (position, newAA)
                    new_input += gpcrs[gpcr].seq[:position] + newAA + gpcrs[gpcr].seq[position+1:] + '\n'

    #print (new_input)
    return (gpcrs, new_input)

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
    parser.add_argument('assay', help='Input what assay is used (all or ebret or shed)')
    parser.add_argument('--file', help='Input File (FASTA/Mechismo format). Applicable for both webApp and command-line versions')
    parser.add_argument('--input', help='Input (Mechismo format). Applicable only when accessed via the webApp')
    parser.add_argument('--numseqs', help='Num of seqs allowed (default: 15)')
    args = parser.parse_args()
    assay = args.assay
    input = args.input
    input_file = args.file
    numseqs = args.numseqs
    if numseqs == None:
        numseqs = 25
    else:
        numseqs = int(numseqs)
    uniq_id = main(numseqs, input, input_file, assay, os.getcwd())
    #uniq_id = main(numseqs, input, input_file, assay, '')
    print ('Done')
