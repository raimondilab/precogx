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

TEST_FASTA_PATH = 'static/predictor/output/86T3G/input.fasta'
TEST_EMB_PATH = 'static/predictor/output/86T3G/embed/'
TEST_ATTN_PATH = 'static/predictor/output/86T3G/attentions/'
ATTN_HEAD = 5

Xtest = []
XtestA = []
gpcr_test = []
num = 0
for header, _seq in esm.data.read_fasta(TEST_FASTA_PATH):
    if header.split('>')[1]:
        gpcr_test.append(header.split('>')[1])
        '''
        fn = TEST_EMB_PATH+header[1:]+'.pt'
        embs = torch.load(fn)
        Xtest.append(embs['mean_representations'][EMB_LAYER])
        '''
        ## attentions
        attns = TEST_ATTN_PATH + header[1:]+'.pt'
        embsA = torch.load(attns)
        #print (embsA.size(dim=0))
        #print (len(embsA.size()))
        if len(embsA.size()) == 5:
            XtestA.append(embsA[num][32][ATTN_HEAD])
            print (len(embsA[num][32][ATTN_HEAD].size()))
            #print (embsA[num][32][ATTN_HEAD])
        else:
            XtestA.append(embsA[32][ATTN_HEAD])
            print (len(embsA[32][ATTN_HEAD].size()))
        num += 1
        #break
        #if len(Xtest) == 50:
         #   break
XtestA = torch.stack(XtestA, dim=0).numpy()
print (XtestA[0][0])
print (XtestA[0][1])
print (len(XtestA))
