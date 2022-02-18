import os, sys
import numpy as np

#path = 'embedding_all/'
path = 'embedding_without_taste/'
for files in os.listdir(path):
    if files.endswith('.npy'):
        X = np.load(path+files)
        #print (X.tolist())
        #print (len(X.tolist()))
        #print (len(X.tolist()[0]))
        l = ''
        for row in X:
            for num, value in enumerate(row):
                if num == len(row) - 1:
                    l += str(value)
                else:
                    l += str(value) + '\t'
            l += '\n'
        name = files.split('.')[0]
        open(path+name+'_embedding.tsv', 'w').write(l)
        #break
#print (l)
