from flask import Flask, render_template, request, jsonify, redirect, url_for
import os, sys, json
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from Bio import SearchIO
from Bio.Blast import NCBIXML
from matplotlib import cm
import pandas as pd
#sys.path.insert(1, 'static/predictor/')
#from precogxb_app.static.predictor import precogx

app = Flask(__name__)

##
path = os.getcwd()
path = app.root_path
sys.path.insert(1, path + '/static/predictor/')
import precogx

## Route to home page
@app.route('/', methods=['GET', 'POST'])
@app.route('/home', methods=['GET', 'POST'])
def home():
    return render_template('index.html')

@app.route('/error', methods=['GET', 'POST'])
def error():
    return render_template('error.html')

def sortPositions(positions):
    data = []
    for position in positions:
        if position not in ['Nterm', 'Cterm']:
            '''
            if '.' in position:
                base = position.split('.')[0]
                #value = int(position.split('.')[1])
                value = int(position.split('.')[1])*(10**(len(position.split('.')[1])*(-1)))
            else:
                base = position
                value = 0
            '''
            if 'x' in position:
                base = position.split('x')[0]
                #value = int(position.split('.')[1])
                value = int(position.split('x')[1])*(10**(len(position.split('x')[1])*(-1)))
            else:
                base = position
                value = 0

            if 'ICL' in base:
                if base == 'ICL1':
                    order = 1.5
                elif base == 'ICL2':
                    order = 3.5
                else:
                    order = 5.5
            elif 'ECL' in base:
                if base == 'ECL1':
                    order = 2.5
                elif base == 'ECL2':
                    order = 4.5
                else:
                    order = 6.5
            else:
                order = int(base)

        else:
            if position == 'Nterm':
                base = 0
                order = 0
                value = 0
            else:
                base = 9
                order = 9
                value = 0
        #print (position, base, order, value)
        row = []
        row.append(position)
        row.append(order)
        row.append(base)
        row.append(value)
        data.append(row)

    df = pd.DataFrame(data, columns = ['position', 'order', 'base', 'value'])
    df = df.sort_values(['order','value'],ascending=[True, True])
    #print (df.to_numpy())
    positions = []
    for row in df.to_numpy():
        positions.append(row[0])

    return (np.array(positions))

def extract_contacts(gprotein_given, cutoff, distance):
    assay = ''
    for line in open(path+'/data/contacts/gprotein_best_layer.txt', 'r'):
        if gprotein_given == line.split('\t')[0]:
            assay = line.split('\t')[1]
    #print (gprotein_given)
    dic = {}; positions = []; pair_positions = []; scores = [];
    '''
    for line in open(path+'/data/contacts/all_positions_count_'+assay+'_scaled_web_old.txt', 'r'):
        gprotein_found = line.split('\t')[0]
        if gprotein_given == gprotein_found:
            #print ('here')
            pos1 = line.split('\t')[-2]
            pos2 = line.replace('\n', '').split('\t')[-1]
            score = float(line.replace('\n','').split('\t')[1])
            if score >= cutoff or score <= (-1.0)*cutoff:
                if pos1 not in dic:
                    dic[pos1] = {}
                dic[pos1][pos2] = score

                if pos2 not in dic:
                    dic[pos2] = {}
                dic[pos2][pos1] = score

                if pos1 not in positions:
                    positions.append(pos1)
                #positions.append(pos1)
                if pos2 not in positions:
                    positions.append(pos2)
                pair_positions.append(pos1+':'+pos2+':'+str(score))
            scores.append(score)
    '''
    #for line in open(path+'/data/contacts/position_'+assay+'_scaled_web_new.txt', 'r'):
    for line in open(path+'/data/contacts/position_'+assay+'_scaled_web_new2.txt', 'r'):
        gprotein_found = line.split('\t')[0]
        if gprotein_given == gprotein_found:
            #print ('here')
            #print (line)
            pos1 = line.split('\t')[1]
            pos2 = line.split('\t')[2]
            score = float(line.split('\t')[3])
            dis = line.replace('\n', '').split('\t')[-1]
            if dis == '-':
                dis = 1000.0 ## set by default a high value so that it is anyway selected
            else:
                dis = float(dis)
            if score >= cutoff or score <= (-1.0)*cutoff and dis >= distance:
                if pos1 not in dic:
                    dic[pos1] = {}
                dic[pos1][pos2] = score

                if pos2 not in dic:
                    dic[pos2] = {}
                dic[pos2][pos1] = score

                if pos1 not in positions:
                    positions.append(pos1)
                #positions.append(pos1)
                if pos2 not in positions:
                    positions.append(pos2)
                pair_positions.append(pos1+':'+pos2+':'+str(score))
            scores.append(score)

    #print ('positions', len(positions))
    #print ('distance', distance, dis)
    scoresMax = max(scores)
    scoresMin = min(scores)
    positions = np.array(positions)
    positions = sortPositions(positions)
    #print ('----------------')
    #print (positions)
    #print ('----------------')
    #sys.exit()
    #positions = list(set(positions))
    #positions = np.array(np.sort(positions))
    data = []
    num_contacts = []
    for pos1 in positions:
        row = []
        for pos2 in positions:
            if pos2 in dic[pos1]:
                row.append(dic[pos1][pos2])
            else:
                #row.append(0)
                row.append(None)
        data.append(row)
        num_contacts.append([round(len(dic[pos1]),2)])
    #print (cutoff, gprotein_given)
    #print ('positions', positions)

    if num_contacts != []:
        scaler = MinMaxScaler(feature_range=(0.35,1.0))
        num_contacts = scaler.fit_transform(num_contacts)
        num_contacts = num_contacts.flatten().tolist()
        for i in range(0, len(num_contacts)):
            num_contacts[i] = str(num_contacts[i])

    #print (num_contacts)

    return scoresMax, scoresMin, data, positions, pair_positions, num_contacts

@app.route('/fetchAttentionMap', methods=['GET', 'POST'])
def fetchAttentionMap():
    if request.method == 'POST':
        data = request.get_json(force=True)
        #print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        #scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff, distance)
        Xtest = np.load(path+'/static/predictor/output/'+uniq_id+'/attentions/'+gpcr_given+'_'+gprotein_given+'.npy')
        seqPositions = [str(i) for i in range(1, len(Xtest[0])+1)]
        #print (seqPositions)
        #return jsonify({'fetch_contactsMin': scoresMin, 'fetch_contactsMax': scoresMax, 'fetch_contacts': scores, 'positions': positions.tolist()})
        return jsonify({'zaxis': Xtest.tolist(),
                        'xaxis': seqPositions,
                        'yaxis': seqPositions,
                        'gpcr_name': '_'.join(gpcr_given.split('_')[:-1]) + '/' + gpcr_given.split('_')[-1]
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")

#
@app.route('/fetchContactsHeatmap', methods=['GET', 'POST'])
def fetchContactsHeatmap():
    if request.method == 'POST':
        data = request.get_json(force=True)
        #print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        cutoff = float(data['cutoff'])
        distance = float(data['distance'])
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff, distance)
        for i in range(0, len(positions)):
            positions[i] = str(positions[i]).replace('.', 'x')
        return jsonify({'fetch_contactsMin': scoresMin,
                        'fetch_contactsMax': scoresMax,
                        'fetch_contacts': scores,
                        'positions': positions.tolist(),
                        'gpcr_name': '_'.join(gpcr_given.split('_')[:-1]) + '/' + gpcr_given.split('_')[-1]
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")

def extract_pca(gprotein, assay, pca_type):
    '''
    if pca_type == 'Best PCA':
        Xs_train_pca = np.load(path+'/static/best_PCA/'+gprotein+'.npy', allow_pickle=True)
    elif pca_type == 'GPCRome':
        Xs_train_pca = np.load(path+'/static/33layer_PCA/33layer.npy', allow_pickle=True)
    else:
    '''
    Xs_train_pca = np.load(path+'/static/pca_all/'+pca_type+'.npy', allow_pickle=True)
    #Xs_train_pca = np.load(path+'/static/best_PCA/GNAZ.npy', allow_pickle=True)
    score_coupling, score_uncoupling, Xs_train_pca_coupling, Xs_train_pca_uncoupling, Xs_train_pca_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = filter_gpcr_list(Xs_train_pca, assay, gprotein)
    #print ('train', Xs_train_pca_coupling)
    score_coupling = score_coupling.tolist()
    score_uncoupling = score_uncoupling.tolist()
    x_train_coupling = Xs_train_pca_coupling[:,0].tolist()
    x_train_uncoupling = Xs_train_pca_uncoupling[:,0].tolist()
    x_train_grey = Xs_train_pca_grey[:,0].tolist()
    y_train_coupling = Xs_train_pca_coupling[:,1].tolist()
    y_train_uncoupling = Xs_train_pca_uncoupling[:,1].tolist()
    y_train_grey = Xs_train_pca_grey[:,1].tolist()
    return score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey

def filter_gpcr_list(X, assay, gprotein):
    genes_to_consider_coupling = []
    score_coupling = []
    score_uncoupling = []
    genes_to_consider_uncoupling = []
    #print (assay)
    #assay = 'ebBRET'
    if assay == 'TGF':
        num = -1
        for line in open(path+'/data/shedding.tsv', 'r'):
            if line[0] != '#':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                id = line.split('\t')[2]
                name = line.split('\t')[3]
                # Make the range -1 to +1
                score = float(line.split('\t')[num+4]) + 1
                #print (score)
                #color = cm.get_cmap('RdYlGn', 100)
                #r,g,b,a = color(score)
                #print (score,r,g,b)
                #if score >= -1.0:
                # Move the center (-1.0) to 0.0
                if score >= 0.0:
                    genes_to_consider_coupling.append(gene+'|'+acc)
                    color = cm.get_cmap('Greens', 100)
                    r,g,b,a = color(score)
                    score_coupling.append('rgb('+str(r)+','+str(g)+','+str(b)+')')
                else:
                    genes_to_consider_uncoupling.append(gene+'|'+acc)
                    color = cm.get_cmap('Greys', 100)
                    score *= (-1.0)
                    # To bring the score from range 0 to 1
                    # to 0.25 to 0.75 so it is neither too white (low value)
                    # not too black (high value)
                    score = score/2 + 0.25
                    r,g,b,a = color(score)
                    score_uncoupling.append('rgb('+str(r)+','+str(g)+','+str(b)+')')
            else:
                flag = 0
                # Replace GNAO1 by GoA
                header = line.replace('\n', '').replace('GNAO1', 'GoA').split('\t')[4:]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        flag = 1
                        break
                #print (num, gprot)
                #print (header)
                if flag == 0:
                    num = -1

    elif assay == 'GEMTA':
        num = 0
        for line in open(path+'/data/ebbret.tsv', 'r', encoding="utf-8"):
            if line[0] != '#':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                id = line.split('\t')[2]
                name = line.split('\t')[3]
                score = float(line.split('\t')[num+4])

                if score > 0.0:
                    genes_to_consider_coupling.append(gene+'|'+acc)
                    color = cm.get_cmap('Greens', 100)
                    r,g,b,a = color(score)
                    score_coupling.append('rgb('+str(r)+','+str(g)+','+str(b)+')')
                else:
                    genes_to_consider_uncoupling.append(gene+'|'+acc)
                    #score_uncoupling.append('rgb('+str(r)+','+str(g)+','+str(b)+')')
                    score_uncoupling.append('grey')
            else:
                header = line.replace('\n', '').split('\t')[4:]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        #print (gprot)
                        break

    elif assay == 'GtoPdb':
        iuphar_map = {
                      'GNAS': 'Gs', 'GNAL': 'Gs',
                      'GNAI1': 'Gi/Go', 'GNAI2': 'Gi/Go', 'GNAI3': 'Gi/Go', 'GNAO1': 'Gi/Go', 'GNAZ': 'Gi/Go', 'GoA': 'Gi/Go', 'GoB': 'Gi/Go',
                      'GNA12': 'G12/G13', 'GNA13': 'G12/G13',
                      'GNAQ': 'Gq/G11', 'GNA11': 'Gq/G11', 'GNA14': 'Gq/G11', 'GNA15': 'Gq/G11'
                      }
        gprotein_fam = iuphar_map[gprotein]
        #print (gprotein_fam)
        for line in open(path+'/data/iuphar.tsv', 'r'):
            if line[0] != '#' and line.split('\t')[1] != '':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                id = line.split('\t')[2]
                name = line.split('\t')[3]
                if gprotein_fam in line:
                    genes_to_consider_coupling.append(gene+'|'+acc)
                    if gprotein_fam in line.split('\t')[4]:
                        score_coupling.append('forestgreen')
                    else:
                        score_coupling.append('lightgreen')
                else:
                    genes_to_consider_uncoupling.append(gene+'|'+acc)
                    score_uncoupling.append('grey')

    elif assay == 'STRING':
        string_map = {
                      'Barr1-GRK2': 'ARRB1',
                      'Barr2': 'ARRB2',
                      'Barr2-GRK2': 'ARRB2'
                      }
        barr = string_map[gprotein]
        num = -1
        for line in open(path+'/data/string.tsv', 'r', encoding="utf-8"):
            gene = line.split('\t')[0]
            acc = line.split('\t')[1]
            id = line.split('\t')[2]
            name = line.split('\t')[3]
            if barr in line:
                genes_to_consider_coupling.append(gene+'|'+acc)
                score_coupling.append('green')
            else:
                genes_to_consider_uncoupling.append(gene+'|'+acc)
                score_uncoupling.append('green')
    #print (genes_to_consider_coupling)
    #print (genes_to_consider_uncoupling)

    gpcr_list = []
    for line in open(path+'/static/pca_all/gpcr_list_unscaled_all_layer_GN.txt', 'r'):
        gene = line.replace('\n', '').split('\t')[1]
        acc = line.replace('\n', '').split('\t')[0]
        gpcr_list.append(gene + '|' + acc)
        #gpcr_list.append(line.replace('\n', '').split('\t')[1])

    X_pos = []
    X_neg = []
    X_grey = []
    genes_to_consider_grey = []
    for gene, row in zip(gpcr_list, X):
        #if 'MC1R' in gene:
        #    print ('couple', gene, row[:2])
        if gene in genes_to_consider_coupling:
            X_pos.append(row)
        elif gene in genes_to_consider_uncoupling:
            X_neg.append(row)
        else:
            X_grey.append(row)
            genes_to_consider_grey.append(gene)

    #print (X_pos)

    return (np.array(score_coupling), np.array(score_uncoupling), np.array(X_pos), np.array(X_neg), np.array(X_grey), genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey)

@app.route('/fetchPCA', methods=['GET', 'POST'])
def fetchPCA():
    if request.method == 'POST':
        data = request.get_json(force=True)
        assay_given = data['assay']
        pca_type = data['pca_type']
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        #print (gprotein_given, gpcr_given)
        uniq_id = data['uniq_id']

        #if assay == '':
        assay = '';
        assayList = []
        gemta = ['GNAS', 'GNAI1', 'GNAI2', 'GoB', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA11', 'GNA14', 'GNA15', 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2']
        tgf = ['GNAS', 'GNAL', 'GNAI1', 'GNAI3', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15']
        both = ['GNAS', 'GNAI1', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15', 'GoA']

        if 'Barr' in gprotein_given:
            assay = 'GEMTA'
            assayList = ['GEMTA', 'STRING', 'Class']
        elif gprotein_given in both:
            assay = 'TGF'
            assayList = ['TGF', 'GEMTA', 'GtoPdb', 'Class']
        elif gprotein_given in tgf:
            assay = 'TGF'
            assayList = ['TGF', 'GtoPdb', 'Class']
        elif gprotein_given in gemta:
            assay = 'GEMTA'
            assayList = ['GEMTA', 'GtoPdb', 'Class']
        elif assay_given == 'Class':
            assay = 'Class'
            assayList = ['TGF', 'GEMTA', 'GtoPdb', 'Class']

        if assay_given in assayList:
            assay = assay_given

        ### MUT
        '''
        if pca_type == 'GPCRome':
            Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/33layer_'+gpcr_given+'.npy', allow_pickle=True)
        elif pca_type == 'Best PCA':
            Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+gpcr_given+'.npy', allow_pickle=True)
        else:
        '''
        '''
        Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+pca_type+'layer_'+gpcr_given+'.npy', allow_pickle=True)
        #print ('test',Xs_test_pca)
        x_test = Xs_test_pca[0].tolist()
        y_test = Xs_test_pca[1].tolist()
        #print (x_test)
        #print (y_test)
        '''

        ## fetch all variants of the given GPCR
        gpcrName = '_'.join(gpcr_given.split('_')[:-1])
        #print (gpcrName)
        variantsToConsider = []; x_test = []; y_test = []; test_names = [];
        if gpcr_given[:-3] != '_WT':
            for files in os.listdir(path+'/static/predictor/output/'+uniq_id+'/PCA/'):
                if files.endswith('.npy') and 'layer' in files  and files.split('_')[-1].split('.')[0] != 'WT':
                    if int(pca_type) == int(files.split('layer')[0]):
                        if gpcrName == '_'.join(files.split('_')[1:-1]):
                            variantsToConsider.append(files)
            
            for files in variantsToConsider:
                Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+files, allow_pickle=True)
                #print ('test',Xs_test_pca)
                x_test.append(Xs_test_pca[0].tolist())
                y_test.append(Xs_test_pca[1].tolist())
                #test_names.append(files.split('_')[1]+'_'+files.split('_')[-1].split('.')[0])
                test_names.append('_'.join(files.split('_')[1:-1]) + '/' + files.split('_')[-1].split('.')[0])


        ### WT
        #wt = gpcr_given.split('_')[0] + '_WT'
        #print ('_'.join(gpcr_given.split('_')[:-1]))
        wt = '_'.join(gpcr_given.split('_')[:-1]) + '_WT'
        wt_name = '_'.join(gpcr_given.split('_')[:-1]) + '/WT'
        #print (wt_name)
        #print (test_names)
        '''
        if pca_type == 'GPCRome':
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/33layer_'+wt+'.npy', allow_pickle=True)
        elif pca_type == 'Best PCA':
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+wt+'.npy', allow_pickle=True)
        else:
        '''
        Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+pca_type+'layer_'+wt+'.npy', allow_pickle=True)
        #print (Xs_wt_pca)
        #x_test = Xs_test_pca[:,0].tolist()
        #y_test = Xs_test_pca[:,1].tolist()
        x_wt = Xs_wt_pca[0].tolist()
        y_wt = Xs_wt_pca[1].tolist()
        #print (x_wt)
        #print (y_wt)
        
        score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = extract_pca(gprotein_given, assay, pca_type)
        #print (x_train, y_train, x_test, y_test)
        #print (assay,genes_to_consider_coupling)
        minX = min(x_train_coupling + x_train_uncoupling + x_train_grey)
        maxX = max(x_train_coupling + x_train_uncoupling + x_train_grey)
        minY = min(y_train_coupling + y_train_uncoupling + x_train_grey)
        maxY = max(y_train_coupling + y_train_uncoupling + x_train_grey)
        #print(minY, maxY)
        return jsonify({'x_train_coupling': x_train_coupling,
                        'x_train_uncoupling': x_train_uncoupling,
                        'y_train_coupling': y_train_coupling,
                        'y_train_uncoupling': y_train_uncoupling,
                        'x_train_grey': x_train_grey,
                        'y_train_grey': y_train_grey,
                        'score_coupling': score_coupling,
                        'score_uncoupling': score_uncoupling,
                        'assay': assay,
                        'genes_to_consider_coupling': genes_to_consider_coupling,
                        'genes_to_consider_uncoupling': genes_to_consider_uncoupling,
                        'genes_to_consider_grey': genes_to_consider_grey,
                        'x_test': x_test,
                        'y_test': y_test,
                        'test_names': test_names,
                        'wt_name': wt_name,
                        'x_wt': x_wt,
                        'y_wt': y_wt,
                        'minX': minX,
                        'maxX': maxX,
                        'minY': minY,
                        'maxY': maxY})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

@app.route('/fetchPCA2', methods=['GET', 'POST'])
def fetchPCA2():
    if request.method == 'POST':
        data = request.get_json(force=True)
        assay_given = data['assay']
        pca_type = data['pca_type']
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        #print (gprotein_given, gpcr_given)
        uniq_id = data['uniq_id']

        #if assay == '':
        assay = assay_given;

        ### MUT
        '''
        if pca_type == 'GPCRome':
            Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/33layer_'+gpcr_given+'.npy', allow_pickle=True)
        elif pca_type == 'Best PCA':
            Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+gpcr_given+'.npy', allow_pickle=True)
        else:
        '''
        '''
        Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+pca_type+'layer_'+gpcr_given+'.npy', allow_pickle=True)
        #print ('test',Xs_test_pca)
        x_test = Xs_test_pca[0].tolist()
        y_test = Xs_test_pca[1].tolist()
        #print (x_test)
        #print (y_test)
        '''
        ## fetch all variants of the given GPCR
        gpcrName = '_'.join(gpcr_given.split('_')[:-1])
        #print (gpcrName)
        variantsToConsider = []; x_test = []; y_test = []; test_names = [];
        if gpcr_given[:-3] != '_WT':
            for files in os.listdir(path+'/static/predictor/output/'+uniq_id+'/PCA/'):
                if files.endswith('.npy') and 'layer' in files  and files.split('_')[-1].split('.')[0] != 'WT':
                    if int(pca_type) == int(files.split('layer')[0]):
                        if gpcrName == '_'.join(files.split('_')[1:-1]):
                            variantsToConsider.append(files)
            
            for files in variantsToConsider:
                Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+files, allow_pickle=True)
                #print ('test',Xs_test_pca)
                x_test.append(Xs_test_pca[0].tolist())
                y_test.append(Xs_test_pca[1].tolist())
                #test_names.append(files.split('_')[1]+'_'+files.split('_')[-1].split('.')[0])
                test_names.append('_'.join(files.split('_')[1:-1]) + '/' + files.split('_')[-1].split('.')[0])


        ### WT
        #wt = gpcr_given.split('_')[0] + '_WT'
        #print ('_'.join(gpcr_given.split('_')[:-1]))
        wt = '_'.join(gpcr_given.split('_')[:-1]) + '_WT'
        wt_name = '_'.join(gpcr_given.split('_')[:-1]) + '/WT'
        '''
        if pca_type == 'GPCRome':
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/33layer_'+wt+'.npy', allow_pickle=True)
        elif pca_type == 'Best PCA':
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+wt+'.npy', allow_pickle=True)
        else:
        '''
        Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+pca_type+'layer_'+wt+'.npy', allow_pickle=True)
        #print (Xs_wt_pca)
        #x_test = Xs_test_pca[:,0].tolist()
        #y_test = Xs_test_pca[:,1].tolist()
        x_wt = Xs_wt_pca[0].tolist()
        y_wt = Xs_wt_pca[1].tolist()
        #print (x_wt)
        #print (y_wt)
        
        #score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = extract_pca(gprotein_given, assay, pca_type)
        '''
        if pca_type == 'Best PCA':
            Xs_train_pca = np.load(path+'/static/best_PCA/'+gprotein_given+'.npy', allow_pickle=True)
        elif pca_type == 'GPCRome':
            Xs_train_pca = np.load(path+'/static/33layer_PCA/33layer.npy', allow_pickle=True)
        else:
            '''
        Xs_train_pca = np.load(path+'/static/pca_all/'+pca_type+'.npy', allow_pickle=True)

        classes = {}
        #for line in open(path+'/data/classification.txt', 'r'):
        for line in open(path+'/data/classification2.txt', 'r'):
            if 'Uniprot_acc' not in line:
                acc = line.split('\t')[0]
                cls = line.split('\t')[-1].replace('\n', '')
                classes[acc] = cls

        gpcr_list = []
        for line in open(path+'/static/pca_all/gpcr_list_unscaled_all_layer_GN.txt', 'r'):
            gene = line.replace('\n', '').split('\t')[1]
            acc = line.replace('\n', '').split('\t')[0]
            gpcr_list.append(gene + '|' + acc)
            #gpcr_list.append(line.replace('\n', '').split('\t')[1])

        X_classA = []
        classA = []
        X_classB1 = []
        classB1 = []
        X_classB2 = []
        classB2 = []
        X_classC = []
        classC = []
        X_frizzeled = []
        frizzeled = []
        X_taste = []
        taste = []
        X_other = []
        other = []
        for gpcr, row in zip(gpcr_list, Xs_train_pca):
            if 'MC1R' in gpcr:
                print ('class', gpcr, row[:2])
            acc = gpcr.split('|')[1]
            if classes[acc] == 'classA':
                X_classA.append(row.tolist())
                classA.append(gpcr)
            elif classes[acc] == 'classB1':
                X_classB1.append(row)
                classB1.append(gpcr)
            elif classes[acc] == 'classB2':
                X_classB2.append(row)
                classB2.append(gpcr)
            elif classes[acc] == 'classC':
                X_classC.append(row)
                classC.append(gpcr)
            elif classes[acc] == 'Frizzeled':
                X_frizzeled.append(row)
                frizzeled.append(gpcr)
            elif classes[acc] == 'Taste':
                X_taste.append(row)
                taste.append(gpcr)
            else:
                X_other.append(row.tolist())
                other.append(gpcr)

        X_classA = np.array(X_classA)
        X_classB1 = np.array(X_classB1)
        X_classB2 = np.array(X_classB2)
        X_classC = np.array(X_classC)
        X_frizzeled = np.array(X_frizzeled)
        X_taste = np.array(X_taste)
        X_other = np.array(X_other)
        #print (X_classA)

        x_classA = X_classA[:,0].tolist()
        y_classA = X_classA[:,1].tolist()
        x_classB1 = X_classB1[:,0].tolist()
        y_classB1 = X_classB1[:,1].tolist()
        x_classB2 = X_classB2[:,0].tolist()
        y_classB2 = X_classB2[:,1].tolist()
        x_classC = X_classC[:,0].tolist()
        y_classC = X_classC[:,1].tolist()
        x_frizzeled = X_frizzeled[:,0].tolist()
        y_frizzeled = X_frizzeled[:,1].tolist()
        x_taste = X_taste[:,0].tolist()
        y_taste = X_taste[:,1].tolist()
        x_other = X_other[:,0].tolist()
        y_other = X_other[:,1].tolist()

        '''
        for gpcr, x, y in zip(classA, x_classA, y_classA):
            if 'MC1R' in gpcr:
                print ('classA', gpcr, x, y)
                print (len(x_classA), len(y_classA))
        '''

        #print (y_classA)
        #print (genes_to_consider_coupling)
        minX = min(x_classA + x_classB1 + x_classB2 + x_classC + x_frizzeled + x_taste + x_other)
        maxX = max(x_classA + x_classB1 + x_classB2 + x_classC + x_frizzeled + x_taste + x_other)
        minY = min(y_classA + y_classB1 + y_classB2 + y_classC + y_frizzeled + y_taste + y_other)
        maxY = max(y_classA + y_classB1 + y_classB2 + y_classC + y_frizzeled + y_taste + y_other)
        #print(minX, maxX)

        #print (x_test, y_test)
        return jsonify({'x_classA': x_classA,
                        'x_classB1': x_classB1,
                        'x_classB1': x_classB2,
                        'x_classC': x_classC,
                        'x_frizzeled': x_frizzeled,
                        'x_taste': x_taste,
                        'x_other': x_other,
                        'y_classA': y_classA,
                        'y_classB1': y_classB1,
                        'y_classB2': y_classB2,
                        'y_classC': y_classC,
                        'y_frizzeled': y_frizzeled,
                        'y_taste': y_taste,
                        'y_other': y_other,
                        'assay': assay,
                        'classA': classA,
                        'classB1': classB1,
                        'classB2': classB2,
                        'classC': classC,
                        'frizzeled': frizzeled,
                        'taste': taste,
                        'other': other,
                        'x_test': x_test,
                        'y_test': y_test,
                        'test_names': test_names,
                        'wt_name': wt_name,
                        'x_wt': x_wt,
                        'y_wt': y_wt,
                        'minX': minX,
                        'maxX': maxX,
                        'minY': minY,
                        'maxY': maxY})
    else:
        return ("<html><h3>It was a GET request</h3></html>")


@app.route('/assignOptions', methods=['GET', 'POST'])
def assignOptions():
    if request.method == 'POST':
        data = request.get_json(force=True)
        assay_given = data['assay']
        pca_type = data['pca_type']
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        #print (gprotein_given, gpcr_given)
        uniq_id = data['uniq_id']

        for line in open(path+'/data/contacts/gprotein_best_layer.txt', 'r'):
            if gprotein_given == line.split('\t')[0]:
                assay = line.split('\t')[1]
                bestPCA = line.split('\t')[2].replace('\n', '')
                break

        layers = [str(i) for i in range(0,34)]
        '''
        layers = []
        for i in range(0,34):
            if str(i) == bestPCA:
                layers.append(str(i)+' ('+'best)')
            else:
                layers.append(str(i))
        '''

        return jsonify({'layers': layers})


def DoBLAST(uniq_id, gpcr_given):
    handle = open(path + "/static/predictor/output/"+uniq_id+"/GPCRDBblast.txt", 'r')
    blast_records = NCBIXML.parse(handle)
    #print (blast_records)

    GPCRDB2SEQ = {}
    SEQ2GPCRDB = {}
    for blast_record in blast_records:
        #print (blast_record.query)
        if gpcr_given == blast_record.query:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    bestHIT = alignment.title.split(' ')[1]
                    q_num = 0
                    s_num = 0
                    for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                        if q!='-' and s!='-':
                            GPCRDB2SEQ[s_num + hsp.sbjct_start] = q_num + hsp.query_start
                            SEQ2GPCRDB[q_num + hsp.query_start] = s_num + hsp.sbjct_start
                            q_num += 1
                            s_num += 1
                        elif q!='-':
                            q_num += 1
                        else:
                            s_num += 1
                    '''
                    for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                        if q!='-' and s!='-':
                            GPCRDB2SEQ[num + hsp.sbjct_start] = num + hsp.query_start
                    '''
                    break
                break
            #print (bestHIT)
            #print (GPCRDB2SEQ)
            break
    return (GPCRDB2SEQ, SEQ2GPCRDB, bestHIT)

@app.route('/fetchContactsSequence', methods=['GET', 'POST'])
def fetchContactsSequence():
    if request.method == 'POST':
        data = request.get_json(force=True)
        #print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        #print (gpcr_given)
        #print (gpcr_given, gpcr_given.split('_')[1], 'sequence')
        path_to_fasta = data['path_to_fasta']
        uniq_id = data['uniq_id']
        cutoff = float(data['cutoff'])
        distance = float(data['distance'])
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff, distance)

        #print ('fetch_seq', positions)
        fasta_sequence = ''; flag = 0
        for line in open(path_to_fasta):
            if line[0] == '>':
                if flag == 1:
                    break
                flag = 0
                #gpcr_found = line.split('>')[1].replace('\n', '').replace(' ', '')
                gpcr_found = line.split('>')[1].replace('\n', '').lstrip().rstrip()
                if gpcr_found == gpcr_given:
                    flag = 1
            elif flag == 1:
                fasta_sequence += line.replace('\n', '')

        GPCRDB2SEQ, SEQ2GPCRDB, bestHIT = DoBLAST(uniq_id, gpcr_given)
        bestHIT_ACC = bestHIT.split('|')[1]
        #print (bestHIT)

        BW2GPCRDB = {}
        GPCRDB2BW = {}
        for line in open(path + '/data/GPCRDB/GPCRDB.tsv', 'r'):
            if 'Name' not in line.split('\t')[0]:
                acc = line.split('\t')[0].split('_')[1]
                if acc == bestHIT_ACC:
                    GPCRDB = int(line.split('\t')[1][1:])
                    #BW = line.split('\t')[2]
                    ## convert . to x in GPCRDB numbering
                    BW = line.split('\t')[2].replace('.', 'x')
                    BW2GPCRDB[BW] = GPCRDB
                    GPCRDB2BW[GPCRDB] = BW

        #print (BW2GPCRDB)
        #print (GPCRDB2SEQ)
        #print (positions)
        seq_positions = []
        bw_positions = []
        for BW in positions:
            if BW in BW2GPCRDB:
                GPCRDB = BW2GPCRDB[BW]
                if GPCRDB in GPCRDB2SEQ:
                    SEQ = GPCRDB2SEQ[GPCRDB]
                    seq_positions.append(int(SEQ))
                    bw_positions.append(BW)
            #else:
            #   print (BW)

        #print (list(set(seq_positions)))
        #seq_positions = list(set(seq_positions))
        #print (seq_positions)
        #print (bw_positions)
        ## Insert mutation
        mutation_position = gpcr_given.split('_')[-1]
        if mutation_position not in seq_positions:
            if mutation_position != 'WT':
                SEQ = int(mutation_position[1:-1])
                if SEQ in SEQ2GPCRDB:
                    GPCRDB = SEQ2GPCRDB[SEQ]
                    if GPCRDB in GPCRDB2BW:
                        BW = GPCRDB2BW[GPCRDB]
                        seq_positions.append(int(SEQ))
                        bw_positions.append(BW)
                    else:
                        seq_positions.append(int(SEQ))
                        bw_positions.append('-')

        seq_positions = np.array(seq_positions)
        bw_positions = np.array(bw_positions)
        indexes = np.unique(seq_positions, return_index=True)[1]
        seq_positions = seq_positions[indexes]
        #print(bw_positions[indexes])
        bw_positions = bw_positions[indexes]
        #print (seq_positions)
        #print (bw_positions)
        '''
        indexes = np.argsort(seq_positions)
        print (indexes)
        seq_positions = seq_positions[indexes]
        print (seq_positions)
        '''
        seq_positions = seq_positions.tolist()
        bw_positions = bw_positions.tolist()

        #print (seq_positions)
        #print (bw_positions)

        '''
        new_seq_positions = []
        new_bw_positions = []
        mutation_position = gpcr_given.split('_')[1]
        if mutation_position != 'WT':
            position = int(mutation_position[1:-1])
            #print (position)
            if position not in seq_positions:
                for i in range(0, len(seq_positions)):
                    #print (seq_positions[i])
                    if seq_positions[i] < position:
                        if i == len(seq_positions)-1:
                            new_seq_positions.append(position)
                            new_bw_positions.append('-')
                        elif seq_positions[i+1] > position:
                            #print (seq_positions[i], position, seq_positions[i+1])
                            new_seq_positions.append(position)
                            new_bw_positions.append('-')
                        else:
                            new_seq_positions.append(seq_positions[i])
                            new_bw_positions.append(bw_positions[i])
                    else:
                        new_seq_positions.append(seq_positions[i])
                        new_bw_positions.append(bw_positions[i])

            seq_positions = new_seq_positions
            bw_positions = new_bw_positions
        '''

        #print (seq_positions)
        #print (bw_positions)

        return jsonify({'fetch_contacts': scores,
                        'seq_positions': seq_positions,
                        'bw_positions': bw_positions,
                        'sequence': fasta_sequence,
                        'gpcr_name': '_'.join(gpcr_given.split('_')[:-1]) + '/' + gpcr_given.split('_')[-1]
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")

# Function to convert mutation position and given positions
# from BW to PDB (based on PDB ID)
@app.route('/convertPositionsBW2PDB', methods=['GET', 'POST'])
def convertPositionsBW2PDB():
    if request.method == 'POST':
        data = request.get_json(force=True)
        pdbID = data['pdbID']
        positions = data['positions']
        pair_positions = data['pair_positions']
        num_contacts = data['num_contacts']
        AMpositions = data['AMpositions']
        AMpair_positions = data['AMpair_positions']
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        #print (pdbID)
        #print (pair_positions)

        GPCRDB2PDB = {}
        for line in open(path + '/data/PDB/GPCRDB/'+pdbID+'.txt', 'r'):
            GPCRDB2PDB[int(line.split('\t')[3].replace('\n', ''))] = int(line.split('\t')[2])
            bestHIT_ACC = line.replace('\n', '').split('\t')[4].split('|')[1]

        BW2GPCRDB = {}
        GPCRDB2BW = {}
        for line in open(path + '/data/GPCRDB/GPCRDB.tsv', 'r'):
            if 'Name' not in line.split('\t')[0]:
                acc = line.split('\t')[0].split('_')[1]
                if acc == bestHIT_ACC:
                    GPCRDB = int(line.split('\t')[1][1:])
                    #BW = line.split('\t')[2]
                    ## convert . to x in GPCRDB numbering
                    BW = line.split('\t')[2].replace('.', 'x')
                    BW2GPCRDB[BW] = GPCRDB
                    GPCRDB2BW[GPCRDB] = BW

        modified_positions = []
        modified_positions_labels = []
        modified_num_contacts = []

        for num, BW in enumerate(positions.split(',')):
            if BW in BW2GPCRDB:
                GPCRDB = BW2GPCRDB[BW]
                if GPCRDB in GPCRDB2PDB:
                    pdbPosition = GPCRDB2PDB[GPCRDB]
                    modified_positions.append(str(pdbPosition))
                    modified_positions_labels.append(str(BW))
                    modified_num_contacts.append(str(num_contacts.split(',')[num]))
        #print (modified_positions_labels)

        modified_pair_positions = []
        if pair_positions.split() != []:
            for pair in pair_positions.split(','):
                BW1 = pair.split(':')[0]
                BW2 = pair.split(':')[1]
                score = pair.split(':')[2]
                if BW1 in BW2GPCRDB and BW2 in BW2GPCRDB:
                    GPCRDB1 = BW2GPCRDB[BW1]
                    GPCRDB2 = BW2GPCRDB[BW2]
                    if GPCRDB1 in GPCRDB2PDB and GPCRDB2 in GPCRDB2PDB:
                        pdbPosition1 = str(GPCRDB2PDB[GPCRDB1])
                        pdbPosition2 = str(GPCRDB2PDB[GPCRDB2])
                        modified_pair_positions.append(pdbPosition1+':'+pdbPosition2+':'+score)

        mutation_position = '-'
        mutation_position_label = '-'
        
        ## Map sequence to GPCRDB
        handle = open(path + "/static/predictor/output/"+uniq_id+"/GPCRDBblast.txt", 'r')
        blast_records = NCBIXML.parse(handle)
        SEQ2GPCRDB = {}
        for blast_record in blast_records:
            #print (blast_record.query)
            if gpcr_given == blast_record.query:
                for alignment in blast_record.alignments:
                    bestHIT = alignment.title.split(' ')[1]
                    #print (bestHIT)
                    if bestHIT_ACC == bestHIT.split('|')[1]:
                        for hsp in alignment.hsps:
                            q_num = 0
                            s_num = 0
                            for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                                if q!='-' and s!='-':
                                    SEQ2GPCRDB[q_num + hsp.query_start] = s_num + hsp.sbjct_start
                                    q_num += 1
                                    s_num += 1
                                elif q!='-':
                                    q_num += 1
                                else:
                                    s_num += 1

                            break
                        break
                break
        
        if gpcr_given[-3:] != '_WT':
            mutation_sequence_position = int(gpcr_given.split('_')[-1][1:-1])
            if mutation_sequence_position in SEQ2GPCRDB:
                mutation_GPCRDB_position = SEQ2GPCRDB[mutation_sequence_position]
                #print (mutation_GPCRDB_position)
                if mutation_GPCRDB_position in GPCRDB2BW:
                    mutation_position_label = GPCRDB2BW[mutation_GPCRDB_position]

            if mutation_GPCRDB_position in GPCRDB2PDB:
                mutation_position = GPCRDB2PDB[mutation_GPCRDB_position]

        #print (AMpositions)
        #print (SEQ2GPCRDB)
        AMmodified_positions = []
        AMmodified_positions_labels = []

        for num, SEQ in enumerate(AMpositions.split(',')):
            SEQ = int(SEQ)
            if SEQ in SEQ2GPCRDB:
                GPCRDB = SEQ2GPCRDB[SEQ]
                if GPCRDB in GPCRDB2PDB:
                    #print (GPCRDB, GPCRDB2BW[GPCRDB])
                    pdbPosition = GPCRDB2PDB[GPCRDB]
                    AMmodified_positions.append(str(pdbPosition))
                    #AMmodified_positions_labels.append(str(SEQ))
                    AMmodified_positions_labels.append(GPCRDB2BW[GPCRDB])
        #print (modified_positions_labels)

        AMmodified_pair_positions = []
        if AMpair_positions.split() != []:
            for pair in AMpair_positions.split(','):
                SEQ1 = int(pair.split(':')[0])
                SEQ2 = int(pair.split(':')[1])
                score = pair.split(':')[2]
                if SEQ1 in SEQ2GPCRDB and SEQ2 in SEQ2GPCRDB:
                    GPCRDB1 = SEQ2GPCRDB[SEQ1]
                    GPCRDB2 = SEQ2GPCRDB[SEQ2]
                    if GPCRDB1 in GPCRDB2PDB and GPCRDB2 in GPCRDB2PDB:
                        pdbPosition1 = str(GPCRDB2PDB[GPCRDB1])
                        pdbPosition2 = str(GPCRDB2PDB[GPCRDB2])
                        AMmodified_pair_positions.append(pdbPosition1+':'+pdbPosition2+':'+score)

        #print (AMmodified_positions)
        #print (AMmodified_pair_positions)

        pdbData = ''
        if 'AF:' in pdbID:
            for line in open(path+'/data/PDB/AlphaFold/'+pdbID+'.pdb', 'r'):
                pdbData += line
            #print (data)

        return jsonify({'modified_positions': '_'.join(modified_positions),
                        'modified_positions_labels': '_'.join(modified_positions_labels),
                        'modified_num_contacts': '_'.join(modified_num_contacts),
                        'modified_pair_positions': '_'.join(modified_pair_positions),
                        'mutation_position': mutation_position,
                        'mutation_position_label': mutation_position_label,
                        'AMmodified_pair_positions': '_'.join(AMmodified_pair_positions),
                        'AMmodified_positions_labels': '_'.join(AMmodified_positions_labels),
                        'AMmodified_positions': '_'.join(AMmodified_positions),
                        'pdbData': pdbData
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")

def extract_attention_contacts(uniq_id, gpcr_given, gprotein_given, AMcutoff):
    AMvector = np.load(path + "/static/predictor/output/"+uniq_id+"/attentions/"+gpcr_given+"_"+gprotein_given+".npy")
    AMpair_positions = []; AMpositions = []
    for i, row in enumerate(AMvector):
        for j, score in enumerate(row):
            if score >= AMcutoff:
                AMpair_positions.append(str(i+1)+':'+str(j+1)+':'+str(score))
                if str(i+1) not in AMpositions:
                    AMpositions.append(str(i+1))
                if str(j+1) not in AMpositions:
                    AMpositions.append(str(j+1))
    #print (AMpositions, AMpair_positions)
    return AMpositions, AMpair_positions

@app.route('/fetchContactsPDBStructure', methods=['GET', 'POST'])
def fetchContactsPDBStructure():
    if request.method == 'POST':
        data = request.get_json(force=True)
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        cutoff = float(data['cutoff'])
        distance = float(data['distance'])
        uniq_id = data['uniq_id']
        AMcutoff = float(data['AMcutoff'])
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff, distance)
        AMpositions, AMpair_positions = extract_attention_contacts(uniq_id, gpcr_given, gprotein_given, AMcutoff)
        ordered_pdbs = reorder_pdbs(uniq_id, gpcr_given, gprotein_given) ## return list of reordered PDB IDs based on GPCR
        #print (ordered_pdbs)
        return jsonify({'try': positions.tolist(),
                        'ordered_pdbs': ordered_pdbs,
                        'positions': ','.join(positions.tolist()),
                        'num_contacts': ','.join(num_contacts),
                        'pair_positions': ','.join(pair_positions),
                        'AMpositions': ','.join(AMpositions),
                        'AMpair_positions': ','.join(AMpair_positions)
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")

## Function to return list of PDB IDs based on GPCR using BLAST
def reorder_pdbs(uniq_id, gpcr, gprotein):
    path_to_fasta = path+"/static/predictor/output/"+uniq_id+"/input.fasta"
    path_to_output = path+"/static/predictor/output/"+uniq_id+"/"

    os.system('blastp -query '+path_to_fasta+' -num_alignments 5000 -db '+ path + '/data/PDB/blastdb/allPDB -out '+path_to_output+'/blastp_output.txt')

    chain_info = {}
    for line in open(path+'/data/PDB/pdblist.txt', 'r'):
        pdbid = line.split(' ')[0]
        gpcr_chain = line.split(' ')[1]
        gprotein_chain = line.split(' ')[2].replace('\n', '')
        barr_chain = line.split(' ')[3].replace('\n', '')

        chain_info[pdbid] = {}
        chain_info[pdbid]['gpcr_chain'] = gpcr_chain
        if gprotein_chain != '-':
            chain_info[pdbid]['gprotein_chain'] = gprotein_chain
        elif barr_chain != '-':
            chain_info[pdbid]['gprotein_chain'] = barr_chain

        '''
        ## Select pdbid depending G-protein or Barrs
        if 'Barr' in gprotein:
            if barr_chain != '-':
                chain_info[pdbid] = {}
                chain_info[pdbid]['gpcr_chain'] = gpcr_chain
                chain_info[pdbid]['gprotein_chain'] = barr_chain
        else:
            if gprotein_chain != '-':
                chain_info[pdbid] = {}
                chain_info[pdbid]['gpcr_chain'] = gpcr_chain
                chain_info[pdbid]['gprotein_chain'] = gprotein_chain
        '''

    dic = {}
    for line in open(path_to_output+'/blastp_output.txt', 'r'):
        if 'Query=' in line:
            #name = line.split('Query=')[1].replace('\n', '').replace(' ', '')
            name = line.split('Query=')[1].replace('\n', '').lstrip().rstrip()
            #print (name)
            dic[name] = []
        elif line[0] == '>':
            #print (line)
            if 'AF:' in line:
                pdbid = line.split('>')[1].split('|')[0].split()[0]
            else:
                pdbid = line.split('>')[1].split('|')[0].split('_')[0].lower().split()[0]
                #print (pdbid)
            if pdbid in chain_info:
                row = []
                #print (chain_info[pdbid])
                row.append(pdbid)
                row.append(chain_info[pdbid]['gpcr_chain'])
                row.append(chain_info[pdbid]['gprotein_chain'])
                #dic[name].append(row)
                dic[name].append(pdbid+'_'+chain_info[pdbid]['gpcr_chain']+'_'+chain_info[pdbid]['gprotein_chain'])
            else:
                print ('not found', pdbid)

    #print ('PDBs', dic[gpcr])
    #print ('PDBs', dic[name])
    return(dic[gpcr])
    #return None

@app.route('/bestGprotein', methods=['GET', 'POST'])
def bestGprotein():
    if request.method == 'POST':
        data = request.get_json(force=True)
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        #print ('Given', gpcr_given)
        gpcr = '_'.join(gpcr_given.split('_')[:-1])
        #print (gpcr, gpcr_given, 'here')
        variant = gpcr_given.split('_')[-1]
        bestGprotein = ''
        colIndex = 2
        for line in open(path + "/static/predictor/output/"+uniq_id+"/out.tsv", 'r'):
            if '#Input' in line:
                header = line.replace('\n', '').replace('#', '').split('\t')
                #print (header)
            elif line[0] != '#':
                row = line.replace('\n','').split('\t')
                if row[0] == gpcr and row[1] == variant:
                    mx = 0
                    bestGprotein = ''
                    colIndex = None
                    for num, (value, head) in enumerate(zip(row[2:], header[2:])):
                        if float(value) > mx:
                            mx = float(value)
                            bestGprotein = head
                            colIndex = num + 2
                    break
        # print (bestGprotein)
        # print (colIndex)
        return jsonify({'bestGprotein': bestGprotein,
                        'colIndex': colIndex
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")

## Route to output page
@app.route('/input', methods=['GET', 'POST'])
def input():
    if request.method == 'POST':
        input = request.form['input']
        input_file = None ## Upload FASTA file
        ## Run the predictor
        try:
            uniq_id, errorCode, flaggedGPCR = precogx.main(15, input, input_file, 'all', app.root_path)
            print (errorCode)
            if errorCode == 1:
                return render_template('error1.html', flaggedGPCR=json.dumps(flaggedGPCR))
                #return render_template('error1.html')
            elif errorCode == 2:
                return render_template('error2.html', flaggedGPCR=json.dumps(flaggedGPCR))
            #uniq_id = 'OXDUB'
            return redirect('/output/'+uniq_id)
        except:
            return render_template('error.html')
    else:
        return ("<html><h3>It was a GET request</h3></html>")

@app.route('/output/<uniq_id>', methods=['GET', 'POST'])
def output(uniq_id):
    #if request.method == 'GET' or request.method == 'POST':
    if request.method == 'GET':
        #print (os.getcwd())
        print ('running', app.root_path)
        #print ('running', app.instance_path)
        path_to_json_output = "/static/predictor/output/"+uniq_id+"/out.json"
        path_to_fasta = path + "/static/predictor/output/"+uniq_id+"/input.fasta"

        ## extract first entry
        with open(path + path_to_json_output) as f:
            d = json.load(f)

        ## Important: here we are passing GPCR_VAR as gpcr name.
        ## In case VAR is IUPHAR, Emax or LogRAi, gpcr name is sent as GPCR_WT
        ## for sequence, contacts, structure and PCA panels
        first_entry = ''
        first_gprotein = ''
        gpcr_list = []
        for key1 in d:
            for num, key2 in enumerate(d[key1]):
                gpcr = key2[0]
                if key2[1] == 'WT' or key2[1] not in ['GtoPdb', 'LogRAi-TGF', 'Emax-GEMTA']:
                    variant = '_' + key2[1]
                else:
                    variant = '_WT'

                if num == 0:
                    mx = 0
                    first_gprotein_index = 2
                    first_entry = gpcr+variant
                    for count, value in enumerate(key2[2:]):
                        if float(value) > mx:
                            mx = float(value)
                            first_gprotein_index = count + 2

                    #print (key1, d[key1])
                    #print (first_gprotein_index)
                    for line in open(path + "/static/predictor/output/"+uniq_id+"/out.tsv", 'r'):
                        if '#Input' in line:
                            header = line.replace('\n', '').replace('#', '').split('\t')
                            #print (header)
                            break
                    #print (header[first_gprotein_index])
                    first_gprotein = header[first_gprotein_index]

                gpcr_list.append(gpcr+variant)
                #break

        #print (first_gprotein_index)
        #path_to_json_output = "/static/predictor/output/"+uniq_id+"/out.json"
        #path_to_fasta = "/static/predictor/output/"+uniq_id+"/input.fasta"
        return render_template('result2.html',
                                path_to_json_output=json.dumps(path_to_json_output),
                                path_to_fasta=json.dumps(path_to_fasta),
                                first_entry=json.dumps(first_entry),
                                first_gprotein=json.dumps(first_gprotein),
                                first_gprotein_index=json.dumps(first_gprotein_index),
                                gpcr_list=json.dumps(gpcr_list),
                                uniq_id=json.dumps(uniq_id))
    else:
        return ("<html><h3>It was a POST request</h3></html>")

# Route to help page
@app.route('/help')
def help():
    return render_template('help.html')


# Route to about page
@app.route('/about')
def about():
    return render_template('about.html')


# Route to about page
@app.route('/faqs')
def faqs():
    return render_template('faqs.html')


if __name__ == '__main__':
    app.run(debug=True)
