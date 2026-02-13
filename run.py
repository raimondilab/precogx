from optparse import Option
from flask import Flask, render_template, request, jsonify, redirect, url_for, send_file, send_from_directory, abort
import os, sys, json
import numpy as np
from sklearn.preprocessing import MinMaxScaler
from Bio import SearchIO
from Bio.Blast import NCBIXML
from matplotlib import cm
import pandas as pd
import pickle
import csv

# sys.path.insert(1, 'static/predictor/')
# from precogxb_app.static.predictor import precogx

app = Flask(__name__)

##
path = os.getcwd()
path = app.root_path
sys.path.insert(1, path + '/static/predictor/')
import precogx

PROTEIN_NAMES = {}


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
                # value = int(position.split('.')[1])
                value = int(position.split('x')[1]) * (10 ** (len(position.split('x')[1]) * (-1)))
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
        # print (position, base, order, value)
        row = []
        row.append(position)
        row.append(order)
        row.append(base)
        row.append(value)
        data.append(row)

    df = pd.DataFrame(data, columns=['position', 'order', 'base', 'value'])
    df = df.sort_values(['order', 'value'], ascending=[True, True])
    # print (df.to_numpy())
    positions = []
    for row in df.to_numpy():
        positions.append(row[0])

    return (np.array(positions))


def extract_contacts(gprotein_given, cutoff, distance):
    assay = ''
    for line in open(path + '/data/contacts/gprotein_best_layer.txt', 'r'):
        if gprotein_given == line.split('\t')[0]:
            assay = line.split('\t')[1]
    # print (gprotein_given)
    dic = {};
    positions = [];
    pair_positions = [];
    scores = [];
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
    # for line in open(path+'/data/contacts/position_'+assay+'_scaled_web_new.txt', 'r'):
    for line in open(path + '/data/contacts/position_' + assay + '_scaled_web_new2.txt', 'r'):
        gprotein_found = line.split('\t')[0]
        if gprotein_given == gprotein_found:
            # print ('here')
            # print (line)
            pos1 = line.split('\t')[1]
            pos2 = line.split('\t')[2]
            score = float(line.split('\t')[3])
            dis = line.replace('\n', '').split('\t')[-1]
            if dis == '-':
                dis = 1000.0  ## set by default a high value so that it is anyway selected
            else:
                dis = float(dis)
            if score >= cutoff or score <= (-1.0) * cutoff and dis >= distance:
                if pos1 not in dic:
                    dic[pos1] = {}
                dic[pos1][pos2] = score

                if pos2 not in dic:
                    dic[pos2] = {}
                dic[pos2][pos1] = score

                if pos1 not in positions:
                    positions.append(pos1)
                # positions.append(pos1)
                if pos2 not in positions:
                    positions.append(pos2)
                pair_positions.append(pos1 + ':' + pos2 + ':' + str(score))
            scores.append(score)

    # print ('positions', len(positions))
    # print ('distance', distance, dis)
    scoresMax = max(scores)
    scoresMin = min(scores)
    positions = np.array(positions)
    positions = sortPositions(positions)
    # print ('----------------')
    # print (positions)
    # print ('----------------')
    # sys.exit()
    # positions = list(set(positions))
    # positions = np.array(np.sort(positions))
    data = []
    num_contacts = []
    for pos1 in positions:
        row = []
        for pos2 in positions:
            if pos2 in dic[pos1]:
                row.append(dic[pos1][pos2])
            else:
                # row.append(0)
                row.append(None)
        data.append(row)
        num_contacts.append([round(len(dic[pos1]), 2)])
    # print (cutoff, gprotein_given)
    # print ('positions', positions)

    if num_contacts != []:
        scaler = MinMaxScaler(feature_range=(0.35, 1.0))
        num_contacts = scaler.fit_transform(num_contacts)
        num_contacts = num_contacts.flatten().tolist()
        for i in range(0, len(num_contacts)):
            num_contacts[i] = str(num_contacts[i])

    # print (num_contacts)

    return scoresMax, scoresMin, data, positions, pair_positions, num_contacts


@app.route('/fetchAttentionMap', methods=['GET', 'POST'])
def fetchAttentionMap():
    if request.method == 'POST':
        data = request.get_json(force=True)
        # print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        # scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff, distance)
        Xtest = np.load(
            path + '/static/predictor/output/' + uniq_id + '/attentions/' + gpcr_given + '_' + gprotein_given + '.npy')
        seqPositions = [str(i) for i in range(1, len(Xtest[0]) + 1)]
        # print (seqPositions)
        # return jsonify({'fetch_contactsMin': scoresMin, 'fetch_contactsMax': scoresMax, 'fetch_contacts': scores, 'positions': positions.tolist()})
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
        # print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        cutoff = float(data['cutoff'])
        distance = float(data['distance'])
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff,
                                                                                                 distance)
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
    Xs_train_pca = np.load(path + '/static/pca_all/' + pca_type + '.npy', allow_pickle=True)
    # Xs_train_pca = np.load(path+'/static/best_PCA/GNAZ.npy', allow_pickle=True)
    score_coupling, score_uncoupling, Xs_train_pca_coupling, Xs_train_pca_uncoupling, Xs_train_pca_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = filter_gpcr_list(
        Xs_train_pca, assay, gprotein)
    # print ('train', Xs_train_pca_coupling)
    score_coupling = score_coupling.tolist()
    score_uncoupling = score_uncoupling.tolist()
    x_train_coupling = Xs_train_pca_coupling[:, 0].tolist()
    x_train_uncoupling = Xs_train_pca_uncoupling[:, 0].tolist()
    x_train_grey = Xs_train_pca_grey[:, 0].tolist()
    y_train_coupling = Xs_train_pca_coupling[:, 1].tolist()
    y_train_uncoupling = Xs_train_pca_uncoupling[:, 1].tolist()
    y_train_grey = Xs_train_pca_grey[:, 1].tolist()
    return score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey


def filter_gpcr_list(X, assay, gprotein):
    genes_to_consider_coupling = []
    score_coupling = []
    score_uncoupling = []
    genes_to_consider_uncoupling = []
    # print (assay)
    # assay = 'ebBRET'
    if assay == 'TGF':
        num = -1
        for line in open(path + '/data/shedding.tsv', 'r'):
            if line[0] != '#':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                id = line.split('\t')[2]
                name = line.split('\t')[3]
                # Make the range -1 to +1
                score = float(line.split('\t')[num + 4]) + 1
                # print (score)
                # color = cm.get_cmap('RdYlGn', 100)
                # r,g,b,a = color(score)
                # print (score,r,g,b)
                # if score >= -1.0:
                # Move the center (-1.0) to 0.0
                if score >= 0.0:
                    genes_to_consider_coupling.append(gene + '|' + acc)
                    color = cm.get_cmap('Greens', 100)
                    r, g, b, a = color(score)
                    score_coupling.append('rgb(' + str(r) + ',' + str(g) + ',' + str(b) + ')')
                else:
                    genes_to_consider_uncoupling.append(gene + '|' + acc)
                    color = cm.get_cmap('Greys', 100)
                    score *= (-1.0)
                    # To bring the score from range 0 to 1
                    # to 0.25 to 0.75 so it is neither too white (low value)
                    # not too black (high value)
                    score = score / 2 + 0.25
                    r, g, b, a = color(score)
                    score_uncoupling.append('rgb(' + str(r) + ',' + str(g) + ',' + str(b) + ')')
            else:
                flag = 0
                # Replace GNAO1 by GoA
                header = line.replace('\n', '').replace('GNAO1', 'GoA').split('\t')[4:]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        flag = 1
                        break
                # print (num, gprot)
                # print (header)
                if flag == 0:
                    num = -1

    elif assay == 'GEMTA':
        num = 0
        for line in open(path + '/data/ebbret.tsv', 'r', encoding="utf-8"):
            if line[0] != '#':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                id = line.split('\t')[2]
                name = line.split('\t')[3]
                score = float(line.split('\t')[num + 4])

                if score > 0.0:
                    genes_to_consider_coupling.append(gene + '|' + acc)
                    color = cm.get_cmap('Greens', 100)
                    r, g, b, a = color(score)
                    score_coupling.append('rgb(' + str(r) + ',' + str(g) + ',' + str(b) + ')')
                else:
                    genes_to_consider_uncoupling.append(gene + '|' + acc)
                    # score_uncoupling.append('rgb('+str(r)+','+str(g)+','+str(b)+')')
                    score_uncoupling.append('grey')
            else:
                header = line.replace('\n', '').split('\t')[4:]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        # print (gprot)
                        break

    elif assay == 'GtoPdb':
        iuphar_map = {
            'GNAS': 'Gs', 'GNAL': 'Gs',
            'GNAI1': 'Gi/Go', 'GNAI2': 'Gi/Go', 'GNAI3': 'Gi/Go', 'GNAO1': 'Gi/Go', 'GNAZ': 'Gi/Go', 'GoA': 'Gi/Go',
            'GoB': 'Gi/Go',
            'GNA12': 'G12/G13', 'GNA13': 'G12/G13',
            'GNAQ': 'Gq/G11', 'GNA11': 'Gq/G11', 'GNA14': 'Gq/G11', 'GNA15': 'Gq/G11'
        }
        gprotein_fam = iuphar_map[gprotein]
        # print (gprotein_fam)
        for line in open(path + '/data/iuphar.tsv', 'r'):
            if line[0] != '#' and line.split('\t')[1] != '':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                id = line.split('\t')[2]
                name = line.split('\t')[3]
                if gprotein_fam in line:
                    genes_to_consider_coupling.append(gene + '|' + acc)
                    if gprotein_fam in line.split('\t')[4]:
                        score_coupling.append('forestgreen')
                    else:
                        score_coupling.append('lightgreen')
                else:
                    genes_to_consider_uncoupling.append(gene + '|' + acc)
                    score_uncoupling.append('grey')

    elif assay == 'STRING':
        string_map = {
            'Barr1-GRK2': 'ARRB1',
            'Barr2': 'ARRB2',
            'Barr2-GRK2': 'ARRB2'
        }
        barr = string_map[gprotein]
        num = -1
        for line in open(path + '/data/string.tsv', 'r', encoding="utf-8"):
            gene = line.split('\t')[0]
            acc = line.split('\t')[1]
            id = line.split('\t')[2]
            name = line.split('\t')[3]
            if barr in line:
                genes_to_consider_coupling.append(gene + '|' + acc)
                score_coupling.append('green')
            else:
                genes_to_consider_uncoupling.append(gene + '|' + acc)
                score_uncoupling.append('green')
    # print (genes_to_consider_coupling)
    # print (genes_to_consider_uncoupling)

    gpcr_list = []
    for line in open(path + '/static/pca_all/gpcr_list_unscaled_all_layer_GN.txt', 'r'):
        gene = line.replace('\n', '').split('\t')[1]
        acc = line.replace('\n', '').split('\t')[0]
        gpcr_list.append(gene + '|' + acc)
        # gpcr_list.append(line.replace('\n', '').split('\t')[1])

    X_pos = []
    X_neg = []
    X_grey = []
    genes_to_consider_grey = []
    for gene, row in zip(gpcr_list, X):
        # if 'MC1R' in gene:
        #    print ('couple', gene, row[:2])
        if gene in genes_to_consider_coupling:
            X_pos.append(row)
        elif gene in genes_to_consider_uncoupling:
            X_neg.append(row)
        else:
            X_grey.append(row)
            genes_to_consider_grey.append(gene)

    # print (X_pos)

    return (np.array(score_coupling), np.array(score_uncoupling), np.array(X_pos), np.array(X_neg), np.array(X_grey),
            genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey)


@app.route('/fetchPCA', methods=['GET', 'POST'])
def fetchPCA():
    if request.method == 'POST':

        data = request.get_json(force=True)
        print(data)
        assay_given = data['assay']
        pca_type = data['pca_type']
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        # print (gprotein_given, gpcr_given)
        uniq_id = data['uniq_id']
        display_option = int(data['displayPCAOption'])
        # print (display_option)

        # if assay == '':
        assay = '';
        assayList = []
        gemta = ['GNAS', 'GNAI1', 'GNAI2', 'GoB', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA11', 'GNA14', 'GNA15',
                 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2']
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

        gpcrName = '_'.join(gpcr_given.split('_')[:-1])  ## Exclude _WT or _MUT
        variantsToConsider = [];
        x_test = [];
        y_test = [];
        test_names = [];

        if display_option == 1:
            '''
            Option 1: Fetch all point mutations of a given GPCR
            '''
            if gpcr_given[:-3] != '_WT':
                for files in os.listdir(path + '/static/predictor/output/' + uniq_id + '/PCA/'):
                    if files.endswith('.npy') and 'layer' in files and files.split('_')[-1].split('.')[0] != 'WT':
                        if int(pca_type) == int(files.split('layer')[0]):
                            if gpcrName == '_'.join(files.split('_')[1:-1]):
                                variantsToConsider.append(files)
            '''
            End of Option 1
            '''
        else:
            '''
            Option 2: Fetch all input sequences
            '''
            for files in os.listdir(path + '/static/predictor/output/' + uniq_id + '/PCA/'):
                if files.endswith('.npy') and 'layer' in files:
                    if int(pca_type) == int(files.split('layer')[0]):
                        if '_'.join(files.split('_')[1:]).split('.')[
                            0] != gpcrName + '_WT':  # Ignore the WT of the selected coz it is cover below
                            variantsToConsider.append(files)
            '''
            End of Option2
            '''

        for files in variantsToConsider:
            Xs_test_pca = np.load(path + '/static/predictor/output/' + uniq_id + '/PCA/' + files, allow_pickle=True)
            # print ('test',Xs_test_pca)
            x_test.append(Xs_test_pca[0].tolist())
            y_test.append(Xs_test_pca[1].tolist())
            # test_names.append(files.split('_')[1]+'_'+files.split('_')[-1].split('.')[0])
            test_names.append('_'.join(files.split('_')[1:-1]) + '/' + files.split('_')[-1].split('.')[0])

        ### WT
        wt = '_'.join(gpcr_given.split('_')[:-1]) + '_WT'
        wt_name = '_'.join(gpcr_given.split('_')[:-1]) + '/WT'

        Xs_wt_pca = np.load(path + '/static/predictor/output/' + uniq_id + '/PCA/' + pca_type + 'layer_' + wt + '.npy',
                            allow_pickle=True)
        # print (Xs_wt_pca)
        # x_test = Xs_test_pca[:,0].tolist()
        # y_test = Xs_test_pca[:,1].tolist()
        x_wt = Xs_wt_pca[0].tolist()
        y_wt = Xs_wt_pca[1].tolist()
        # print (x_wt)
        # print (y_wt)

        score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = extract_pca(
            gprotein_given, assay, pca_type)
        # print (x_train, y_train, x_test, y_test)
        # print (assay,genes_to_consider_coupling)
        minX = min(x_train_coupling + x_train_uncoupling + x_train_grey)
        maxX = max(x_train_coupling + x_train_uncoupling + x_train_grey)
        minY = min(y_train_coupling + y_train_uncoupling + x_train_grey)
        maxY = max(y_train_coupling + y_train_uncoupling + x_train_grey)
        # print(minY, maxY)
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
        # print (gprotein_given, gpcr_given)
        uniq_id = data['uniq_id']
        display_option = int(data['displayPCAOption'])

        # if assay == '':
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
        gpcrName = '_'.join(gpcr_given.split('_')[:-1])  ## Exclude _WT or _MUT
        variantsToConsider = [];
        x_test = [];
        y_test = [];
        test_names = [];

        if display_option == 1:
            '''
            Option 1: Fetch all point mutations of a given GPCR
            '''
            if gpcr_given[:-3] != '_WT':
                for files in os.listdir(path + '/static/predictor/output/' + uniq_id + '/PCA/'):
                    if files.endswith('.npy') and 'layer' in files and files.split('_')[-1].split('.')[0] != 'WT':
                        if int(pca_type) == int(files.split('layer')[0]):
                            if gpcrName == '_'.join(files.split('_')[1:-1]):
                                variantsToConsider.append(files)
            '''
            End of Option 1
            '''
        else:
            '''
            Option 2: Fetch all input sequences
            '''
            for files in os.listdir(path + '/static/predictor/output/' + uniq_id + '/PCA/'):
                if files.endswith('.npy') and 'layer' in files:
                    if int(pca_type) == int(files.split('layer')[0]):
                        if '_'.join(files.split('_')[1:]).split('.')[
                            0] != gpcrName + '_WT':  # Ignore the WT of the selected coz it is cover below
                            variantsToConsider.append(files)
            '''
            End of Option2
            '''

        for files in variantsToConsider:
            Xs_test_pca = np.load(path + '/static/predictor/output/' + uniq_id + '/PCA/' + files, allow_pickle=True)
            # print ('test',Xs_test_pca)
            x_test.append(Xs_test_pca[0].tolist())
            y_test.append(Xs_test_pca[1].tolist())
            # test_names.append(files.split('_')[1]+'_'+files.split('_')[-1].split('.')[0])
            test_names.append('_'.join(files.split('_')[1:-1]) + '/' + files.split('_')[-1].split('.')[0])

        ### WT
        wt = '_'.join(gpcr_given.split('_')[:-1]) + '_WT'
        wt_name = '_'.join(gpcr_given.split('_')[:-1]) + '/WT'
        # print (test_names)
        '''
        if pca_type == 'GPCRome':
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/33layer_'+wt+'.npy', allow_pickle=True)
        elif pca_type == 'Best PCA':
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+wt+'.npy', allow_pickle=True)
        else:
        '''
        Xs_wt_pca = np.load(path + '/static/predictor/output/' + uniq_id + '/PCA/' + pca_type + 'layer_' + wt + '.npy',
                            allow_pickle=True)
        # print (Xs_wt_pca)
        # x_test = Xs_test_pca[:,0].tolist()
        # y_test = Xs_test_pca[:,1].tolist()
        x_wt = Xs_wt_pca[0].tolist()
        y_wt = Xs_wt_pca[1].tolist()
        # print (x_wt)
        # print (y_wt)

        # score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = extract_pca(gprotein_given, assay, pca_type)
        '''
        if pca_type == 'Best PCA':
            Xs_train_pca = np.load(path+'/static/best_PCA/'+gprotein_given+'.npy', allow_pickle=True)
        elif pca_type == 'GPCRome':
            Xs_train_pca = np.load(path+'/static/33layer_PCA/33layer.npy', allow_pickle=True)
        else:
            '''
        Xs_train_pca = np.load(path + '/static/pca_all/' + pca_type + '.npy', allow_pickle=True)

        classes = {}
        # for line in open(path+'/data/classification.txt', 'r'):
        for line in open(path + '/data/classification2.txt', 'r'):
            if 'Uniprot_acc' not in line:
                acc = line.split('\t')[0]
                cls = line.split('\t')[-1].replace('\n', '')
                classes[acc] = cls

        gpcr_list = []
        for line in open(path + '/static/pca_all/gpcr_list_unscaled_all_layer_GN.txt', 'r'):
            gene = line.replace('\n', '').split('\t')[1]
            acc = line.replace('\n', '').split('\t')[0]
            gpcr_list.append(gene + '|' + acc)
            # gpcr_list.append(line.replace('\n', '').split('\t')[1])

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
            '''
            if 'MC1R' in gpcr:
                print ('class', gpcr, row[:2])
            '''
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
        # print (X_classA)

        x_classA = X_classA[:, 0].tolist()
        y_classA = X_classA[:, 1].tolist()
        x_classB1 = X_classB1[:, 0].tolist()
        y_classB1 = X_classB1[:, 1].tolist()
        x_classB2 = X_classB2[:, 0].tolist()
        y_classB2 = X_classB2[:, 1].tolist()
        x_classC = X_classC[:, 0].tolist()
        y_classC = X_classC[:, 1].tolist()
        x_frizzeled = X_frizzeled[:, 0].tolist()
        y_frizzeled = X_frizzeled[:, 1].tolist()
        x_taste = X_taste[:, 0].tolist()
        y_taste = X_taste[:, 1].tolist()
        x_other = X_other[:, 0].tolist()
        y_other = X_other[:, 1].tolist()

        '''
        for gpcr, x, y in zip(classA, x_classA, y_classA):
            if 'MC1R' in gpcr:
                print ('classA', gpcr, x, y)
                print (len(x_classA), len(y_classA))
        '''

        # print (y_classA)
        # print (genes_to_consider_coupling)
        minX = min(x_classA + x_classB1 + x_classB2 + x_classC + x_frizzeled + x_taste + x_other)
        maxX = max(x_classA + x_classB1 + x_classB2 + x_classC + x_frizzeled + x_taste + x_other)
        minY = min(y_classA + y_classB1 + y_classB2 + y_classC + y_frizzeled + y_taste + y_other)
        maxY = max(y_classA + y_classB1 + y_classB2 + y_classC + y_frizzeled + y_taste + y_other)
        # print(minX, maxX)

        # print (x_test, y_test)
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
        # print (gprotein_given, gpcr_given)
        uniq_id = data['uniq_id']

        for line in open(path + '/data/contacts/gprotein_best_layer.txt', 'r'):
            if gprotein_given == line.split('\t')[0]:
                assay = line.split('\t')[1]
                bestPCA = line.split('\t')[2].replace('\n', '')
                break

        layers = [str(i) for i in range(0, 34)]
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
    handle = open(path + "/static/predictor/output/" + uniq_id + "/GPCRDBblast.txt", 'r')
    blast_records = NCBIXML.parse(handle)
    # print (blast_records)

    GPCRDB2SEQ = {}
    SEQ2GPCRDB = {}
    for blast_record in blast_records:
        # print (blast_record.query)
        if gpcr_given == blast_record.query:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    bestHIT = alignment.title.split(' ')[1]
                    q_num = 0
                    s_num = 0
                    for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                        if q != '-' and s != '-':
                            GPCRDB2SEQ[s_num + hsp.sbjct_start] = q_num + hsp.query_start
                            SEQ2GPCRDB[q_num + hsp.query_start] = s_num + hsp.sbjct_start
                            q_num += 1
                            s_num += 1
                        elif q != '-':
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
            # print (bestHIT)
            # print (GPCRDB2SEQ)
            break
    return (GPCRDB2SEQ, SEQ2GPCRDB, bestHIT)


@app.route('/fetchContactsSequence', methods=['GET', 'POST'])
def fetchContactsSequence():
    if request.method == 'POST':
        data = request.get_json(force=True)
        # print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        # print (gpcr_given)
        # print (gpcr_given, gpcr_given.split('_')[1], 'sequence')
        path_to_fasta = data['path_to_fasta']
        uniq_id = data['uniq_id']
        cutoff = float(data['cutoff'])
        distance = float(data['distance'])
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff,
                                                                                                 distance)

        # print ('fetch_seq', positions)
        fasta_sequence = '';
        flag = 0
        for line in open(path_to_fasta):
            if line[0] == '>':
                if flag == 1:
                    break
                flag = 0
                # gpcr_found = line.split('>')[1].replace('\n', '').replace(' ', '')
                gpcr_found = line.split('>')[1].replace('\n', '').lstrip().rstrip()
                if gpcr_found == gpcr_given:
                    flag = 1
            elif flag == 1:
                fasta_sequence += line.replace('\n', '')

        GPCRDB2SEQ, SEQ2GPCRDB, bestHIT = DoBLAST(uniq_id, gpcr_given)
        bestHIT_ACC = bestHIT.split('|')[1]
        # print (bestHIT)

        BW2GPCRDB = {}
        GPCRDB2BW = {}
        for line in open(path + '/data/GPCRDB/GPCRDB.tsv', 'r'):
            if 'Name' not in line.split('\t')[0]:
                acc = line.split('\t')[0].split('_')[1]
                if acc == bestHIT_ACC:
                    GPCRDB = int(line.split('\t')[1][1:])
                    # BW = line.split('\t')[2]
                    ## convert . to x in GPCRDB numbering
                    BW = line.split('\t')[2].replace('.', 'x')
                    BW2GPCRDB[BW] = GPCRDB
                    GPCRDB2BW[GPCRDB] = BW

        # print (BW2GPCRDB)
        # print (GPCRDB2SEQ)
        # print (positions)
        seq_positions = []
        bw_positions = []
        for BW in positions:
            if BW in BW2GPCRDB:
                GPCRDB = BW2GPCRDB[BW]
                if GPCRDB in GPCRDB2SEQ:
                    SEQ = GPCRDB2SEQ[GPCRDB]
                    seq_positions.append(int(SEQ))
                    bw_positions.append(BW)
            # else:
            #   print (BW)

        # print (list(set(seq_positions)))
        # seq_positions = list(set(seq_positions))
        # print (seq_positions)
        # print (bw_positions)
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
        # print(bw_positions[indexes])
        bw_positions = bw_positions[indexes]
        # print (seq_positions)
        # print (bw_positions)
        '''
        indexes = np.argsort(seq_positions)
        print (indexes)
        seq_positions = seq_positions[indexes]
        print (seq_positions)
        '''
        seq_positions = seq_positions.tolist()
        bw_positions = bw_positions.tolist()

        # print (seq_positions)
        # print (bw_positions)

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

        # print (seq_positions)
        # print (bw_positions)

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
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        # print (pdbID)
        # print (pair_positions)

        GPCRDB2PDB = {}
        for line in open(path + '/data/PDB/GPCRDB/' + pdbID + '.txt', 'r'):
            GPCRDB2PDB[int(line.split('\t')[3].replace('\n', ''))] = int(line.split('\t')[2])
            bestHIT_ACC = line.replace('\n', '').split('\t')[4].split('|')[1]

        BW2GPCRDB = {}
        GPCRDB2BW = {}
        for line in open(path + '/data/GPCRDB/GPCRDB.tsv', 'r'):
            if 'Name' not in line.split('\t')[0]:
                acc = line.split('\t')[0].split('_')[1]
                if acc == bestHIT_ACC:
                    GPCRDB = int(line.split('\t')[1][1:])
                    # BW = line.split('\t')[2]
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
        # print (modified_positions_labels)

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
                        modified_pair_positions.append(pdbPosition1 + ':' + pdbPosition2 + ':' + score)

        mutation_position = '-'
        mutation_position_label = '-'
        if '_WT' not in gpcr_given:
            mutation_sequence_position = int(gpcr_given.split('_')[-1][1:-1])
            handle = open(path + "/static/predictor/output/" + uniq_id + "/GPCRDBblast.txt", 'r')
            blast_records = NCBIXML.parse(handle)
            SEQ2GPCRDB = {}
            for blast_record in blast_records:
                # print (blast_record.query)
                if gpcr_given == blast_record.query:
                    for alignment in blast_record.alignments:
                        bestHIT = alignment.title.split(' ')[1]
                        # print (bestHIT)
                        if bestHIT_ACC == bestHIT.split('|')[1]:
                            for hsp in alignment.hsps:
                                q_num = 0
                                s_num = 0
                                for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                                    if q != '-' and s != '-':
                                        SEQ2GPCRDB[q_num + hsp.query_start] = s_num + hsp.sbjct_start
                                        q_num += 1
                                        s_num += 1
                                    elif q != '-':
                                        q_num += 1
                                    else:
                                        s_num += 1

                                break
                            break
                    break

            if mutation_sequence_position in SEQ2GPCRDB:
                mutation_GPCRDB_position = SEQ2GPCRDB[mutation_sequence_position]
                # print (mutation_GPCRDB_position)
                if mutation_GPCRDB_position in GPCRDB2BW:
                    mutation_position_label = GPCRDB2BW[mutation_GPCRDB_position]

            if mutation_GPCRDB_position in GPCRDB2PDB:
                mutation_position = GPCRDB2PDB[mutation_GPCRDB_position]

        # print (modified_positions)
        # print (modified_pair_positions)

        # print (mutation_position)
        # print (mutation_position_label)

        pdbData = ''
        if 'AF:' in pdbID:
            for line in open(path + '/data/PDB/AlphaFold/' + pdbID + '.pdb', 'r'):
                pdbData += line
            # print (data)

        return jsonify({'modified_positions': '_'.join(modified_positions),
                        'modified_positions_labels': '_'.join(modified_positions_labels),
                        'modified_num_contacts': '_'.join(modified_num_contacts),
                        'modified_pair_positions': '_'.join(modified_pair_positions),
                        'mutation_position': mutation_position,
                        'mutation_position_label': mutation_position_label,
                        'pdbData': pdbData
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")


def extract_attention_contacts(uniq_id, gpcr_given, gprotein_given, AMcutoff):
    AMvector = np.load(
        path + "/static/predictor/output/" + uniq_id + "/attentions/" + gpcr_given + "_" + gprotein_given + ".npy")
    AMpair_positions = [];
    AMpositions = []
    for i, row in enumerate(AMvector):
        for j, score in enumerate(row):
            if score >= AMcutoff:
                AMpair_positions.append(str(i + 1) + ':' + str(j + 1) + ':' + str(score))
                if str(i + 1) not in AMpositions:
                    AMpositions.append(str(i + 1))
                if str(j + 1) not in AMpositions:
                    AMpositions.append(str(j + 1))
    # print (AMpositions, AMpair_positions)
    return AMpositions, AMpair_positions


@app.route('/fetchContactsPDBStructure', methods=['GET', 'POST'])
def fetchContactsPDBStructure():
    if request.method == 'POST':
        data = request.get_json(force=True)
        if data:
            gprotein_given = data['gprotein']
            gpcr_given = data['gpcr']
            cutoff = float(data['cutoff'])
            distance = float(data['distance'])
            uniq_id = data['uniq_id']
            scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given,
                                                                                                     cutoff, distance)
            ordered_pdbs = reorder_pdbs(uniq_id, gpcr_given,
                                        gprotein_given)  ## return list of reordered PDB IDs based on GPCR
            # print (ordered_pdbs)
            return jsonify({'try': positions.tolist(),
                            'ordered_pdbs': ordered_pdbs,
                            'positions': ','.join(positions.tolist()),
                            'num_contacts': ','.join(num_contacts),
                            'pair_positions': ','.join(pair_positions)})
    else:
        return ("<html><h3>It was a GET request</h3></html>")


## Function to return list of PDB IDs based on GPCR using BLAST
def reorder_pdbs(uniq_id, gpcr, gprotein):
    path_to_fasta = path + "/static/predictor/output/" + uniq_id + "/input.fasta"
    path_to_output = path + "/static/predictor/output/" + uniq_id + "/"

    os.system(
        'blastp -query ' + path_to_fasta + ' -num_alignments 5000 -db ' + path + '/data/PDB/blastdb/allPDB -out ' + path_to_output + '/blastp_output.txt')

    chain_info = {}
    for line in open(path + '/data/PDB/pdblist.txt', 'r'):
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
    for line in open(path_to_output + '/blastp_output.txt', 'r'):
        if 'Query=' in line:
            # name = line.split('Query=')[1].replace('\n', '').replace(' ', '')
            name = line.split('Query=')[1].replace('\n', '').lstrip().rstrip()
            # print (name)
            dic[name] = []
        elif line[0] == '>':
            # print (line)
            if 'AF:' in line:
                pdbid = line.split('>')[1].split('|')[0].split()[0]
            else:
                pdbid = line.split('>')[1].split('|')[0].split('_')[0].lower().split()[0]
                # print (pdbid)
            if pdbid in chain_info:
                row = []
                # print (chain_info[pdbid])
                row.append(pdbid)
                row.append(chain_info[pdbid]['gpcr_chain'])
                row.append(chain_info[pdbid]['gprotein_chain'])
                # dic[name].append(row)
                dic[name].append(
                    pdbid + '_' + chain_info[pdbid]['gpcr_chain'] + '_' + chain_info[pdbid]['gprotein_chain'])
            else:
                print('not found', pdbid)

    # print ('PDBs', dic[gpcr])
    # print ('PDBs', dic[name])
    if gpcr in dic.keys():
        return (dic[gpcr])
    return None


@app.route('/bestGprotein', methods=['GET', 'POST'])
def bestGprotein():
    if request.method == 'POST':
        data = request.get_json(force=True)
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        # print ('Given', gpcr_given)
        gpcr = '_'.join(gpcr_given.split('_')[:-1])
        # print (gpcr, gpcr_given, 'here')
        variant = gpcr_given.split('_')[-1]
        bestGprotein = ''
        colIndex = 2
        for line in open(path + "/static/predictor/output/" + uniq_id + "/out.tsv", 'r'):
            if '#Input' in line:
                header = line.replace('\n', '').replace('#', '').split('\t')
                # print (header)
            elif line[0] != '#':
                row = line.replace('\n', '').split('\t')
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
        input_file = None  ## Upload FASTA file
        ## Run the predictor
        try:
            uniq_id, errorCode, flaggedGPCR = precogx.main(15, input, input_file, 'all', app.root_path)
            print(errorCode)
            if errorCode == 1:
                return render_template('error1.html', flaggedGPCR=json.dumps(flaggedGPCR))
                # return render_template('error1.html')
            elif errorCode == 2:
                return render_template('error2.html', flaggedGPCR=json.dumps(flaggedGPCR))
            # uniq_id = 'OXDUB'
            return redirect('/output/' + uniq_id)
        except:
            return render_template('error.html')
    else:
        return ("<html><h3>It was a GET request</h3></html>")


@app.route('/output/<uniq_id>', methods=['GET', 'POST'])
def output(uniq_id):
    # if request.method == 'GET' or request.method == 'POST':
    if request.method == 'GET':
        # print (os.getcwd())
        print('running', app.root_path)
        # print ('running', app.instance_path)
        path_to_json_output = "/static/predictor/output/" + uniq_id + "/out.json"
        path_to_fasta = path + "/static/predictor/output/" + uniq_id + "/input.fasta"

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
                    first_entry = gpcr + variant
                    for count, value in enumerate(key2[2:]):
                        if float(value) > mx:
                            mx = float(value)
                            first_gprotein_index = count + 2

                    # print (key1, d[key1])
                    # print (first_gprotein_index)
                    for line in open(path + "/static/predictor/output/" + uniq_id + "/out.tsv", 'r'):
                        if '#Input' in line:
                            header = line.replace('\n', '').replace('#', '').split('\t')
                            # print (header)
                            break
                    # print (header[first_gprotein_index])
                    first_gprotein = header[first_gprotein_index]

                gpcr_list.append(gpcr + variant)
                # break

        # print (first_gprotein_index)
        # path_to_json_output = "/static/predictor/output/"+uniq_id+"/out.json"
        # path_to_fasta = "/static/predictor/output/"+uniq_id+"/input.fasta"
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


PDF_FOLDER = path + '/pdfs'
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
PKL_PATH = os.path.join(BASE_DIR, 'static', 'info.pkl')
CSV_PATH = os.path.join(BASE_DIR, 'static', 'Reported_coupling.csv')

# Cache dos dados
GPCR_METADATA = {}
COUPLING_DATA = {}


def load_data():
    global GPCR_METADATA, COUPLING_DATA, PROTEIN_NAMES

    if os.path.exists(PKL_PATH):
        try:
            with open(PKL_PATH, 'rb') as f:
                GPCR_METADATA = pickle.load(f)

                GPCR_METADATA['TRHR'] = {'acc': 'P34981', 'cn': ' Thyrotropin-releasing hormone receptor '}
                GPCR_METADATA['NTSR1'] = {'acc': 'P30989', 'cn': ' Neurotensin receptor type 1 '}

        except Exception as e:
            print(f"Erro ao carregar pickle: {e}")

    if os.path.exists(CSV_PATH):
        try:
            with open(CSV_PATH, 'r', encoding='utf-8-sig', errors='ignore') as f:
                reader = csv.reader(f)
                for row in reader:
                    if len(row) > 2:
                        receptor = row[0].strip()
                        coupling_info = row[2:-1]
                        prot_name = row[1].strip()

                        if prot_name:
                            PROTEIN_NAMES[receptor] = prot_name

                        formatted_coupling = []
                        familias_g = ["Gs", "Gi/o", "Gq/11", "G12/13"]

                        for i, val in enumerate(coupling_info):
                            if val and val.strip() in ['1', '2']:
                                if i < len(familias_g):
                                    type_c = "Primary" if val.strip() == '1' else "Secondary"
                                    formatted_coupling.append({"family": familias_g[i], "type": type_c})

                        COUPLING_DATA[receptor] = formatted_coupling
        except Exception as e:
            print(f"Loading error: {e}")


load_data()


# Route to shedding page
@app.route('/shedding')
def shedding():
    return render_template('shedding.html')


@app.route('/precog3D')
def precog3d():
    return render_template('precog3D.html')


GPCR_FAMILIES = {
    "5-Hydroxytryptamine receptors": ["HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR4",
                                      "HTR6", "HTR7"],
    "Acetylcholine receptors (Muscarinic)": ["CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"],
    "Adrenoceptors": ["ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3"],
    "Adenosine receptors": ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"],
    "Alpha-Ketoglutarate receptor": ["OXGR1"],
    "Anaphylatoxin chemotactic receptors": ["C3AR1", "C5AR1", "C5AR2"],
    "Angiotensin receptors": ["AGTR1", "AGTR2"],
    "Apelin receptor": ["APLNR"],
    "Bile acid receptor": ["GPBAR1"],
    "Bombesin receptors": ["NMBR", "GRPR", "BB3"],
    "Bradykinin receptors": ["BDKRB1", "BDKRB2"],
    "Calcitonin receptors": ["CALCR", "CALCRL"],
    "Cannabinoid receptors": ["CNR1", "CNR2"],
    "Chemokine receptors": ["CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CCR9", "CCR10", "CXCR1",
                            "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CX3CR1", "XCR1"],
    "Cholecystokinin receptors": ["CCKAR", "CCKBR"],
    "Dopamine receptors": ["DRD1", "DRD2", "DRD3", "DRD4", "DRD5"],
    "Endothelin receptors": ["EDNRA", "EDNRB"],
    "Formylpeptide receptors": ["FPR1", "FPR2", "FPR3"],
    "Free fatty acid receptors": ["FFAR1", "FFAR2", "FFAR3", "FFAR4"],
    "GABA receptors": ["GABBR1", "GABBR2"],
    "Galanin receptors": ["GALR1", "GALR2", "GALR3"],
    "Ghrelin receptor": ["GHSR"],
    "Glucagon receptor family": ["GCGR", "GLP1R", "GLP2R", "GIPR", "GHRHR"],
    "Glycoprotein hormone receptors": ["FSHR", "LHCGR", "TSHR"],
    "Gonadotrophin-releasing hormone receptors": ["GNRHR", "GNRHR2"],
    "Histamine receptors": ["HRH1", "HRH2", "HRH3", "HRH4"],
    "Hydroxycarboxylic acid receptors": ["HCAR1", "HCAR2", "HCAR3"],
    "Kisspeptin receptor": ["KISS1R"],
    "Leukotriene receptors": ["LTB4R", "LTB4R2", "CYSLTR1", "CYSLTR2"],
    "Lysophospholipid (LPA) receptors": ["LPAR1", "LPAR2", "LPAR3", "LPAR4", "LPAR5", "LPAR6"],
    "Lysophospholipid (S1P) receptors": ["S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5"],
    "Melanin-concentrating hormone receptors": ["MCHR1", "MCHR2"],
    "Melanocortin receptors": ["MC1R", "MC2R", "MC3R", "MC4R", "MC5R"],
    "Melatonin receptors": ["MTNR1A", "MTNR1B"],
    "Motilin receptor": ["MLNR"],
    "Neuropeptide FF/neuropeptide AF receptors": ["NPFFR1", "NPFFR2"],
    "Neuropeptide S receptor": ["NPSR1"],
    "Neuropeptide W/neuropeptide B receptors": ["NPBWR1", "NPBWR2"],
    "Neuropeptide Y receptors": ["NPY1R", "NPY2R", "NPY4R", "NPY5R"],
    "Neurotensin receptors": ["NTSR1", "NTSR2"],
    "Opioid receptors": ["OPRD1", "OPRK1", "OPRM1", "OPRL1"],
    "Orexin receptors": ["HCRTR1", "HCRTR2"],
    "P2Y receptors": ["P2RY1", "P2RY2", "P2RY4", "P2RY6", "P2RY10", "P2RY11", "P2RY12", "P2RY13", "P2RY14"],
    "Parathyroid hormone receptors": ["PTH1R", "PTH2R"],
    "Platelet-activating factor receptor": ["PTAFR"],
    "Prokineticin receptors": ["PROKR1", "PROKR2"],
    "Prolactin-releasing peptide receptor": ["PRLHR"],
    "Prostanoid receptors": ["PTGDR", "PTGDR2", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGFR", "PTGIR", "TBXA2R"],
    "Protease-activated receptors": ["F2R", "F2RL1", "F2RL2", "F2RL3"],
    "Relaxin family peptide receptors": ["RXFP1", "RXFP2", "RXFP3", "RXFP4"],
    "Somatostatin receptors": ["SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5"],
    "Tachykinin receptors": ["TACR1", "TACR2", "TACR3"],
    "Thyrotropin-releasing hormone receptors": ["TRHR"],
    "Urotensin receptor": ["UTS2R"],
    "VIP and PACAP receptors": ["ADCYAP1R1", "VIPR1", "VIPR2"],
    "Vasopressin and oxytocin receptors": ["AVPR1A", "AVPR1B", "AVPR2", "OXTR"],
    "Other / Orphan GPCRs": [
        "GPR17", "GPR34", "GPR35", "GPR55", "GPR84",
        "GPR119", "GPR132", "GPR174", "GPR183",
        "MRGPRX1", "MRGPRX2",
        "GPR156", "GPR158", "GPR179", "GPRC5A", "GPRC5B", "GPRC5C", "GPRC5D"
    ],
    "Class Frizzled": ["FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "SMO"]
}

GPCR_FAMILIES_COMPLEX = {
    "5-Hydroxytryptamine receptors": ["HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7"],
    "Acetylcholine receptors (Muscarinic)": ["CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5"],
    "Adrenoceptors": ["ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3"],
    "Adenosine receptors": ["ADORA1", "ADORA2A", "ADORA2B", "ADORA3"],
    "Alpha-Ketoglutarate receptor": ["OXGR1"],
    "Anaphylatoxin chemotactic receptors": ["C3AR1", "C5AR1", "C5AR2"],
    "Angiotensin receptors": ["AGTR1", "AGTR2"],
    "Apelin receptor": ["APLNR"],
    "Bile acid receptor": ["GPBAR1"],
    "Bombesin receptors": ["NMBR", "GRPR", "BRS3"],
    "Bradykinin receptors": ["BDKRB1", "BDKRB2"],
    "Calcitonin receptors": ["CALCR", "CALCRL"],
    "Cannabinoid receptors": ["CNR1", "CNR2"],
    "Chemokine receptors": [
        "CCR1", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CCR9", "CCR10",
        "CXCR1", "CXCR2", "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CX3CR1", "XCR1",
        "ACKR1", "ACKR2", "ACKR3", "ACKR4", "CCRL2"
    ],
    "Cholecystokinin receptors": ["CCKAR", "CCKBR"],
    "Dopamine receptors": ["DRD1", "DRD2", "DRD3", "DRD4", "DRD5"],
    "Endothelin receptors": ["EDNRA", "EDNRB"],
    "Formylpeptide receptors": ["FPR1", "FPR2", "FPR3"],
    "Free fatty acid receptors": ["FFAR1", "FFAR2", "FFAR3", "FFAR4"],
    "GABA receptors": ["GABBR1", "GABBR2"],
    "Galanin receptors": ["GALR1", "GALR2", "GALR3"],
    "Ghrelin receptor": ["GHSR"],
    "Glucagon receptor family": ["GCGR", "GLP1R", "GLP2R", "GIPR", "GHRHR"],
    "Glycoprotein hormone receptors": ["FSHR", "LHCGR", "TSHR"],
    "Gonadotrophin-releasing hormone receptors": ["GNRHR"],
    "Histamine receptors": ["HRH1", "HRH2", "HRH3", "HRH4"],
    "Hydroxycarboxylic acid receptors": ["HCAR1", "HCAR2", "HCAR3"],
    "Kisspeptin receptor": ["KISS1R"],
    "Leukotriene receptors": ["LTB4R", "LTB4R2", "CYSLTR1", "CYSLTR2"],
    "Lysophospholipid (LPA) receptors": ["LPAR1", "LPAR2", "LPAR3", "LPAR4", "LPAR5", "LPAR6"],
    "Lysophospholipid (S1P) receptors": ["S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5"],
    "Melanin-concentrating hormone receptors": ["MCHR1", "MCHR2"],
    "Melanocortin receptors": ["MC1R", "MC2R", "MC3R", "MC4R", "MC5R"],
    "Melatonin receptors": ["MTNR1A", "MTNR1B"],
    "Motilin receptor": ["MLNR"],
    "Neuropeptide FF/neuropeptide AF receptors": ["NPFFR1", "NPFFR2"],
    "Neuropeptide S receptor": ["NPSR1"],
    "Neuropeptide W/neuropeptide B receptors": ["NPBWR1", "NPBWR2"],
    "Neuropeptide Y receptors": ["NPY1R", "NPY2R", "NPY4R", "NPY5R"],
    "Neurotensin receptors": ["NTSR1", "NTSR2"],
    "Opioid receptors": ["OPRD1", "OPRK1", "OPRM1", "OPRL1"],
    "Orexin receptors": ["HCRTR1", "HCRTR2"],
    "P2Y receptors": ["P2RY1", "P2RY2", "P2RY4", "P2RY6", "P2RY8", "P2RY10", "P2RY11", "P2RY12", "P2RY13", "P2RY14"],
    "Parathyroid hormone receptors": ["PTH1R", "PTH2R"],
    "Platelet-activating factor receptor": ["PTAFR"],
    "Prokineticin receptors": ["PROKR1", "PROKR2"],
    "Prolactin-releasing peptide receptor": ["PRLHR"],
    "Prostanoid receptors": ["PTGDR", "PTGDR2", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGFR", "PTGIR", "TBXA2R"],
    "Protease-activated receptors": ["F2R", "F2RL1", "F2RL2", "F2RL3"],
    "Relaxin family peptide receptors": ["RXFP1", "RXFP2", "RXFP3", "RXFP4"],
    "Somatostatin receptors": ["SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5"],
    "Tachykinin receptors": ["TACR1", "TACR2", "TACR3"],
    "Thyrotropin-releasing hormone receptors": ["TRHR"],
    "Urotensin receptor": ["UTS2R"],
    "VIP and PACAP receptors": ["ADCYAP1R1", "VIPR1", "VIPR2"],
    "Vasopressin and oxytocin receptors": ["AVPR1A", "AVPR1B", "AVPR2", "OXTR"],
    "Metabotropic glutamate receptors": ["GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8"],
    "Taste receptors": [
        "TAS1R1", "TAS1R2", "TAS1R3", "TAS2R1", "TAS2R10", "TAS2R13", "TAS2R14", "TAS2R16", "TAS2R19",
        "TAS2R20", "TAS2R3", "TAS2R30", "TAS2R31", "TAS2R38", "TAS2R39", "TAS2R4", "TAS2R40", "TAS2R41",
        "TAS2R42", "TAS2R43", "TAS2R45", "TAS2R46", "TAS2R5", "TAS2R50", "TAS2R60", "TAS2R7", "TAS2R8", "TAS2R9"
    ],
    "Opsins": ["OPN1LW", "OPN1MW", "OPN1SW", "OPN3", "OPN4", "OPN5", "RHO"],
    "Trace amine-associated receptors": ["TAAR1", "TAAR2", "TAAR5", "TAAR6", "TAAR8", "TAAR9"],
    "Chelmerin receptors": ["CMKLR1", "CMKLR2"],
    "Calcium-sensing receptor": ["CASR"],
    "Corticotropin-releasing factor receptors": ["CRHR1", "CRHR2"],
    "Neuromedin U receptors": ["NMUR1", "NMUR2"],
    "Secretin receptor": ["SCTR"],
    "Leucine-rich repeat-containing GPCRs": ["LGR4", "LGR5", "LGR6"],
    "Pyroglutamylated RFamide peptide receptor": ["QRFPR"],
    "Oxoeicosanoid receptor": ["OXER1"],
    "Succinate receptor": ["SUCNR1"],
    "Class Frizzled": ["FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "FZD10", "SMO"],
    "Adhesion GPCRs": [
        "ADGRA1", "ADGRA2", "ADGRA3", "ADGRB1", "ADGRB2", "ADGRB3", "ADGRD1", "ADGRD2",
        "ADGRE1", "ADGRE2", "ADGRE3", "ADGRE5", "ADGRF1", "ADGRF3", "ADGRF4", "ADGRF5",
        "ADGRG1", "ADGRG2", "ADGRG3", "ADGRG4", "ADGRG5", "ADGRG6", "ADGRG7",
        "ADGRL1", "ADGRL2", "ADGRL3", "ADGRL4", "ADGRV1", "CELSR1", "CELSR2", "CELSR3"
    ],
    "Other / Orphan GPCRs": [
        "GPER1", "GPR101", "GPR107", "GPR119", "GPR12", "GPR132", "GPR135", "GPR137", "GPR139", "GPR141",
        "GPR142", "GPR143", "GPR146", "GPR148", "GPR149", "GPR15", "GPR150", "GPR151", "GPR152", "GPR153",
        "GPR156", "GPR157", "GPR158", "GPR160", "GPR161", "GPR162", "GPR17", "GPR171", "GPR173", "GPR174",
        "GPR176", "GPR179", "GPR18", "GPR182", "GPR183", "GPR19", "GPR20", "GPR21", "GPR22", "GPR25",
        "GPR26", "GPR27", "GPR3", "GPR31", "GPR32", "GPR33", "GPR34", "GPR35", "GPR37", "GPR37L1",
        "GPR39", "GPR4", "GPR42", "GPR45", "GPR50", "GPR52", "GPR55", "GPR6", "GPR61", "GPR62",
        "GPR63", "GPR65", "GPR68", "GPR75", "GPR78", "GPR82", "GPR83", "GPR84", "GPR85", "GPR87",
        "GPR88", "GPRC5A", "GPRC5B", "GPRC5C", "GPRC5D", "GPRC6A", "MAS1", "MAS1L", "MRGPRD",
        "MRGPRE", "MRGPRF", "MRGPRG", "MRGPRX1", "MRGPRX2", "MRGPRX3", "MRGPRX4", "TPRA1"
    ],
    "Olfactory receptors": [
        "OR10A2", "OR10A3", "OR10A4", "OR10A5", "OR10A6", "OR10A7", "OR10AC1", "OR10AD1", "OR10AG1", "OR10C1",
        "OR10G2", "OR10G3", "OR10G4", "OR10G6", "OR10G7", "OR10G8", "OR10G9", "OR10H1", "OR10H2", "OR10H3",
        "OR10H4", "OR10H5", "OR10J1", "OR10J3", "OR10J4", "OR10J5", "OR10K1", "OR10K2", "OR10P1", "OR10Q1",
        "OR10R2", "OR10S1", "OR10T2", "OR10V1", "OR10W1", "OR10X1", "OR10Z1", "OR11A1", "OR11G2", "OR11H1",
        "OR11H12", "OR11H2", "OR11H4", "OR11H6", "OR11H7", "OR11L1", "OR12D1", "OR12D2", "OR12D3", "OR13A1",
        "OR13C2", "OR13C3", "OR13C4", "OR13C5", "OR13C7", "OR13C8", "OR13C9", "OR13D1", "OR13F1", "OR13G1",
        "OR13H1", "OR13J1", "OR14A16", "OR14A2", "OR14C36", "OR14I1", "OR14J1", "OR14K1", "OR14L1P", "OR1A1",
        "OR1A2", "OR1B1", "OR1C1", "OR1D2", "OR1D4", "OR1D5", "OR1E1", "OR1E2", "OR1E3", "OR1F1", "OR1G1",
        "OR1I1", "OR1J1", "OR1J2", "OR1J4", "OR1K1", "OR1L1", "OR1L3", "OR1L4", "OR1L6", "OR1L8", "OR1M1",
        "OR1N1", "OR1N2", "OR1P1", "OR1Q1", "OR1S1", "OR1S2", "OR2A12", "OR2A14", "OR2A2", "OR2A25", "OR2A4",
        "OR2A42", "OR2A5", "OR2A7", "OR2AE1", "OR2AG1", "OR2AG2", "OR2AJ1", "OR2AK2", "OR2AP1", "OR2AT4",
        "OR2B11", "OR2B2", "OR2B3", "OR2B6", "OR2C1", "OR2C3", "OR2D2", "OR2D3", "OR2F1", "OR2F2", "OR2G2",
        "OR2G3", "OR2G6", "OR2H1", "OR2H2", "OR2J1", "OR2J2", "OR2J3", "OR2K2", "OR2L13", "OR2L2", "OR2L3",
        "OR2L5", "OR2L8", "OR2M2", "OR2M3", "OR2M4", "OR2M5", "OR2M7", "OR2S2", "OR2T1", "OR2T10", "OR2T11",
        "OR2T12", "OR2T2", "OR2T27", "OR2T29", "OR2T3", "OR2T33", "OR2T34", "OR2T35", "OR2T4", "OR2T5",
        "OR2T6", "OR2T7", "OR2T8", "OR2V1", "OR2V2", "OR2W1", "OR2W3", "OR2Y1", "OR2Z1", "OR3A1", "OR3A2",
        "OR3A3", "OR4A15", "OR4A16", "OR4A47", "OR4A5", "OR4A8", "OR4B1", "OR4C11", "OR4C12", "OR4C13",
        "OR4C15", "OR4C16", "OR4C3", "OR4C45", "OR4C46", "OR4C5", "OR4C6", "OR4D1", "OR4D10", "OR4D11",
        "OR4D2", "OR4D5", "OR4D6", "OR4D9", "OR4E1", "OR4E2", "OR4F15", "OR4F17", "OR4F21", "OR4F29",
        "OR4F4", "OR4F5", "OR4F6", "OR4K1", "OR4K13", "OR4K14", "OR4K15", "OR4K17", "OR4K2", "OR4K3",
        "OR4K5", "OR4L1", "OR4M2", "OR4M2B", "OR4N2", "OR4N4", "OR4N4C", "OR4N5", "OR4P4", "OR4Q2",
        "OR4Q3", "OR4S1", "OR4S2", "OR4X1", "OR4X2", "OR51A2", "OR51A4", "OR51A7", "OR51B2", "OR51B4",
        "OR51B5", "OR51B6", "OR51D1", "OR51E1", "OR51E2", "OR51F1", "OR51F2", "OR51G1", "OR51G2",
        "OR51H1", "OR51I1", "OR51I2", "OR51J1", "OR51L1", "OR51M1", "OR51Q1", "OR51S1", "OR51T1",
        "OR51V1", "OR52A1", "OR52A4P", "OR52A5", "OR52B2", "OR52B4", "OR52B6", "OR52D1", "OR52E1",
        "OR52E2", "OR52E4", "OR52E5", "OR52E6", "OR52E8", "OR52H1", "OR52I1", "OR52I2", "OR52J3",
        "OR52K1", "OR52K2", "OR52L1", "OR52M1", "OR52N1", "OR52N2", "OR52N4", "OR52N5", "OR52P1",
        "OR52R1", "OR52W1", "OR52Z1P", "OR56A1", "OR56A3", "OR56A4", "OR56A5", "OR56B1", "OR56B4",
        "OR5A1", "OR5A2", "OR5AC1", "OR5AC2", "OR5AK2", "OR5AL1", "OR5AN1", "OR5AP2", "OR5AR1",
        "OR5AS1", "OR5AU1", "OR5B12", "OR5B17", "OR5B2", "OR5B21", "OR5B3", "OR5C1", "OR5D13",
        "OR5D14", "OR5D16", "OR5D18", "OR5F1", "OR5G3", "OR5H1", "OR5H14", "OR5H15", "OR5H2",
        "OR5H6", "OR5H8", "OR5I1", "OR5J2", "OR5K1", "OR5K2", "OR5K3", "OR5K4", "OR5L1", "OR5L2",
        "OR5M1", "OR5M10", "OR5M11", "OR5M3", "OR5M8", "OR5M9", "OR5P2", "OR5P3", "OR5T1", "OR5T2",
        "OR5T3", "OR5V1", "OR5W2", "OR6A2", "OR6B1", "OR6B2", "OR6B3", "OR6C1", "OR6C2", "OR6C3",
        "OR6C4", "OR6C6", "OR6C65", "OR6C68", "OR6C70", "OR6C74", "OR6C75", "OR6C76", "OR6F1",
        "OR6J1", "OR6K2", "OR6K3", "OR6K6", "OR6M1", "OR6N1", "OR6N2", "OR6P1", "OR6Q1", "OR6S1",
        "OR6T1", "OR6V1", "OR6X1", "OR6Y1", "OR7A10", "OR7A17", "OR7A5", "OR7C1", "OR7C2", "OR7D2",
        "OR7D4", "OR7E24", "OR7G1", "OR7G2", "OR7G3", "OR8A1", "OR8B12", "OR8B2", "OR8B3", "OR8B4",
        "OR8B8", "OR8D1", "OR8D2", "OR8D4", "OR8G1", "OR8G5", "OR8H1", "OR8H2", "OR8H3", "OR8I2",
        "OR8J1", "OR8J2", "OR8J3", "OR8K1", "OR8K3", "OR8K5", "OR8S1", "OR8U1", "OR8U3", "OR8U8",
        "OR8U9", "OR9A1P", "OR9A2", "OR9A4", "OR9G1", "OR9G4", "OR9G9", "OR9I1", "OR9K2", "OR9Q1",
        "OR9Q2"
    ]
}


@app.route('/api/receptors', methods=['GET'])
def get_receptors():
    try:
        if not os.path.exists(PDF_FOLDER):
            return jsonify({"error": "PDF folder not found"}), 500

        real_files = set([f.replace('.pdf', '') for f in os.listdir(PDF_FOLDER) if f.endswith('.pdf')])

        organized_data = {}

        for family, members in GPCR_FAMILIES.items():
            present_members = []
            for m in members:
                if m in real_files:

                    common_name = ""
                    if m in GPCR_METADATA and 'cn' in GPCR_METADATA[m]:
                        common_name = str(GPCR_METADATA[m]['cn']).strip()

                    short_name = PROTEIN_NAMES.get(m, "")

                    full_search_text = f"{m} {short_name} {common_name}"

                    coupling_info = COUPLING_DATA.get(m, [])

                    present_members.append({
                        "id": m,
                        "text": m,
                        "protein": short_name,
                        "desc": common_name,
                        "search_text": full_search_text,
                        "coupling": coupling_info
                    })

            if present_members:
                organized_data[family] = present_members

        return jsonify(organized_data)

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/api/receptors-complex', methods=['GET'])
def get_receptors_complex():
    try:
        GPCR_METADATA['ACKR1'] = {'acc': 'Q16570', 'cn': ' Atypical chemokine receptor 1 '}
        GPCR_METADATA['ACKR2'] = {'acc': 'O00590', 'cn': ' Atypical chemokine receptor 2 '}
        GPCR_METADATA['ACKR3'] = {'acc': 'P25106', 'cn': ' Atypical chemokine receptor 3 '}
        GPCR_METADATA['ACKR4'] = {'acc': 'Q9NPB9', 'cn': ' Atypical chemokine receptor 4 '}
        GPCR_METADATA['ADGRA1'] = {'acc': 'Q86SQ6', 'cn': ' Adhesion G protein-coupled receptor A1 '}
        GPCR_METADATA['ADGRA2'] = {'acc': 'Q96PE1', 'cn': ' Adhesion G protein-coupled receptor A2 '}
        GPCR_METADATA['ADGRA3'] = {'acc': 'Q8IWK6', 'cn': ' Adhesion G protein-coupled receptor A3 '}
        GPCR_METADATA['ADGRB1'] = {'acc': 'O14514', 'cn': ' Adhesion G protein-coupled receptor B1 '}
        GPCR_METADATA['ADGRB2'] = {'acc': 'O60241', 'cn': ' Adhesion G protein-coupled receptor B2 '}
        GPCR_METADATA['ADGRB3'] = {'acc': 'O60242', 'cn': ' Adhesion G protein-coupled receptor B3 '}
        GPCR_METADATA['ADGRD1'] = {'acc': 'Q6QNK2', 'cn': ' Adhesion G-protein coupled receptor D1 '}
        GPCR_METADATA['ADGRD2'] = {'acc': 'Q7Z7M1', 'cn': ' Adhesion G protein-coupled receptor D2 '}
        GPCR_METADATA['ADGRE1'] = {'acc': 'Q14246', 'cn': ' Adhesion G protein-coupled receptor E1 '}
        GPCR_METADATA['ADGRE2'] = {'acc': 'Q9UHX3', 'cn': ' Adhesion G protein-coupled receptor E2 '}
        GPCR_METADATA['ADGRE3'] = {'acc': 'Q9BY15', 'cn': ' Adhesion G protein-coupled receptor E3 '}
        GPCR_METADATA['ADGRE5'] = {'acc': 'P48960', 'cn': ' Adhesion G protein-coupled receptor E5 '}
        GPCR_METADATA['ADGRF1'] = {'acc': 'Q5T601', 'cn': ' Adhesion G-protein coupled receptor F1 '}
        GPCR_METADATA['ADGRF3'] = {'acc': 'Q8IZF5', 'cn': ' Adhesion G-protein coupled receptor F3 '}
        GPCR_METADATA['ADGRF4'] = {'acc': 'Q8IZF3', 'cn': ' Adhesion G protein-coupled receptor F4 '}
        GPCR_METADATA['ADGRF5'] = {'acc': 'Q8IZF2', 'cn': ' Adhesion G protein-coupled receptor F5 '}
        GPCR_METADATA['ADGRG1'] = {'acc': 'Q9Y653', 'cn': ' Adhesion G-protein coupled receptor G1 '}
        GPCR_METADATA['ADGRG2'] = {'acc': 'Q8IZP9', 'cn': ' Adhesion G-protein coupled receptor G2 '}
        GPCR_METADATA['ADGRG3'] = {'acc': 'Q86Y34', 'cn': ' Adhesion G protein-coupled receptor G3 '}
        GPCR_METADATA['ADGRG4'] = {'acc': 'Q8IZF6', 'cn': ' Adhesion G-protein coupled receptor G4 '}
        GPCR_METADATA['ADGRG5'] = {'acc': 'Q8IZF4', 'cn': ' Adhesion G-protein coupled receptor G5 '}
        GPCR_METADATA['ADGRG6'] = {'acc': 'Q86SQ4', 'cn': ' Adhesion G-protein coupled receptor G6 '}
        GPCR_METADATA['ADGRG7'] = {'acc': 'Q96K78', 'cn': ' Adhesion G-protein coupled receptor G7 '}
        GPCR_METADATA['ADGRL1'] = {'acc': 'O94910', 'cn': ' Adhesion G protein-coupled receptor L1 '}
        GPCR_METADATA['ADGRL2'] = {'acc': 'O95490', 'cn': ' Adhesion G protein-coupled receptor L2 '}
        GPCR_METADATA['ADGRL3'] = {'acc': 'Q9HAR2', 'cn': ' Adhesion G protein-coupled receptor L3 '}
        GPCR_METADATA['ADGRL4'] = {'acc': 'Q9HBW9', 'cn': ' Adhesion G protein-coupled receptor L4 '}
        GPCR_METADATA['ADGRV1'] = {'acc': 'Q8WXG9', 'cn': ' Adhesion G-protein coupled receptor V1 '}
        GPCR_METADATA['AGTR2'] = {'acc': 'P50052', 'cn': ' Type-2 angiotensin II receptor '}
        GPCR_METADATA['BRS3'] = {'acc': 'P32247', 'cn': ' Bombesin receptor subtype-3 '}
        GPCR_METADATA['C3AR1'] = {'acc': 'Q16581', 'cn': ' C3a anaphylatoxin chemotactic receptor '}
        GPCR_METADATA['C5AR1'] = {'acc': 'P21730', 'cn': ' C5a anaphylatoxin chemotactic receptor 1 '}
        GPCR_METADATA['C5AR2'] = {'acc': 'Q9P296', 'cn': ' C5a anaphylatoxin chemotactic receptor 2 '}
        GPCR_METADATA['CALCR'] = {'acc': 'P30988', 'cn': ' Calcitonin receptor '}
        GPCR_METADATA['CALCRL'] = {'acc': 'Q16602', 'cn': ' Calcitonin gene-related peptide type 1 receptor '}
        GPCR_METADATA['CASR'] = {'acc': 'P41180', 'cn': ' Extracellular calcium-sensing receptor '}
        GPCR_METADATA['CCR1'] = {'acc': 'P32246', 'cn': ' C-C chemokine receptor type 1 '}
        GPCR_METADATA['CCR10'] = {'acc': 'P46092', 'cn': ' C-C chemokine receptor type 10 '}
        GPCR_METADATA['CCR2'] = {'acc': 'P41597', 'cn': ' C-C chemokine receptor type 2 '}
        GPCR_METADATA['CCR3'] = {'acc': 'P51677', 'cn': ' C-C chemokine receptor type 3 '}
        GPCR_METADATA['CCR4'] = {'acc': 'P51679', 'cn': ' C-C chemokine receptor type 4 '}
        GPCR_METADATA['CCR5'] = {'acc': 'P51681', 'cn': ' C-C chemokine receptor type 5 '}
        GPCR_METADATA['CCR6'] = {'acc': 'P51684', 'cn': ' C-C chemokine receptor type 6 '}
        GPCR_METADATA['CCR7'] = {'acc': 'P32248', 'cn': ' C-C chemokine receptor type 7 '}
        GPCR_METADATA['CCR8'] = {'acc': 'P51685', 'cn': ' C-C chemokine receptor type 8 '}
        GPCR_METADATA['CCR9'] = {'acc': 'P51686', 'cn': ' C-C chemokine receptor type 9 '}
        GPCR_METADATA['CCRL2'] = {'acc': 'O00421', 'cn': ' C-C chemokine receptor-like 2 '}
        GPCR_METADATA['CELSR1'] = {'acc': 'Q9NYQ6', 'cn': ' Cadherin EGF LAG seven-pass G-type receptor 1 '}
        GPCR_METADATA['CELSR2'] = {'acc': 'Q9HCU4', 'cn': ' Cadherin EGF LAG seven-pass G-type receptor 2 '}
        GPCR_METADATA['CELSR3'] = {'acc': 'Q9NYQ7', 'cn': ' Cadherin EGF LAG seven-pass G-type receptor 3 '}
        GPCR_METADATA['CMKLR1'] = {'acc': 'Q99788', 'cn': ' Chemerin-like receptor 1 '}
        GPCR_METADATA['CMKLR2'] = {'acc': 'P46091', 'cn': ' Chemerin-like receptor 2 '}
        GPCR_METADATA['CRHR1'] = {'acc': 'P34998', 'cn': ' Corticotropin-releasing factor receptor 1 '}
        GPCR_METADATA['CRHR2'] = {'acc': 'Q13324', 'cn': ' Corticotropin-releasing factor receptor 2 '}
        GPCR_METADATA['CX3CR1'] = {'acc': 'P49238', 'cn': ' CX3C chemokine receptor 1 '}
        GPCR_METADATA['CXCR1'] = {'acc': 'P25024', 'cn': ' C-X-C chemokine receptor type 1 '}
        GPCR_METADATA['CXCR2'] = {'acc': 'P25025', 'cn': ' C-X-C chemokine receptor type 2 '}
        GPCR_METADATA['CXCR3'] = {'acc': 'P49682', 'cn': ' C-X-C chemokine receptor type 3 '}
        GPCR_METADATA['CXCR4'] = {'acc': 'P61073', 'cn': ' C-X-C chemokine receptor type 4 '}
        GPCR_METADATA['CXCR5'] = {'acc': 'P32302', 'cn': ' C-X-C chemokine receptor type 5 '}
        GPCR_METADATA['CXCR6'] = {'acc': 'O00574', 'cn': ' C-X-C chemokine receptor type 6 '}
        GPCR_METADATA['FPR3'] = {'acc': 'P25089', 'cn': ' N-formyl peptide receptor 3 '}
        GPCR_METADATA['FSHR'] = {'acc': 'P23945', 'cn': ' Follicle-stimulating hormone receptor '}
        GPCR_METADATA['FZD1'] = {'acc': 'Q9UP38', 'cn': ' Frizzled-1 '}
        GPCR_METADATA['FZD10'] = {'acc': 'Q9ULW2', 'cn': ' Frizzled-10 '}
        GPCR_METADATA['FZD2'] = {'acc': 'Q14332', 'cn': ' Frizzled-2 '}
        GPCR_METADATA['FZD3'] = {'acc': 'Q9NPG1', 'cn': ' Frizzled-3 '}
        GPCR_METADATA['FZD4'] = {'acc': 'Q9ULV1', 'cn': ' Frizzled-4 '}
        GPCR_METADATA['FZD5'] = {'acc': 'Q13467', 'cn': ' Frizzled-5 '}
        GPCR_METADATA['FZD6'] = {'acc': 'O60353', 'cn': ' Frizzled-6 '}
        GPCR_METADATA['FZD7'] = {'acc': 'O75084', 'cn': ' Frizzled-7 '}
        GPCR_METADATA['FZD8'] = {'acc': 'Q9H461', 'cn': ' Frizzled-8 '}
        GPCR_METADATA['FZD9'] = {'acc': 'O00144', 'cn': ' Frizzled-9 '}
        GPCR_METADATA['GABBR1'] = {'acc': 'Q9UBS5', 'cn': ' Gamma-aminobutyric acid type B receptor subunit 1 '}
        GPCR_METADATA['GABBR2'] = {'acc': 'O75899', 'cn': ' Gamma-aminobutyric acid type B receptor subunit 2 '}
        GPCR_METADATA['GCGR'] = {'acc': 'P47871', 'cn': ' Glucagon receptor '}
        GPCR_METADATA['GHRHR'] = {'acc': 'Q02643', 'cn': ' Growth hormone-releasing hormone receptor '}
        GPCR_METADATA['GLP1R'] = {'acc': 'P43220', 'cn': ' Glucagon-like peptide 1 receptor '}
        GPCR_METADATA['GLP2R'] = {'acc': 'O95838', 'cn': ' Glucagon-like peptide 2 receptor '}
        GPCR_METADATA['GPBAR1'] = {'acc': 'Q8TDU6', 'cn': ' G-protein coupled bile acid receptor 1 '}
        GPCR_METADATA['GPER1'] = {'acc': 'Q99527', 'cn': ' G-protein coupled estrogen receptor 1 '}
        GPCR_METADATA['GPR101'] = {'acc': 'Q96P66', 'cn': ' Probable G-protein coupled receptor 101 '}
        GPCR_METADATA['GPR107'] = {'acc': 'Q5VW38', 'cn': ' Protein GPR107 '}
        GPCR_METADATA['GPR12'] = {'acc': 'P47775', 'cn': ' G-protein coupled receptor 12 '}
        GPCR_METADATA['GPR135'] = {'acc': 'Q8IZ08', 'cn': ' G-protein coupled receptor 135 '}
        GPCR_METADATA['GPR137'] = {'acc': 'Q96N19', 'cn': ' Integral membrane protein GPR137 '}
        GPCR_METADATA['GPR139'] = {'acc': 'Q6DWJ6', 'cn': ' Probable G-protein coupled receptor 139 '}
        GPCR_METADATA['GPR141'] = {'acc': 'Q7Z602', 'cn': ' Probable G-protein coupled receptor 141 '}
        GPCR_METADATA['GPR142'] = {'acc': 'Q7Z601', 'cn': ' G protein-coupled receptor 142 '}
        GPCR_METADATA['GPR143'] = {'acc': 'P51810', 'cn': ' G-protein coupled receptor 143 '}
        GPCR_METADATA['GPR146'] = {'acc': 'Q96CH1', 'cn': ' G-protein coupled receptor 146 '}
        GPCR_METADATA['GPR148'] = {'acc': 'Q8TDV2', 'cn': ' Probable G-protein coupled receptor 148 '}
        GPCR_METADATA['GPR149'] = {'acc': 'Q86SP6', 'cn': ' Probable G-protein coupled receptor 149 '}
        GPCR_METADATA['GPR15'] = {'acc': 'P49685', 'cn': ' G-protein coupled receptor 15 '}
        GPCR_METADATA['GPR150'] = {'acc': 'Q8NGU9', 'cn': ' Probable G-protein coupled receptor 150 '}
        GPCR_METADATA['GPR151'] = {'acc': 'Q8TDV0', 'cn': ' G-protein coupled receptor 151 '}
        GPCR_METADATA['GPR152'] = {'acc': 'Q8TDT2', 'cn': ' Probable G-protein coupled receptor 152 '}
        GPCR_METADATA['GPR153'] = {'acc': 'Q6NV75', 'cn': ' Probable G-protein coupled receptor 153 '}
        GPCR_METADATA['GPR156'] = {'acc': 'Q8NFN8', 'cn': ' Probable G-protein coupled receptor 156 '}
        GPCR_METADATA['GPR157'] = {'acc': 'Q5UAW9', 'cn': ' G-protein coupled receptor 157 '}
        GPCR_METADATA['GPR158'] = {'acc': 'Q5T848', 'cn': ' Metabotropic glycine receptor '}
        GPCR_METADATA['GPR160'] = {'acc': 'Q9UJ42', 'cn': ' Probable G-protein coupled receptor 160 '}
        GPCR_METADATA['GPR161'] = {'acc': 'Q8N6U8', 'cn': ' G-protein coupled receptor 161 '}
        GPCR_METADATA['GPR162'] = {'acc': 'Q16538', 'cn': ' Probable G-protein coupled receptor 162 '}
        GPCR_METADATA['GPR171'] = {'acc': 'O14626', 'cn': ' G-protein coupled receptor 171 '}
        GPCR_METADATA['GPR173'] = {'acc': 'Q9NS66', 'cn': ' Probable G-protein coupled receptor 173 '}
        GPCR_METADATA['GPR176'] = {'acc': 'Q14439', 'cn': ' G-protein coupled receptor 176 '}
        GPCR_METADATA['GPR179'] = {'acc': 'Q6PRD1', 'cn': ' Probable G-protein coupled receptor 179 '}
        GPCR_METADATA['GPR18'] = {'acc': 'Q14330', 'cn': ' N-arachidonyl glycine receptor '}
        GPCR_METADATA['GPR182'] = {'acc': 'O15218', 'cn': ' Atypical chemokine receptor 5 '}
        GPCR_METADATA['GPR19'] = {'acc': 'Q15760', 'cn': ' Probable G-protein coupled receptor 19 '}
        GPCR_METADATA['GPR20'] = {'acc': 'Q99678', 'cn': ' G-protein coupled receptor 20 '}
        GPCR_METADATA['GPR21'] = {'acc': 'Q99679', 'cn': ' Probable G-protein coupled receptor 21 '}
        GPCR_METADATA['GPR22'] = {'acc': 'Q99680', 'cn': ' G-protein coupled receptor 22 '}
        GPCR_METADATA['GPR25'] = {'acc': 'O00155', 'cn': ' C-X-C chemokine receptor GPR25 '}
        GPCR_METADATA['GPR26'] = {'acc': 'Q8NDV2', 'cn': ' G-protein coupled receptor 26 '}
        GPCR_METADATA['GPR27'] = {'acc': 'Q9NS67', 'cn': ' Probable G-protein coupled receptor 27 '}
        GPCR_METADATA['GPR3'] = {'acc': 'P46089', 'cn': ' G-protein coupled receptor 3 '}
        GPCR_METADATA['GPR31'] = {'acc': 'O00270', 'cn': ' 12-(S)-hydroxy-5,8,10,14-eicosatetraenoic acid receptor '}
        GPCR_METADATA['GPR32'] = {'acc': 'O75388', 'cn': ' Probable G-protein coupled receptor 32 '}
        GPCR_METADATA['GPR33'] = {'acc': 'Q49SQ1', 'cn': ' Probable G-protein coupled receptor 33 '}
        GPCR_METADATA['GPR37'] = {'acc': 'O15354', 'cn': ' Prosaposin receptor GPR37 '}
        GPCR_METADATA['GPR37L1'] = {'acc': 'O60883', 'cn': ' G-protein coupled receptor 37-like 1 '}
        GPCR_METADATA['GPR39'] = {'acc': 'O43194', 'cn': ' G-protein coupled receptor 39 '}
        GPCR_METADATA['GPR4'] = {'acc': 'P46093', 'cn': ' G-protein coupled receptor 4 '}
        GPCR_METADATA['GPR42'] = {'acc': 'O15529', 'cn': ' G-protein coupled receptor 42 '}
        GPCR_METADATA['GPR45'] = {'acc': 'Q9Y5Y3', 'cn': ' Probable G-protein coupled receptor 45 '}
        GPCR_METADATA['GPR50'] = {'acc': 'Q13585', 'cn': ' Melatonin-related receptor '}
        GPCR_METADATA['GPR52'] = {'acc': 'Q9Y2T5', 'cn': ' G-protein coupled receptor 52 '}
        GPCR_METADATA['GPR6'] = {'acc': 'P46095', 'cn': ' G-protein coupled receptor 6 '}
        GPCR_METADATA['GPR61'] = {'acc': 'Q9BZJ8', 'cn': ' G-protein coupled receptor 61 '}
        GPCR_METADATA['GPR62'] = {'acc': 'Q9BZJ7', 'cn': ' G-protein coupled receptor 62 '}
        GPCR_METADATA['GPR63'] = {'acc': 'Q9BZJ6', 'cn': ' Probable G-protein coupled receptor 63 '}
        GPCR_METADATA['GPR65'] = {'acc': 'Q8IYL9', 'cn': ' G-protein coupled receptor 65 '}
        GPCR_METADATA['GPR68'] = {'acc': 'Q15743', 'cn': ' G-protein coupled receptor 68 '}
        GPCR_METADATA['GPR75'] = {'acc': 'O95800', 'cn': ' Probable G-protein coupled receptor 75 '}
        GPCR_METADATA['GPR78'] = {'acc': 'Q96P69', 'cn': ' G-protein coupled receptor 78 '}
        GPCR_METADATA['GPR82'] = {'acc': 'Q96P67', 'cn': ' Probable G-protein coupled receptor 82 '}
        GPCR_METADATA['GPR83'] = {'acc': 'Q9NYM4', 'cn': ' G-protein coupled receptor 83 '}
        GPCR_METADATA['GPR85'] = {'acc': 'P60893', 'cn': ' Probable G-protein coupled receptor 85 '}
        GPCR_METADATA['GPR87'] = {'acc': 'Q9BY21', 'cn': ' G-protein coupled receptor 87 '}
        GPCR_METADATA['GPR88'] = {'acc': 'Q9GZN0', 'cn': ' G protein-coupled receptor 88 '}
        GPCR_METADATA['GPRC5A'] = {'acc': 'Q8NFJ5', 'cn': ' Retinoic acid-induced protein 3 '}
        GPCR_METADATA['GPRC5B'] = {'acc': 'Q9NZH0', 'cn': ' G-protein coupled receptor family C group 5 member B '}
        GPCR_METADATA['GPRC5C'] = {'acc': 'Q9NQ84', 'cn': ' G-protein coupled receptor family C group 5 member C '}
        GPCR_METADATA['GPRC5D'] = {'acc': 'Q9NZD1', 'cn': ' G-protein coupled receptor family C group 5 member D '}
        GPCR_METADATA['GPRC6A'] = {'acc': 'Q5T6X5', 'cn': ' G-protein coupled receptor family C group 6 member A '}
        GPCR_METADATA['GRM1'] = {'acc': 'Q13255', 'cn': ' Metabotropic glutamate receptor 1 '}
        GPCR_METADATA['GRM2'] = {'acc': 'Q14416', 'cn': ' Metabotropic glutamate receptor 2 '}
        GPCR_METADATA['GRM3'] = {'acc': 'Q14832', 'cn': ' Metabotropic glutamate receptor 3 '}
        GPCR_METADATA['GRM4'] = {'acc': 'Q14833', 'cn': ' Metabotropic glutamate receptor 4 '}
        GPCR_METADATA['GRM5'] = {'acc': 'P41594', 'cn': ' Metabotropic glutamate receptor 5 '}
        GPCR_METADATA['GRM6'] = {'acc': 'O15303', 'cn': ' Metabotropic glutamate receptor 6 '}
        GPCR_METADATA['GRM7'] = {'acc': 'Q14831', 'cn': ' Metabotropic glutamate receptor 7 '}
        GPCR_METADATA['GRM8'] = {'acc': 'O00222', 'cn': ' Metabotropic glutamate receptor 8 '}
        GPCR_METADATA['HCAR1'] = {'acc': 'Q9BXC0', 'cn': ' Hydroxycarboxylic acid receptor 1 '}
        GPCR_METADATA['HCAR2'] = {'acc': 'Q8TDS4', 'cn': ' Hydroxycarboxylic acid receptor 2 '}
        GPCR_METADATA['HCAR3'] = {'acc': 'P49019', 'cn': ' Hydroxycarboxylic acid receptor 3 '}
        GPCR_METADATA['HTR5A'] = {'acc': 'P47898', 'cn': ' 5-hydroxytryptamine receptor 5A '}
        GPCR_METADATA['LGR4'] = {'acc': 'Q9BXB1', 'cn': ' Leucine-rich repeat-containing G-protein coupled receptor 4 '}
        GPCR_METADATA['LGR5'] = {'acc': 'O75473', 'cn': ' Leucine-rich repeat-containing G-protein coupled receptor 5 '}
        GPCR_METADATA['LGR6'] = {'acc': 'Q9HBX8', 'cn': ' Leucine-rich repeat-containing G-protein coupled receptor 6 '}
        GPCR_METADATA['LHCGR'] = {'acc': 'P22888', 'cn': ' Lutropin-choriogonadotropic hormone receptor '}
        GPCR_METADATA['MAS1'] = {'acc': 'P04201', 'cn': ' Proto-oncogene Mas '}
        GPCR_METADATA['MAS1L'] = {'acc': 'P35410', 'cn': ' Mas-related G-protein coupled receptor MRG '}
        GPCR_METADATA['MC2R'] = {'acc': 'Q01718', 'cn': ' Adrenocorticotropic hormone receptor '}
        GPCR_METADATA['MRGPRD'] = {'acc': 'Q8TDS7', 'cn': ' Mas-related G-protein coupled receptor member D '}
        GPCR_METADATA['MRGPRE'] = {'acc': 'Q86SM8', 'cn': ' Mas-related G-protein coupled receptor member E '}
        GPCR_METADATA['MRGPRF'] = {'acc': 'Q96AM1', 'cn': ' Mas-related G-protein coupled receptor member F '}
        GPCR_METADATA['MRGPRG'] = {'acc': 'Q86SM5', 'cn': ' Mas-related G-protein coupled receptor member G '}
        GPCR_METADATA['MRGPRX3'] = {'acc': 'Q96LB0', 'cn': ' Mas-related G-protein coupled receptor member X3 '}
        GPCR_METADATA['MRGPRX4'] = {'acc': 'Q96LA9', 'cn': ' Mas-related G-protein coupled receptor member X4 '}
        GPCR_METADATA['NPBWR2'] = {'acc': 'P48146', 'cn': ' Neuropeptides B/W receptor type 2 '}
        GPCR_METADATA['NPSR1'] = {'acc': 'Q6W5P4', 'cn': ' Neuropeptide S receptor '}
        GPCR_METADATA['NPY1R'] = {'acc': 'P25929', 'cn': ' Neuropeptide Y receptor type 1 '}
        GPCR_METADATA['NPY2R'] = {'acc': 'P49146', 'cn': ' Neuropeptide Y receptor type 2 '}
        GPCR_METADATA['NPY4R'] = {'acc': 'P50391', 'cn': ' Neuropeptide Y receptor type 4 '}
        GPCR_METADATA['NPY5R'] = {'acc': 'Q15761', 'cn': ' Neuropeptide Y receptor type 5 '}
        GPCR_METADATA['NTSR1'] = {'acc': 'P30989', 'cn': ' Neurotensin receptor type 1 '}
        GPCR_METADATA['NTSR2'] = {'acc': 'O95665', 'cn': ' Neurotensin receptor type 2 '}
        GPCR_METADATA['OPN1LW'] = {'acc': 'P04000', 'cn': ' Long-wave-sensitive opsin 1 '}
        GPCR_METADATA['OPN1MW'] = {'acc': 'P04001', 'cn': ' Medium-wave-sensitive opsin 1 '}
        GPCR_METADATA['OPN1SW'] = {'acc': 'P03999', 'cn': ' Short-wave-sensitive opsin 1 '}
        GPCR_METADATA['OPN3'] = {'acc': 'Q9H1Y3', 'cn': ' Opsin-3 '}
        GPCR_METADATA['OPN4'] = {'acc': 'Q9UHM6', 'cn': ' Melanopsin '}
        GPCR_METADATA['OPN5'] = {'acc': 'Q6U736', 'cn': ' Opsin-5 '}
        GPCR_METADATA['OR10A2'] = {'acc': 'Q9H208', 'cn': ' Olfactory receptor 10A2 '}
        GPCR_METADATA['OR10A3'] = {'acc': 'P58181', 'cn': ' Olfactory receptor 10A3 '}
        GPCR_METADATA['OR10A4'] = {'acc': 'Q9H209', 'cn': ' Olfactory receptor 10A4 '}
        GPCR_METADATA['OR10A5'] = {'acc': 'Q9H207', 'cn': ' Olfactory receptor 10A5 '}
        GPCR_METADATA['OR10A6'] = {'acc': 'Q8NH74', 'cn': ' Olfactory receptor 10A6 '}
        GPCR_METADATA['OR10A7'] = {'acc': 'Q8NGE5', 'cn': ' Olfactory receptor 10A7 '}
        GPCR_METADATA['OR10AC1'] = {'acc': 'Q8NH08', 'cn': ' Olfactory receptor 10AC1 '}
        GPCR_METADATA['OR10AD1'] = {'acc': 'Q8NGE0', 'cn': ' Olfactory receptor 10AD1 '}
        GPCR_METADATA['OR10AG1'] = {'acc': 'Q8NH19', 'cn': ' Olfactory receptor 10AG1 '}
        GPCR_METADATA['OR10C1'] = {'acc': 'Q96KK4', 'cn': ' Olfactory receptor 10C1 '}
        GPCR_METADATA['OR10G2'] = {'acc': 'Q8NGC3', 'cn': ' Olfactory receptor 10G2 '}
        GPCR_METADATA['OR10G3'] = {'acc': 'Q8NGC4', 'cn': ' Olfactory receptor 10G3 '}
        GPCR_METADATA['OR10G4'] = {'acc': 'Q8NGN3', 'cn': ' Olfactory receptor 10G4 '}
        GPCR_METADATA['OR10G6'] = {'acc': 'Q8NH81', 'cn': ' Olfactory receptor 10G6 '}
        GPCR_METADATA['OR10G7'] = {'acc': 'Q8NGN6', 'cn': ' Olfactory receptor 10G7 '}
        GPCR_METADATA['OR10G8'] = {'acc': 'Q8NGN5', 'cn': ' Olfactory receptor 10G8 '}
        GPCR_METADATA['OR10G9'] = {'acc': 'Q8NGN4', 'cn': ' Olfactory receptor 10G9 '}
        GPCR_METADATA['OR10H1'] = {'acc': 'Q9Y4A9', 'cn': ' Olfactory receptor 10H1 '}
        GPCR_METADATA['OR10H2'] = {'acc': 'O60403', 'cn': ' Olfactory receptor 10H2 '}
        GPCR_METADATA['OR10H3'] = {'acc': 'O60404', 'cn': ' Olfactory receptor 10H3 '}
        GPCR_METADATA['OR10H4'] = {'acc': 'Q8NGA5', 'cn': ' Olfactory receptor 10H4 '}
        GPCR_METADATA['OR10H5'] = {'acc': 'Q8NGA6', 'cn': ' Olfactory receptor 10H5 '}
        GPCR_METADATA['OR10J1'] = {'acc': 'P30954', 'cn': ' Olfactory receptor 10J1 '}
        GPCR_METADATA['OR10J3'] = {'acc': 'Q5JRS4', 'cn': ' Olfactory receptor 10J3 '}
        GPCR_METADATA['OR10J4'] = {'acc': 'P0C629', 'cn': ' Olfactory receptor 10J4 '}
        GPCR_METADATA['OR10J5'] = {'acc': 'Q8NHC4', 'cn': ' Olfactory receptor 10J5 '}
        GPCR_METADATA['OR10K1'] = {'acc': 'Q8NGX5', 'cn': ' Olfactory receptor 10K1 '}
        GPCR_METADATA['OR10K2'] = {'acc': 'Q6IF99', 'cn': ' Olfactory receptor 10K2 '}
        GPCR_METADATA['OR10P1'] = {'acc': 'Q8NGE3', 'cn': ' Olfactory receptor 10P1 '}
        GPCR_METADATA['OR10Q1'] = {'acc': 'Q8NGQ4', 'cn': ' Olfactory receptor 10Q1 '}
        GPCR_METADATA['OR10R2'] = {'acc': 'Q8NGX6', 'cn': ' Olfactory receptor 10R2 '}
        GPCR_METADATA['OR10S1'] = {'acc': 'Q8NGN2', 'cn': ' Olfactory receptor 10S1 '}
        GPCR_METADATA['OR10T2'] = {'acc': 'Q8NGX3', 'cn': ' Olfactory receptor 10T2 '}
        GPCR_METADATA['OR10V1'] = {'acc': 'Q8NGI7', 'cn': ' Olfactory receptor 10V1 '}
        GPCR_METADATA['OR10W1'] = {'acc': 'Q8NGF6', 'cn': ' Olfactory receptor 10W1 '}
        GPCR_METADATA['OR10X1'] = {'acc': 'Q8NGY0', 'cn': ' Olfactory receptor 10X1 '}
        GPCR_METADATA['OR10Z1'] = {'acc': 'Q8NGY1', 'cn': ' Olfactory receptor 10Z1 '}
        GPCR_METADATA['OR11A1'] = {'acc': 'Q9GZK7', 'cn': ' Olfactory receptor 11A1 '}
        GPCR_METADATA['OR11G2'] = {'acc': 'Q8NGC1', 'cn': ' Olfactory receptor 11G2 '}
        GPCR_METADATA['OR11H1'] = {'acc': 'Q8NG94', 'cn': ' Olfactory receptor 11H1 '}
        GPCR_METADATA['OR11H12'] = {'acc': 'B2RN74', 'cn': ' Olfactory receptor 11H12 '}
        GPCR_METADATA['OR11H2'] = {'acc': 'Q8NH07', 'cn': ' Olfactory receptor 11H2 '}
        GPCR_METADATA['OR11H4'] = {'acc': 'Q8NGC9', 'cn': ' Olfactory receptor 11H4 '}
        GPCR_METADATA['OR11H6'] = {'acc': 'Q8NGC7', 'cn': ' Olfactory receptor 11H6 '}
        GPCR_METADATA['OR11H7'] = {'acc': 'Q8NGC8', 'cn': ' Olfactory receptor 11H7 '}
        GPCR_METADATA['OR11L1'] = {'acc': 'Q8NGX0', 'cn': ' Olfactory receptor 11L1 '}
        GPCR_METADATA['OR12D1'] = {'acc': 'P0DN82', 'cn': ' Olfactory receptor 12D1 '}
        GPCR_METADATA['OR12D2'] = {'acc': 'P58182', 'cn': ' Olfactory receptor 12D2 '}
        GPCR_METADATA['OR12D3'] = {'acc': 'Q9UGF7', 'cn': ' Olfactory receptor 12D3 '}
        GPCR_METADATA['OR13A1'] = {'acc': 'Q8NGR1', 'cn': ' Olfactory receptor 13A1 '}
        GPCR_METADATA['OR13C2'] = {'acc': 'Q8NGS9', 'cn': ' Olfactory receptor 13C2 '}
        GPCR_METADATA['OR13C3'] = {'acc': 'Q8NGS6', 'cn': ' Olfactory receptor 13C3 '}
        GPCR_METADATA['OR13C4'] = {'acc': 'Q8NGS5', 'cn': ' Olfactory receptor 13C4 '}
        GPCR_METADATA['OR13C5'] = {'acc': 'Q8NGS8', 'cn': ' Olfactory receptor 13C5 '}
        GPCR_METADATA['OR13C7'] = {'acc': 'P0DN81', 'cn': ' Olfactory receptor 13C7 '}
        GPCR_METADATA['OR13C8'] = {'acc': 'Q8NGS7', 'cn': ' Olfactory receptor 13C8 '}
        GPCR_METADATA['OR13C9'] = {'acc': 'Q8NGT0', 'cn': ' Olfactory receptor 13C9 '}
        GPCR_METADATA['OR13D1'] = {'acc': 'Q8NGV5', 'cn': ' Olfactory receptor 13D1 '}
        GPCR_METADATA['OR13F1'] = {'acc': 'Q8NGS4', 'cn': ' Olfactory receptor 13F1 '}
        GPCR_METADATA['OR13G1'] = {'acc': 'Q8NGZ3', 'cn': ' Olfactory receptor 13G1 '}
        GPCR_METADATA['OR13H1'] = {'acc': 'Q8NG92', 'cn': ' Olfactory receptor 13H1 '}
        GPCR_METADATA['OR13J1'] = {'acc': 'Q8NGT2', 'cn': ' Olfactory receptor 13J1 '}
        GPCR_METADATA['OR14A16'] = {'acc': 'Q8NHC5', 'cn': ' Olfactory receptor 14A16 '}
        GPCR_METADATA['OR14A2'] = {'acc': 'Q96R54', 'cn': ' Olfactory receptor 14A2 '}
        GPCR_METADATA['OR14C36'] = {'acc': 'Q8NHC7', 'cn': ' Olfactory receptor 14C36 '}
        GPCR_METADATA['OR14I1'] = {'acc': 'A6ND48', 'cn': ' Olfactory receptor 14I1 '}
        GPCR_METADATA['OR14J1'] = {'acc': 'Q9UGF5', 'cn': ' Olfactory receptor 14J1 '}
        GPCR_METADATA['OR14K1'] = {'acc': 'Q8NGZ2', 'cn': ' Olfactory receptor 14K1 '}
        GPCR_METADATA['OR14L1P'] = {'acc': 'Q8NHC6', 'cn': ' Olfactory receptor 14L1 '}
        GPCR_METADATA['OR1A1'] = {'acc': 'Q9P1Q5', 'cn': ' Olfactory receptor 1A1 '}
        GPCR_METADATA['OR1A2'] = {'acc': 'Q9Y585', 'cn': ' Olfactory receptor 1A2 '}
        GPCR_METADATA['OR1B1'] = {'acc': 'Q8NGR6', 'cn': ' Olfactory receptor 1B1 '}
        GPCR_METADATA['OR1C1'] = {'acc': 'Q15619', 'cn': ' Olfactory receptor 1C1 '}
        GPCR_METADATA['OR1D2'] = {'acc': 'P34982', 'cn': ' Olfactory receptor 1D2 '}
        GPCR_METADATA['OR1D4'] = {'acc': 'P47884', 'cn': ' Olfactory receptor 1D4 '}
        GPCR_METADATA['OR1D5'] = {'acc': 'P58170', 'cn': ' Olfactory receptor 1D5 '}
        GPCR_METADATA['OR1E1'] = {'acc': 'P30953', 'cn': ' Olfactory receptor 1E1 '}
        GPCR_METADATA['OR1E2'] = {'acc': 'P47887', 'cn': ' Olfactory receptor 1E2 '}
        GPCR_METADATA['OR1E3'] = {'acc': 'Q8WZA6', 'cn': ' Olfactory receptor 1E3 '}
        GPCR_METADATA['OR1F1'] = {'acc': 'O43749', 'cn': ' Olfactory receptor 1F1 '}
        GPCR_METADATA['OR1G1'] = {'acc': 'P47890', 'cn': ' Olfactory receptor 1G1 '}
        GPCR_METADATA['OR1I1'] = {'acc': 'O60431', 'cn': ' Olfactory receptor 1I1 '}
        GPCR_METADATA['OR1J1'] = {'acc': 'Q8NGS3', 'cn': ' Olfactory receptor 1J1 '}
        GPCR_METADATA['OR1J2'] = {'acc': 'Q8NGS2', 'cn': ' Olfactory receptor 1J2 '}
        GPCR_METADATA['OR1J4'] = {'acc': 'Q8NGS1', 'cn': ' Olfactory receptor 1J4 '}
        GPCR_METADATA['OR1K1'] = {'acc': 'Q8NGR3', 'cn': ' Olfactory receptor 1K1 '}
        GPCR_METADATA['OR1L1'] = {'acc': 'Q8NH94', 'cn': ' Olfactory receptor 1L1 '}
        GPCR_METADATA['OR1L3'] = {'acc': 'Q8NH93', 'cn': ' Olfactory receptor 1L3 '}
        GPCR_METADATA['OR1L4'] = {'acc': 'Q8NGR5', 'cn': ' Olfactory receptor 1L4 '}
        GPCR_METADATA['OR1L6'] = {'acc': 'Q8NGR2', 'cn': ' Olfactory receptor 1L6 '}
        GPCR_METADATA['OR1L8'] = {'acc': 'Q8NGR8', 'cn': ' Olfactory receptor 1L8 '}
        GPCR_METADATA['OR1M1'] = {'acc': 'Q8NGA1', 'cn': ' Olfactory receptor 1M1 '}
        GPCR_METADATA['OR1N1'] = {'acc': 'Q8NGS0', 'cn': ' Olfactory receptor 1N1 '}
        GPCR_METADATA['OR1N2'] = {'acc': 'Q8NGR9', 'cn': ' Olfactory receptor 1N2 '}
        GPCR_METADATA['OR1P1'] = {'acc': 'Q8NH06', 'cn': ' Olfactory receptor 1P1 '}
        GPCR_METADATA['OR1Q1'] = {'acc': 'Q15612', 'cn': ' Olfactory receptor 1Q1 '}
        GPCR_METADATA['OR1S1'] = {'acc': 'Q8NH92', 'cn': ' Olfactory receptor 1S1 '}
        GPCR_METADATA['OR1S2'] = {'acc': 'Q8NGQ3', 'cn': ' Olfactory receptor 1S2 '}
        GPCR_METADATA['OR2A12'] = {'acc': 'Q8NGT7', 'cn': ' Olfactory receptor 2A12 '}
        GPCR_METADATA['OR2A14'] = {'acc': 'Q96R47', 'cn': ' Olfactory receptor 2A14 '}
        GPCR_METADATA['OR2A2'] = {'acc': 'Q6IF42', 'cn': ' Olfactory receptor 2A2 '}
        GPCR_METADATA['OR2A25'] = {'acc': 'A4D2G3', 'cn': ' Olfactory receptor 2A25 '}
        GPCR_METADATA['OR2A4'] = {'acc': 'O95047', 'cn': ' Olfactory receptor 2A4 '}
        GPCR_METADATA['OR2A42'] = {'acc': 'Q8NGT9', 'cn': ' Olfactory receptor 2A1/2A42 '}
        GPCR_METADATA['OR2A5'] = {'acc': 'Q96R48', 'cn': ' Olfactory receptor 2A5 '}
        GPCR_METADATA['OR2A7'] = {'acc': 'Q96R45', 'cn': ' Olfactory receptor 2A7 '}
        GPCR_METADATA['OR2AE1'] = {'acc': 'Q8NHA4', 'cn': ' Olfactory receptor 2AE1 '}
        GPCR_METADATA['OR2AG1'] = {'acc': 'Q9H205', 'cn': ' Olfactory receptor 2AG1 '}
        GPCR_METADATA['OR2AG2'] = {'acc': 'A6NM03', 'cn': ' Olfactory receptor 2AG2 '}
        GPCR_METADATA['OR2AJ1'] = {'acc': 'Q8NGZ0', 'cn': ' Olfactory receptor 2AJ1 '}
        GPCR_METADATA['OR2AK2'] = {'acc': 'Q8NG84', 'cn': ' Olfactory receptor 2AK2 '}
        GPCR_METADATA['OR2AP1'] = {'acc': 'Q8NGE2', 'cn': ' Olfactory receptor 2AP1 '}
        GPCR_METADATA['OR2AT4'] = {'acc': 'A6NND4', 'cn': ' Olfactory receptor 2AT4 '}
        GPCR_METADATA['OR2B11'] = {'acc': 'Q5JQS5', 'cn': ' Olfactory receptor 2B11 '}
        GPCR_METADATA['OR2B2'] = {'acc': 'Q9GZK3', 'cn': ' Olfactory receptor 2B2 '}
        GPCR_METADATA['OR2B3'] = {'acc': 'O76000', 'cn': ' Putative olfactory receptor 2B3 '}
        GPCR_METADATA['OR2B6'] = {'acc': 'P58173', 'cn': ' Olfactory receptor 2B6 '}
        GPCR_METADATA['OR2C1'] = {'acc': 'O95371', 'cn': ' Olfactory receptor 2C1 '}
        GPCR_METADATA['OR2C3'] = {'acc': 'Q8N628', 'cn': ' Olfactory receptor 2C3 '}
        GPCR_METADATA['OR2D2'] = {'acc': 'Q9H210', 'cn': ' Olfactory receptor 2D2 '}
        GPCR_METADATA['OR2D3'] = {'acc': 'Q8NGH3', 'cn': ' Olfactory receptor 2D3 '}
        GPCR_METADATA['OR2F1'] = {'acc': 'Q13607', 'cn': ' Olfactory receptor 2F1 '}
        GPCR_METADATA['OR2F2'] = {'acc': 'O95006', 'cn': ' Olfactory receptor 2F2 '}
        GPCR_METADATA['OR2G2'] = {'acc': 'Q8NGZ5', 'cn': ' Olfactory receptor 2G2 '}
        GPCR_METADATA['OR2G3'] = {'acc': 'Q8NGZ4', 'cn': ' Olfactory receptor 2G3 '}
        GPCR_METADATA['OR2G6'] = {'acc': 'Q5TZ20', 'cn': ' Olfactory receptor 2G6 '}
        GPCR_METADATA['OR2H1'] = {'acc': 'Q9GZK4', 'cn': ' Olfactory receptor 2H1 '}
        GPCR_METADATA['OR2H2'] = {'acc': 'O95918', 'cn': ' Olfactory receptor 2H2 '}
        GPCR_METADATA['OR2J1'] = {'acc': 'Q9GZK6', 'cn': ' Olfactory receptor 2J1 '}
        GPCR_METADATA['OR2J2'] = {'acc': 'O76002', 'cn': ' Olfactory receptor 2J2 '}
        GPCR_METADATA['OR2J3'] = {'acc': 'O76001', 'cn': ' Olfactory receptor 2J3 '}
        GPCR_METADATA['OR2K2'] = {'acc': 'Q8NGT1', 'cn': ' Olfactory receptor 2K2 '}
        GPCR_METADATA['OR2L13'] = {'acc': 'Q8N349', 'cn': ' Olfactory receptor 2L13 '}
        GPCR_METADATA['OR2L2'] = {'acc': 'Q8NH16', 'cn': ' Olfactory receptor 2L2 '}
        GPCR_METADATA['OR2L3'] = {'acc': 'Q8NG85', 'cn': ' Olfactory receptor 2L3 '}
        GPCR_METADATA['OR2L5'] = {'acc': 'Q8NG80', 'cn': ' Olfactory receptor 2L5 '}
        GPCR_METADATA['OR2L8'] = {'acc': 'Q8NGY9', 'cn': ' Olfactory receptor 2L8 '}
        GPCR_METADATA['OR2M2'] = {'acc': 'Q96R28', 'cn': ' Olfactory receptor 2M2 '}
        GPCR_METADATA['OR2M3'] = {'acc': 'Q8NG83', 'cn': ' Olfactory receptor 2M3 '}
        GPCR_METADATA['OR2M4'] = {'acc': 'Q96R27', 'cn': ' Olfactory receptor 2M4 '}
        GPCR_METADATA['OR2M5'] = {'acc': 'A3KFT3', 'cn': ' Olfactory receptor 2M5 '}
        GPCR_METADATA['OR2M7'] = {'acc': 'Q8NG81', 'cn': ' Olfactory receptor 2M7 '}
        GPCR_METADATA['OR2S2'] = {'acc': 'Q9NQN1', 'cn': ' Olfactory receptor 2S2 '}
        GPCR_METADATA['OR2T1'] = {'acc': 'O43869', 'cn': ' Olfactory receptor 2T1 '}
        GPCR_METADATA['OR2T10'] = {'acc': 'Q8NGZ9', 'cn': ' Olfactory receptor 2T10 '}
        GPCR_METADATA['OR2T11'] = {'acc': 'Q8NH01', 'cn': ' Olfactory receptor 2T11 '}
        GPCR_METADATA['OR2T12'] = {'acc': 'Q8NG77', 'cn': ' Olfactory receptor 2T12 '}
        GPCR_METADATA['OR2T2'] = {'acc': 'Q6IF00', 'cn': ' Olfactory receptor 2T2 '}
        GPCR_METADATA['OR2T27'] = {'acc': 'Q8NH04', 'cn': ' Olfactory receptor 2T27 '}
        GPCR_METADATA['OR2T29'] = {'acc': 'Q8NH02', 'cn': ' Olfactory receptor 2T29 '}
        GPCR_METADATA['OR2T3'] = {'acc': 'Q8NH03', 'cn': ' Olfactory receptor 2T3 '}
        GPCR_METADATA['OR2T33'] = {'acc': 'Q8NG76', 'cn': ' Olfactory receptor 2T33 '}
        GPCR_METADATA['OR2T34'] = {'acc': 'Q8NGX1', 'cn': ' Olfactory receptor 2T34 '}
        GPCR_METADATA['OR2T35'] = {'acc': 'Q8NGX2', 'cn': ' Olfactory receptor 2T35 '}
        GPCR_METADATA['OR2T4'] = {'acc': 'Q8NH00', 'cn': ' Olfactory receptor 2T4 '}
        GPCR_METADATA['OR2T5'] = {'acc': 'Q6IEZ7', 'cn': ' Olfactory receptor 2T5 '}
        GPCR_METADATA['OR2T6'] = {'acc': 'Q8NHC8', 'cn': ' Olfactory receptor 2T6 '}
        GPCR_METADATA['OR2T7'] = {'acc': 'P0C7T2', 'cn': ' Olfactory receptor 2T7 '}
        GPCR_METADATA['OR2T8'] = {'acc': 'A6NH00', 'cn': ' Olfactory receptor 2T8 '}
        GPCR_METADATA['OR2V1'] = {'acc': 'Q8NHB1', 'cn': ' Olfactory receptor 2V1 '}
        GPCR_METADATA['OR2V2'] = {'acc': 'Q96R30', 'cn': ' Olfactory receptor 2V2 '}
        GPCR_METADATA['OR2W1'] = {'acc': 'Q9Y3N9', 'cn': ' Olfactory receptor 2W1 '}
        GPCR_METADATA['OR2W3'] = {'acc': 'Q7Z3T1', 'cn': ' Olfactory receptor 2W3 '}
        GPCR_METADATA['OR2Y1'] = {'acc': 'Q8NGV0', 'cn': ' Olfactory receptor 2Y1 '}
        GPCR_METADATA['OR2Z1'] = {'acc': 'Q8NG97', 'cn': ' Olfactory receptor 2Z1 '}
        GPCR_METADATA['OR3A1'] = {'acc': 'P47881', 'cn': ' Olfactory receptor 3A1 '}
        GPCR_METADATA['OR3A2'] = {'acc': 'P47893', 'cn': ' Olfactory receptor 3A2 '}
        GPCR_METADATA['OR3A3'] = {'acc': 'P47888', 'cn': ' Olfactory receptor 3A3 '}
        GPCR_METADATA['OR4A15'] = {'acc': 'Q8NGL6', 'cn': ' Olfactory receptor 4A15 '}
        GPCR_METADATA['OR4A16'] = {'acc': 'Q8NH70', 'cn': ' Olfactory receptor 4A16 '}
        GPCR_METADATA['OR4A47'] = {'acc': 'Q6IF82', 'cn': ' Olfactory receptor 4A47 '}
        GPCR_METADATA['OR4A5'] = {'acc': 'Q8NH83', 'cn': ' Olfactory receptor 4A5 '}
        GPCR_METADATA['OR4A8'] = {'acc': 'P0C604', 'cn': ' Olfactory receptor 4A8 '}
        GPCR_METADATA['OR4B1'] = {'acc': 'Q8NGF8', 'cn': ' Olfactory receptor 4B1 '}
        GPCR_METADATA['OR4C11'] = {'acc': 'Q6IEV9', 'cn': ' Olfactory receptor 4C11 '}
        GPCR_METADATA['OR4C12'] = {'acc': 'Q96R67', 'cn': ' Olfactory receptor 4C12 '}
        GPCR_METADATA['OR4C13'] = {'acc': 'Q8NGP0', 'cn': ' Olfactory receptor 4C13 '}
        GPCR_METADATA['OR4C15'] = {'acc': 'Q8NGM1', 'cn': ' Olfactory receptor 4C15 '}
        GPCR_METADATA['OR4C16'] = {'acc': 'Q8NGL9', 'cn': ' Olfactory receptor 4C16 '}
        GPCR_METADATA['OR4C3'] = {'acc': 'Q8NH37', 'cn': ' Olfactory receptor 4C3 '}
        GPCR_METADATA['OR4C45'] = {'acc': 'A6NMZ5', 'cn': ' Olfactory receptor 4C45 '}
        GPCR_METADATA['OR4C46'] = {'acc': 'A6NHA9', 'cn': ' Olfactory receptor 4C46 '}
        GPCR_METADATA['OR4C5'] = {'acc': 'Q8NGB2', 'cn': ' Olfactory receptor 4C5 '}
        GPCR_METADATA['OR4C6'] = {'acc': 'Q8NH72', 'cn': ' Olfactory receptor 4C6 '}
        GPCR_METADATA['OR4D1'] = {'acc': 'Q15615', 'cn': ' Olfactory receptor 4D1 '}
        GPCR_METADATA['OR4D10'] = {'acc': 'Q8NGI6', 'cn': ' Olfactory receptor 4D10 '}
        GPCR_METADATA['OR4D11'] = {'acc': 'Q8NGI4', 'cn': ' Olfactory receptor 4D11 '}
        GPCR_METADATA['OR4D2'] = {'acc': 'P58180', 'cn': ' Olfactory receptor 4D2 '}
        GPCR_METADATA['OR4D5'] = {'acc': 'Q8NGN0', 'cn': ' Olfactory receptor 4D5 '}
        GPCR_METADATA['OR4D6'] = {'acc': 'Q8NGJ1', 'cn': ' Olfactory receptor 4D6 '}
        GPCR_METADATA['OR4D9'] = {'acc': 'Q8NGE8', 'cn': ' Olfactory receptor 4D9 '}
        GPCR_METADATA['OR4E1'] = {'acc': 'P0C645', 'cn': ' Olfactory receptor 4E1 '}
        GPCR_METADATA['OR4E2'] = {'acc': 'Q8NGC2', 'cn': ' Olfactory receptor 4E2 '}
        GPCR_METADATA['OR4F15'] = {'acc': 'Q8NGB8', 'cn': ' Olfactory receptor 4F15 '}
        GPCR_METADATA['OR4F17'] = {'acc': 'Q8NGA8', 'cn': ' Olfactory receptor 4F17 '}
        GPCR_METADATA['OR4F21'] = {'acc': 'O95013', 'cn': ' Olfactory receptor 4F21 '}
        GPCR_METADATA['OR4F29'] = {'acc': 'Q6IEY1', 'cn': ' Olfactory receptor 4F3/4F16/4F29 '}
        GPCR_METADATA['OR4F4'] = {'acc': 'Q96R69', 'cn': ' Olfactory receptor 4F4 '}
        GPCR_METADATA['OR4F5'] = {'acc': 'Q8NH21', 'cn': ' Olfactory receptor 4F5 '}
        GPCR_METADATA['OR4F6'] = {'acc': 'Q8NGB9', 'cn': ' Olfactory receptor 4F6 '}
        GPCR_METADATA['OR4K1'] = {'acc': 'Q8NGD4', 'cn': ' Olfactory receptor 4K1 '}
        GPCR_METADATA['OR4K13'] = {'acc': 'Q8NH42', 'cn': ' Olfactory receptor 4K13 '}
        GPCR_METADATA['OR4K14'] = {'acc': 'Q8NGD5', 'cn': ' Olfactory receptor 4K14 '}
        GPCR_METADATA['OR4K15'] = {'acc': 'Q8NH41', 'cn': ' Olfactory receptor 4K15 '}
        GPCR_METADATA['OR4K17'] = {'acc': 'Q8NGC6', 'cn': ' Olfactory receptor 4K17 '}
        GPCR_METADATA['OR4K2'] = {'acc': 'Q8NGD2', 'cn': ' Olfactory receptor 4K2 '}
        GPCR_METADATA['OR4K3'] = {'acc': 'Q96R72', 'cn': ' Olfactory receptor 4K3 '}
        GPCR_METADATA['OR4K5'] = {'acc': 'Q8NGD3', 'cn': ' Olfactory receptor 4K5 '}
        GPCR_METADATA['OR4L1'] = {'acc': 'Q8NH43', 'cn': ' Olfactory receptor 4L1 '}
        GPCR_METADATA['OR4M2'] = {'acc': 'Q8NGB6', 'cn': ' Olfactory receptor 4M2 '}
        GPCR_METADATA['OR4M2B'] = {'acc': 'A0A0X1KG70', 'cn': ' Olfactory receptor 4M2B '}
        GPCR_METADATA['OR4N2'] = {'acc': 'Q8NGD1', 'cn': ' Olfactory receptor 4N2 '}
        GPCR_METADATA['OR4N4'] = {'acc': 'Q8N0Y3', 'cn': ' Olfactory receptor 4N4 '}
        GPCR_METADATA['OR4N4C'] = {'acc': 'A0A096LPK9', 'cn': ' Olfactory receptor 4N4C '}
        GPCR_METADATA['OR4N5'] = {'acc': 'Q8IXE1', 'cn': ' Olfactory receptor 4N5 '}
        GPCR_METADATA['OR4P4'] = {'acc': 'Q8NGL7', 'cn': ' Olfactory receptor 4P4 '}
        GPCR_METADATA['OR4Q2'] = {'acc': 'P0C623', 'cn': ' Olfactory receptor 4Q2 '}
        GPCR_METADATA['OR4Q3'] = {'acc': 'Q8NH05', 'cn': ' Olfactory receptor 4Q3 '}
        GPCR_METADATA['OR4S1'] = {'acc': 'Q8NGB4', 'cn': ' Olfactory receptor 4S1 '}
        GPCR_METADATA['OR4S2'] = {'acc': 'Q8NH73', 'cn': ' Olfactory receptor 4S2 '}
        GPCR_METADATA['OR4X1'] = {'acc': 'Q8NH49', 'cn': ' Olfactory receptor 4X1 '}
        GPCR_METADATA['OR4X2'] = {'acc': 'Q8NGF9', 'cn': ' Olfactory receptor 4X2 '}
        GPCR_METADATA['OR51A2'] = {'acc': 'Q8NGJ7', 'cn': ' Olfactory receptor 51A2 '}
        GPCR_METADATA['OR51A4'] = {'acc': 'Q8NGJ6', 'cn': ' Olfactory receptor 51A4 '}
        GPCR_METADATA['OR51A7'] = {'acc': 'Q8NH64', 'cn': ' Olfactory receptor 51A7 '}
        GPCR_METADATA['OR51B2'] = {'acc': 'Q9Y5P1', 'cn': ' Olfactory receptor 51B2 '}
        GPCR_METADATA['OR51B4'] = {'acc': 'Q9Y5P0', 'cn': ' Olfactory receptor 51B4 '}
        GPCR_METADATA['OR51B5'] = {'acc': 'Q9H339', 'cn': ' Olfactory receptor 51B5 '}
        GPCR_METADATA['OR51B6'] = {'acc': 'Q9H340', 'cn': ' Olfactory receptor 51B6 '}
        GPCR_METADATA['OR51D1'] = {'acc': 'Q8NGF3', 'cn': ' Olfactory receptor 51D1 '}
        GPCR_METADATA['OR51E1'] = {'acc': 'Q8TCB6', 'cn': ' Olfactory receptor 51E1 '}
        GPCR_METADATA['OR51E2'] = {'acc': 'Q9H255', 'cn': ' Olfactory receptor 51E2 '}
        GPCR_METADATA['OR51F1'] = {'acc': 'A6NGY5', 'cn': ' Olfactory receptor 51F1 '}
        GPCR_METADATA['OR51F2'] = {'acc': 'Q8NH61', 'cn': ' Olfactory receptor 51F2 '}
        GPCR_METADATA['OR51G1'] = {'acc': 'Q8NGK1', 'cn': ' Olfactory receptor 51G1 '}
        GPCR_METADATA['OR51G2'] = {'acc': 'Q8NGK0', 'cn': ' Olfactory receptor 51G2 '}
        GPCR_METADATA['OR51H1'] = {'acc': 'Q8NH63', 'cn': ' Olfactory receptor 51H1 '}
        GPCR_METADATA['OR51I1'] = {'acc': 'Q9H343', 'cn': ' Olfactory receptor 51I1 '}
        GPCR_METADATA['OR51I2'] = {'acc': 'Q9H344', 'cn': ' Olfactory receptor 51I2 '}
        GPCR_METADATA['OR51J1'] = {'acc': 'Q9H342', 'cn': ' Olfactory receptor 51J1 '}
        GPCR_METADATA['OR51L1'] = {'acc': 'Q8NGJ5', 'cn': ' Olfactory receptor 51L1 '}
        GPCR_METADATA['OR51M1'] = {'acc': 'Q9H341', 'cn': ' Olfactory receptor 51M1 '}
        GPCR_METADATA['OR51Q1'] = {'acc': 'Q8NH59', 'cn': ' Olfactory receptor 51Q1 '}
        GPCR_METADATA['OR51S1'] = {'acc': 'Q8NGJ8', 'cn': ' Olfactory receptor 51S1 '}
        GPCR_METADATA['OR51T1'] = {'acc': 'Q8NGJ9', 'cn': ' Olfactory receptor 51T1 '}
        GPCR_METADATA['OR51V1'] = {'acc': 'Q9H2C8', 'cn': ' Olfactory receptor 51V1 '}
        GPCR_METADATA['OR52A1'] = {'acc': 'Q9UKL2', 'cn': ' Olfactory receptor 52A1 '}
        GPCR_METADATA['OR52A4P'] = {'acc': 'A6NMU1', 'cn': ' Olfactory receptor 52A4 '}
        GPCR_METADATA['OR52A5'] = {'acc': 'Q9H2C5', 'cn': ' Olfactory receptor 52A5 '}
        GPCR_METADATA['OR52B2'] = {'acc': 'Q96RD2', 'cn': ' Olfactory receptor 52B2 '}
        GPCR_METADATA['OR52B4'] = {'acc': 'Q8NGK2', 'cn': ' Olfactory receptor 52B4 '}
        GPCR_METADATA['OR52B6'] = {'acc': 'Q8NGF0', 'cn': ' Olfactory receptor 52B6 '}
        GPCR_METADATA['OR52D1'] = {'acc': 'Q9H346', 'cn': ' Olfactory receptor 52D1 '}
        GPCR_METADATA['OR52E1'] = {'acc': 'Q8NGJ3', 'cn': ' Olfactory receptor 52E1 '}
        GPCR_METADATA['OR52E2'] = {'acc': 'Q8NGJ4', 'cn': ' Olfactory receptor 52E2 '}
        GPCR_METADATA['OR52E4'] = {'acc': 'Q8NGH9', 'cn': ' Olfactory receptor 52E4 '}
        GPCR_METADATA['OR52E5'] = {'acc': 'Q8NH55', 'cn': ' Olfactory receptor 52E5 '}
        GPCR_METADATA['OR52E6'] = {'acc': 'Q96RD3', 'cn': ' Olfactory receptor 52E6 '}
        GPCR_METADATA['OR52E8'] = {'acc': 'Q6IFG1', 'cn': ' Olfactory receptor 52E8 '}
        GPCR_METADATA['OR52H1'] = {'acc': 'Q8NGJ2', 'cn': ' Olfactory receptor 52H1 '}
        GPCR_METADATA['OR52I1'] = {'acc': 'Q8NGK6', 'cn': ' Olfactory receptor 52I1 '}
        GPCR_METADATA['OR52I2'] = {'acc': 'Q8NH67', 'cn': ' Olfactory receptor 52I2 '}
        GPCR_METADATA['OR52J3'] = {'acc': 'Q8NH60', 'cn': ' Olfactory receptor 52J3 '}
        GPCR_METADATA['OR52K1'] = {'acc': 'Q8NGK4', 'cn': ' Olfactory receptor 52K1 '}
        GPCR_METADATA['OR52K2'] = {'acc': 'Q8NGK3', 'cn': ' Olfactory receptor 52K2 '}
        GPCR_METADATA['OR52L1'] = {'acc': 'Q8NGH7', 'cn': ' Olfactory receptor 52L1 '}
        GPCR_METADATA['OR52M1'] = {'acc': 'Q8NGK5', 'cn': ' Olfactory receptor 52M1 '}
        GPCR_METADATA['OR52N1'] = {'acc': 'Q8NH53', 'cn': ' Olfactory receptor 52N1 '}
        GPCR_METADATA['OR52N2'] = {'acc': 'Q8NGI0', 'cn': ' Olfactory receptor 52N2 '}
        GPCR_METADATA['OR52N4'] = {'acc': 'Q8NGI2', 'cn': ' Olfactory receptor 52N4 '}
        GPCR_METADATA['OR52N5'] = {'acc': 'Q8NH56', 'cn': ' Olfactory receptor 52N5 '}
        GPCR_METADATA['OR52P1'] = {'acc': 'Q8NH57', 'cn': ' Olfactory receptor 52P1 '}
        GPCR_METADATA['OR52R1'] = {'acc': 'Q8NGF1', 'cn': ' Olfactory receptor 52R1 '}
        GPCR_METADATA['OR52W1'] = {'acc': 'Q6IF63', 'cn': ' Olfactory receptor 52W1 '}
        GPCR_METADATA['OR52Z1P'] = {'acc': 'P0C646', 'cn': ' Olfactory receptor 52Z1P '}
        GPCR_METADATA['OR56A1'] = {'acc': 'Q8NGH5', 'cn': ' Olfactory receptor 56A1 '}
        GPCR_METADATA['OR56A3'] = {'acc': 'Q8NH54', 'cn': ' Olfactory receptor 56A3 '}
        GPCR_METADATA['OR56A4'] = {'acc': 'Q8NGH8', 'cn': ' Olfactory receptor 56A4 '}
        GPCR_METADATA['OR56A5'] = {'acc': 'P0C7T3', 'cn': ' Olfactory receptor 56A5 '}
        GPCR_METADATA['OR56B1'] = {'acc': 'Q8NGI3', 'cn': ' Olfactory receptor 56B1 '}
        GPCR_METADATA['OR56B4'] = {'acc': 'Q8NH76', 'cn': ' Olfactory receptor 56B4 '}
        GPCR_METADATA['OR5A1'] = {'acc': 'Q8NGJ0', 'cn': ' Olfactory receptor 5A1 '}
        GPCR_METADATA['OR5A2'] = {'acc': 'Q8NGI9', 'cn': ' Olfactory receptor 5A2 '}
        GPCR_METADATA['OR5AC1'] = {'acc': 'P0C628', 'cn': ' Olfactory receptor 5AC1 '}
        GPCR_METADATA['OR5AC2'] = {'acc': 'Q9NZP5', 'cn': ' Olfactory receptor 5AC2 '}
        GPCR_METADATA['OR5AK2'] = {'acc': 'Q8NH90', 'cn': ' Olfactory receptor 5AK2 '}
        GPCR_METADATA['OR5AL1'] = {'acc': 'P0C617', 'cn': ' Olfactory receptor 5AL1 '}
        GPCR_METADATA['OR5AN1'] = {'acc': 'Q8NGI8', 'cn': ' Olfactory receptor 5AN1 '}
        GPCR_METADATA['OR5AP2'] = {'acc': 'Q8NGF4', 'cn': ' Olfactory receptor 5AP2 '}
        GPCR_METADATA['OR5AR1'] = {'acc': 'Q8NGP9', 'cn': ' Olfactory receptor 5AR1 '}
        GPCR_METADATA['OR5AS1'] = {'acc': 'Q8N127', 'cn': ' Olfactory receptor 5AS1 '}
        GPCR_METADATA['OR5AU1'] = {'acc': 'Q8NGC0', 'cn': ' Olfactory receptor 5AU1 '}
        GPCR_METADATA['OR5B12'] = {'acc': 'Q96R08', 'cn': ' Olfactory receptor 5B12 '}
        GPCR_METADATA['OR5B17'] = {'acc': 'Q8NGF7', 'cn': ' Olfactory receptor 5B17 '}
        GPCR_METADATA['OR5B2'] = {'acc': 'Q96R09', 'cn': ' Olfactory receptor 5B2 '}
        GPCR_METADATA['OR5B21'] = {'acc': 'A6NL26', 'cn': ' Olfactory receptor 5B21 '}
        GPCR_METADATA['OR5B3'] = {'acc': 'Q8NH48', 'cn': ' Olfactory receptor 5B3 '}
        GPCR_METADATA['OR5C1'] = {'acc': 'Q8NGR4', 'cn': ' Olfactory receptor 5C1 '}
        GPCR_METADATA['OR5D13'] = {'acc': 'Q8NGL4', 'cn': ' Olfactory receptor 5D13 '}
        GPCR_METADATA['OR5D14'] = {'acc': 'Q8NGL3', 'cn': ' Olfactory receptor 5D14 '}
        GPCR_METADATA['OR5D16'] = {'acc': 'Q8NGK9', 'cn': ' Olfactory receptor 5D16 '}
        GPCR_METADATA['OR5D18'] = {'acc': 'Q8NGL1', 'cn': ' Olfactory receptor 5D18 '}
        GPCR_METADATA['OR5F1'] = {'acc': 'O95221', 'cn': ' Olfactory receptor 5F1 '}
        GPCR_METADATA['OR5G3'] = {'acc': 'P0C626', 'cn': ' Olfactory receptor 5G3 '}
        GPCR_METADATA['OR5H1'] = {'acc': 'A6NKK0', 'cn': ' Olfactory receptor 5H1 '}
        GPCR_METADATA['OR5H14'] = {'acc': 'A6NHG9', 'cn': ' Olfactory receptor 5H14 '}
        GPCR_METADATA['OR5H15'] = {'acc': 'A6NDH6', 'cn': ' Olfactory receptor 5H15 '}
        GPCR_METADATA['OR5H2'] = {'acc': 'Q8NGV7', 'cn': ' Olfactory receptor 5H2 '}
        GPCR_METADATA['OR5H6'] = {'acc': 'Q8NGV6', 'cn': ' Olfactory receptor 5H6 '}
        GPCR_METADATA['OR5H8'] = {'acc': 'P0DN80', 'cn': ' Olfactory receptor 5H8 '}
        GPCR_METADATA['OR5I1'] = {'acc': 'Q13606', 'cn': ' Olfactory receptor 5I1 '}
        GPCR_METADATA['OR5J2'] = {'acc': 'Q8NH18', 'cn': ' Olfactory receptor 5J2 '}
        GPCR_METADATA['OR5K1'] = {'acc': 'Q8NHB7', 'cn': ' Olfactory receptor 5K1 '}
        GPCR_METADATA['OR5K2'] = {'acc': 'Q8NHB8', 'cn': ' Olfactory receptor 5K2 '}
        GPCR_METADATA['OR5K3'] = {'acc': 'A6NET4', 'cn': ' Olfactory receptor 5K3 '}
        GPCR_METADATA['OR5K4'] = {'acc': 'A6NMS3', 'cn': ' Olfactory receptor 5K4 '}
        GPCR_METADATA['OR5L1'] = {'acc': 'Q8NGL2', 'cn': ' Olfactory receptor 5L1 '}
        GPCR_METADATA['OR5L2'] = {'acc': 'Q8NGL0', 'cn': ' Olfactory receptor 5L2 '}
        GPCR_METADATA['OR5M1'] = {'acc': 'Q8NGP8', 'cn': ' Olfactory receptor 5M1 '}
        GPCR_METADATA['OR5M10'] = {'acc': 'Q6IEU7', 'cn': ' Olfactory receptor 5M10 '}
        GPCR_METADATA['OR5M11'] = {'acc': 'Q96RB7', 'cn': ' Olfactory receptor 5M11 '}
        GPCR_METADATA['OR5M3'] = {'acc': 'Q8NGP4', 'cn': ' Olfactory receptor 5M3 '}
        GPCR_METADATA['OR5M8'] = {'acc': 'Q8NGP6', 'cn': ' Olfactory receptor 5M8 '}
        GPCR_METADATA['OR5M9'] = {'acc': 'Q8NGP3', 'cn': ' Olfactory receptor 5M9 '}
        GPCR_METADATA['OR5P2'] = {'acc': 'Q8WZ92', 'cn': ' Olfactory receptor 5P2 '}
        GPCR_METADATA['OR5P3'] = {'acc': 'Q8WZ94', 'cn': ' Olfactory receptor 5P3 '}
        GPCR_METADATA['OR5T1'] = {'acc': 'Q8NG75', 'cn': ' Olfactory receptor 5T1 '}
        GPCR_METADATA['OR5T2'] = {'acc': 'Q8NGG2', 'cn': ' Olfactory receptor 5T2 '}
        GPCR_METADATA['OR5T3'] = {'acc': 'Q8NGG3', 'cn': ' Olfactory receptor 5T3 '}
        GPCR_METADATA['OR5V1'] = {'acc': 'Q9UGF6', 'cn': ' Olfactory receptor 5V1 '}
        GPCR_METADATA['OR5W2'] = {'acc': 'Q8NH69', 'cn': ' Olfactory receptor 5W2 '}
        GPCR_METADATA['OR6A2'] = {'acc': 'O95222', 'cn': ' Olfactory receptor 6A2 '}
        GPCR_METADATA['OR6B1'] = {'acc': 'O95007', 'cn': ' Olfactory receptor 6B1 '}
        GPCR_METADATA['OR6B2'] = {'acc': 'Q6IFH4', 'cn': ' Olfactory receptor 6B2 '}
        GPCR_METADATA['OR6B3'] = {'acc': 'Q8NGW1', 'cn': ' Olfactory receptor 6B3 '}
        GPCR_METADATA['OR6C1'] = {'acc': 'Q96RD1', 'cn': ' Olfactory receptor 6C1 '}
        GPCR_METADATA['OR6C2'] = {'acc': 'Q9NZP2', 'cn': ' Olfactory receptor 6C2 '}
        GPCR_METADATA['OR6C3'] = {'acc': 'Q9NZP0', 'cn': ' Olfactory receptor 6C3 '}
        GPCR_METADATA['OR6C4'] = {'acc': 'Q8NGE1', 'cn': ' Olfactory receptor 6C4 '}
        GPCR_METADATA['OR6C6'] = {'acc': 'A6NF89', 'cn': ' Olfactory receptor 6C6 '}
        GPCR_METADATA['OR6C65'] = {'acc': 'A6NJZ3', 'cn': ' Olfactory receptor 6C65 '}
        GPCR_METADATA['OR6C68'] = {'acc': 'A6NDL8', 'cn': ' Olfactory receptor 6C68 '}
        GPCR_METADATA['OR6C70'] = {'acc': 'A6NIJ9', 'cn': ' Olfactory receptor 6C70 '}
        GPCR_METADATA['OR6C74'] = {'acc': 'A6NCV1', 'cn': ' Olfactory receptor 6C74 '}
        GPCR_METADATA['OR6C75'] = {'acc': 'A6NL08', 'cn': ' Olfactory receptor 6C75 '}
        GPCR_METADATA['OR6C76'] = {'acc': 'A6NM76', 'cn': ' Olfactory receptor 6C76 '}
        GPCR_METADATA['OR6F1'] = {'acc': 'Q8NGZ6', 'cn': ' Olfactory receptor 6F1 '}
        GPCR_METADATA['OR6J1'] = {'acc': 'Q8NGC5', 'cn': ' Olfactory receptor 6J1 '}
        GPCR_METADATA['OR6K2'] = {'acc': 'Q8NGY2', 'cn': ' Olfactory receptor 6K2 '}
        GPCR_METADATA['OR6K3'] = {'acc': 'Q8NGY3', 'cn': ' Olfactory receptor 6K3 '}
        GPCR_METADATA['OR6K6'] = {'acc': 'Q8NGW6', 'cn': ' Olfactory receptor 6K6 '}
        GPCR_METADATA['OR6M1'] = {'acc': 'Q8NGM8', 'cn': ' Olfactory receptor 6M1 '}
        GPCR_METADATA['OR6N1'] = {'acc': 'Q8NGY5', 'cn': ' Olfactory receptor 6N1 '}
        GPCR_METADATA['OR6N2'] = {'acc': 'Q8NGY6', 'cn': ' Olfactory receptor 6N2 '}
        GPCR_METADATA['OR6P1'] = {'acc': 'Q8NGX9', 'cn': ' Olfactory receptor 6P1 '}
        GPCR_METADATA['OR6Q1'] = {'acc': 'Q8NGQ2', 'cn': ' Olfactory receptor 6Q1 '}
        GPCR_METADATA['OR6S1'] = {'acc': 'Q8NH40', 'cn': ' Olfactory receptor 6S1 '}
        GPCR_METADATA['OR6T1'] = {'acc': 'Q8NGN1', 'cn': ' Olfactory receptor 6T1 '}
        GPCR_METADATA['OR6V1'] = {'acc': 'Q8N148', 'cn': ' Olfactory receptor 6V1 '}
        GPCR_METADATA['OR6X1'] = {'acc': 'Q8NH79', 'cn': ' Olfactory receptor 6X1 '}
        GPCR_METADATA['OR6Y1'] = {'acc': 'Q8NGX8', 'cn': ' Olfactory receptor 6Y1 '}
        GPCR_METADATA['OR7A10'] = {'acc': 'O76100', 'cn': ' Olfactory receptor 7A10 '}
        GPCR_METADATA['OR7A17'] = {'acc': 'O14581', 'cn': ' Olfactory receptor 7A17 '}
        GPCR_METADATA['OR7A5'] = {'acc': 'Q15622', 'cn': ' Olfactory receptor 7A5 '}
        GPCR_METADATA['OR7C1'] = {'acc': 'O76099', 'cn': ' Olfactory receptor 7C1 '}
        GPCR_METADATA['OR7C2'] = {'acc': 'O60412', 'cn': ' Olfactory receptor 7C2 '}
        GPCR_METADATA['OR7D2'] = {'acc': 'Q96RA2', 'cn': ' Olfactory receptor 7D2 '}
        GPCR_METADATA['OR7D4'] = {'acc': 'Q8NG98', 'cn': ' Olfactory receptor 7D4 '}
        GPCR_METADATA['OR7E24'] = {'acc': 'Q6IFN5', 'cn': ' Olfactory receptor 7E24 '}
        GPCR_METADATA['OR7G1'] = {'acc': 'Q8NGA0', 'cn': ' Olfactory receptor 7G1 '}
        GPCR_METADATA['OR7G2'] = {'acc': 'Q8NG99', 'cn': ' Olfactory receptor 7G2 '}
        GPCR_METADATA['OR7G3'] = {'acc': 'Q8NG95', 'cn': ' Olfactory receptor 7G3 '}
        GPCR_METADATA['OR8A1'] = {'acc': 'Q8NGG7', 'cn': ' Olfactory receptor 8A1 '}
        GPCR_METADATA['OR8B12'] = {'acc': 'Q8NGG6', 'cn': ' Olfactory receptor 8B12 '}
        GPCR_METADATA['OR8B2'] = {'acc': 'Q96RD0', 'cn': ' Olfactory receptor 8B2 '}
        GPCR_METADATA['OR8B3'] = {'acc': 'Q8NGG8', 'cn': ' Olfactory receptor 8B3 '}
        GPCR_METADATA['OR8B4'] = {'acc': 'Q96RC9', 'cn': ' Olfactory receptor 8B4 '}
        GPCR_METADATA['OR8B8'] = {'acc': 'Q15620', 'cn': ' Olfactory receptor 8B8 '}
        GPCR_METADATA['OR8D1'] = {'acc': 'Q8WZ84', 'cn': ' Olfactory receptor 8D1 '}
        GPCR_METADATA['OR8D2'] = {'acc': 'Q9GZM6', 'cn': ' Olfactory receptor 8D2 '}
        GPCR_METADATA['OR8D4'] = {'acc': 'Q8NGM9', 'cn': ' Olfactory receptor 8D4 '}
        GPCR_METADATA['OR8G1'] = {'acc': 'Q15617', 'cn': ' Olfactory receptor 8G1 '}
        GPCR_METADATA['OR8G5'] = {'acc': 'Q8NG78', 'cn': ' Olfactory receptor 8G5 '}
        GPCR_METADATA['OR8H1'] = {'acc': 'Q8NGG4', 'cn': ' Olfactory receptor 8H1 '}
        GPCR_METADATA['OR8H2'] = {'acc': 'Q8N162', 'cn': ' Olfactory receptor 8H2 '}
        GPCR_METADATA['OR8H3'] = {'acc': 'Q8N146', 'cn': ' Olfactory receptor 8H3 '}
        GPCR_METADATA['OR8I2'] = {'acc': 'Q8N0Y5', 'cn': ' Olfactory receptor 8I2 '}
        GPCR_METADATA['OR8J1'] = {'acc': 'Q8NGP2', 'cn': ' Olfactory receptor 8J1 '}
        GPCR_METADATA['OR8J2'] = {'acc': 'Q8NGG1', 'cn': ' Olfactory receptor 8J2 '}
        GPCR_METADATA['OR8J3'] = {'acc': 'Q8NGG0', 'cn': ' Olfactory receptor 8J3 '}
        GPCR_METADATA['OR8K1'] = {'acc': 'Q8NGG5', 'cn': ' Olfactory receptor 8K1 '}
        GPCR_METADATA['OR8K3'] = {'acc': 'Q8NH51', 'cn': ' Olfactory receptor 8K3 '}
        GPCR_METADATA['OR8K5'] = {'acc': 'Q8NH50', 'cn': ' Olfactory receptor 8K5 '}
        GPCR_METADATA['OR8S1'] = {'acc': 'Q8NH09', 'cn': ' Olfactory receptor 8S1 '}
        GPCR_METADATA['OR8U1'] = {'acc': 'Q8NH10', 'cn': ' Olfactory receptor 8U1 '}
        GPCR_METADATA['OR8U3'] = {'acc': 'Q8NH85', 'cn': ' Olfactory receptor 8U3 '}
        GPCR_METADATA['OR8U8'] = {'acc': 'P0C7N1', 'cn': ' Olfactory receptor 8U8 '}
        GPCR_METADATA['OR8U9'] = {'acc': 'P0C7N5', 'cn': ' Olfactory receptor 8U9 '}
        GPCR_METADATA['OR9A1P'] = {'acc': 'Q8NGU1', 'cn': ' Olfactory receptor 9A1 '}
        GPCR_METADATA['OR9A2'] = {'acc': 'Q8NGT5', 'cn': ' Olfactory receptor 9A2 '}
        GPCR_METADATA['OR9A4'] = {'acc': 'Q8NGU2', 'cn': ' Olfactory receptor 9A4 '}
        GPCR_METADATA['OR9G1'] = {'acc': 'Q8NH87', 'cn': ' Olfactory receptor 9G1 '}
        GPCR_METADATA['OR9G4'] = {'acc': 'Q8NGQ1', 'cn': ' Olfactory receptor 9G4 '}
        GPCR_METADATA['OR9G9'] = {'acc': 'P0C7N8', 'cn': ' Olfactory receptor 9G9 '}
        GPCR_METADATA['OR9I1'] = {'acc': 'Q8NGQ6', 'cn': ' Olfactory receptor 9I1 '}
        GPCR_METADATA['OR9K2'] = {'acc': 'Q8NGE7', 'cn': ' Olfactory receptor 9K2 '}
        GPCR_METADATA['OR9Q1'] = {'acc': 'Q8NGQ5', 'cn': ' Olfactory receptor 9Q1 '}
        GPCR_METADATA['OR9Q2'] = {'acc': 'Q8NGE9', 'cn': ' Olfactory receptor 9Q2 '}
        GPCR_METADATA['OXER1'] = {'acc': 'Q8TDS5', 'cn': ' Oxoeicosanoid receptor 1 '}
        GPCR_METADATA['P2RY8'] = {'acc': 'Q86VZ1', 'cn': ' S-geranylgeranyl-glutathione receptor P2RY8 '}
        GPCR_METADATA['PROKR1'] = {'acc': 'Q8TCW9', 'cn': ' Prokineticin receptor 1 '}
        GPCR_METADATA['PROKR2'] = {'acc': 'Q8NFJ6', 'cn': ' Prokineticin receptor 2 '}
        GPCR_METADATA['PTGDR2'] = {'acc': 'Q9Y5Y4', 'cn': ' Prostaglandin D2 receptor 2 '}
        GPCR_METADATA['QRFPR'] = {'acc': 'Q96P65', 'cn': ' Pyroglutamylated RF-amide peptide receptor '}
        GPCR_METADATA['RHO'] = {'acc': 'P08100', 'cn': ' Rhodopsin '}
        GPCR_METADATA['RXFP1'] = {'acc': 'Q9HBX9', 'cn': ' Relaxin receptor 1 '}
        GPCR_METADATA['RXFP2'] = {'acc': 'Q8WXD0', 'cn': ' Relaxin receptor 2 '}
        GPCR_METADATA['RXFP3'] = {'acc': 'Q9NSD7', 'cn': ' Relaxin-3 receptor 1 '}
        GPCR_METADATA['RXFP4'] = {'acc': 'Q8TDU9', 'cn': ' Relaxin-3 receptor 2 '}
        GPCR_METADATA['S1PR4'] = {'acc': 'O95977', 'cn': ' Sphingosine 1-phosphate receptor 4 '}
        GPCR_METADATA['SCTR'] = {'acc': 'P47872', 'cn': ' Secretin receptor '}
        GPCR_METADATA['SMO'] = {'acc': 'Q99835', 'cn': ' Protein smoothened '}
        GPCR_METADATA['SUCNR1'] = {'acc': 'Q9BXA5', 'cn': ' Succinate receptor 1 '}
        GPCR_METADATA['TAAR1'] = {'acc': 'Q96RJ0', 'cn': ' Trace amine-associated receptor 1 '}
        GPCR_METADATA['TAAR2'] = {'acc': 'Q9P1P5', 'cn': ' Trace amine-associated receptor 2 '}
        GPCR_METADATA['TAAR5'] = {'acc': 'O14804', 'cn': ' Trace amine-associated receptor 5 '}
        GPCR_METADATA['TAAR6'] = {'acc': 'Q96RI8', 'cn': ' Trace amine-associated receptor 6 '}
        GPCR_METADATA['TAAR8'] = {'acc': 'Q969N4', 'cn': ' Trace amine-associated receptor 8 '}
        GPCR_METADATA['TAAR9'] = {'acc': 'Q96RI9', 'cn': ' Trace amine-associated receptor 9 '}
        GPCR_METADATA['TAS1R1'] = {'acc': 'Q7RTX1', 'cn': ' Taste receptor type 1 member 1 '}
        GPCR_METADATA['TAS1R2'] = {'acc': 'Q8TE23', 'cn': ' Taste receptor type 1 member 2 '}
        GPCR_METADATA['TAS1R3'] = {'acc': 'Q7RTX0', 'cn': ' Taste receptor type 1 member 3 '}
        GPCR_METADATA['TAS2R1'] = {'acc': 'Q9NYW7', 'cn': ' Taste receptor type 2 member 1 '}
        GPCR_METADATA['TAS2R10'] = {'acc': 'Q9NYW0', 'cn': ' Taste receptor type 2 member 10 '}
        GPCR_METADATA['TAS2R13'] = {'acc': 'Q9NYV9', 'cn': ' Taste receptor type 2 member 13 '}
        GPCR_METADATA['TAS2R14'] = {'acc': 'Q9NYV8', 'cn': ' Taste receptor type 2 member 14 '}
        GPCR_METADATA['TAS2R16'] = {'acc': 'Q9NYV7', 'cn': ' Taste receptor type 2 member 16 '}
        GPCR_METADATA['TAS2R19'] = {'acc': 'P59542', 'cn': ' Taste receptor type 2 member 19 '}
        GPCR_METADATA['TAS2R20'] = {'acc': 'P59543', 'cn': ' Taste receptor type 2 member 20 '}
        GPCR_METADATA['TAS2R3'] = {'acc': 'Q9NYW6', 'cn': ' Taste receptor type 2 member 3 '}
        GPCR_METADATA['TAS2R30'] = {'acc': 'P59541', 'cn': ' Taste receptor type 2 member 30 '}
        GPCR_METADATA['TAS2R31'] = {'acc': 'P59538', 'cn': ' Taste receptor type 2 member 31 '}
        GPCR_METADATA['TAS2R38'] = {'acc': 'P59533', 'cn': ' Taste receptor type 2 member 38 '}
        GPCR_METADATA['TAS2R39'] = {'acc': 'P59534', 'cn': ' Taste receptor type 2 member 39 '}
        GPCR_METADATA['TAS2R4'] = {'acc': 'Q9NYW5', 'cn': ' Taste receptor type 2 member 4 '}
        GPCR_METADATA['TAS2R40'] = {'acc': 'P59535', 'cn': ' Taste receptor type 2 member 40 '}
        GPCR_METADATA['TAS2R41'] = {'acc': 'P59536', 'cn': ' Taste receptor type 2 member 41 '}
        GPCR_METADATA['TAS2R42'] = {'acc': 'Q7RTR8', 'cn': ' Taste receptor type 2 member 42 '}
        GPCR_METADATA['TAS2R43'] = {'acc': 'P59537', 'cn': ' Taste receptor type 2 member 43 '}
        GPCR_METADATA['TAS2R45'] = {'acc': 'P59539', 'cn': ' Taste receptor type 2 member 45 '}
        GPCR_METADATA['TAS2R46'] = {'acc': 'P59540', 'cn': ' Taste receptor type 2 member 46 '}
        GPCR_METADATA['TAS2R5'] = {'acc': 'Q9NYW4', 'cn': ' Taste receptor type 2 member 5 '}
        GPCR_METADATA['TAS2R50'] = {'acc': 'P59544', 'cn': ' Taste receptor type 2 member 50 '}
        GPCR_METADATA['TAS2R60'] = {'acc': 'P59551', 'cn': ' Taste receptor type 2 member 60 '}
        GPCR_METADATA['TAS2R7'] = {'acc': 'Q9NYW3', 'cn': ' Taste receptor type 2 member 7 '}
        GPCR_METADATA['TAS2R8'] = {'acc': 'Q9NYW2', 'cn': ' Taste receptor type 2 member 8 '}
        GPCR_METADATA['TAS2R9'] = {'acc': 'Q9NYW1', 'cn': ' Taste receptor type 2 member 9 '}
        GPCR_METADATA['TPRA1'] = {'acc': 'Q86W33', 'cn': ' Transmembrane protein adipocyte-associated 1 '}
        GPCR_METADATA['TRHR'] = {'acc': 'P34981', 'cn': ' Thyrotropin-releasing hormone receptor '}
        GPCR_METADATA['TSHR'] = {'acc': 'P16473', 'cn': ' Thyrotropin receptor '}
        GPCR_METADATA['VIPR1'] = {'acc': 'P32241', 'cn': ' Vasoactive intestinal polypeptide receptor 1 '}
        GPCR_METADATA['VIPR2'] = {'acc': 'P41587', 'cn': ' Vasoactive intestinal polypeptide receptor 2 '}
        GPCR_METADATA['XCR1'] = {'acc': 'P46094', 'cn': ' Chemokine XC receptor 1 '}

        organized_data = {}

        for family, members in GPCR_FAMILIES_COMPLEX.items():
            present_members = []
            for m in members:
                common_name = ""
                if m in GPCR_METADATA and 'cn' in GPCR_METADATA[m]:
                    common_name = str(GPCR_METADATA[m]['cn']).strip()

                    short_name = PROTEIN_NAMES.get(m, "")

                    full_search_text = f"{m} {short_name} {common_name}"

                    coupling_info = COUPLING_DATA.get(m, [])

                    present_members.append({
                        "id": m,
                        "text": m,
                        "protein": short_name,
                        "desc": common_name,
                        "search_text": full_search_text,
                        "coupling": coupling_info
                    })

            if present_members:
                organized_data[family] = present_members

        return jsonify(organized_data)

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route('/api/get-pdf/<receptor_name>', methods=['GET'])
def get_pdf(receptor_name):
    filename = f"{receptor_name}.pdf"
    file_path = os.path.join(PDF_FOLDER, filename)

    if os.path.exists(file_path):
        return send_file(file_path, mimetype='application/pdf')
    else:
        return "File not found", 404


@app.route('/api/receptor-details/<receptor_name>')
def get_receptor_details(receptor_name):
    data = {
        "common_name": "",
        "uniprot": "",
        "protein": "",
        "coupling": []
    }

    if receptor_name in GPCR_METADATA:
        info = GPCR_METADATA[receptor_name]
        data["common_name"] = info.get('cn', '').strip()
        data["uniprot"] = info.get('acc', '').strip()

    if receptor_name in PROTEIN_NAMES:
        data["protein"] = PROTEIN_NAMES.get(receptor_name, "")

    if receptor_name in COUPLING_DATA:
        data["coupling"] = COUPLING_DATA[receptor_name]

    return jsonify(data)


mapped_genes = pd.read_csv(path + "/static/mapped_genes.tsv", sep='\t')


@app.route('/api/structure/<receptor_id>')
def get_structure(receptor_id):
    try:

        if "_" not in receptor_id:
            abort(400)

        parts = receptor_id.split("_")
        receptor_raw = parts[0]
        gprotein_raw = parts[1]

        r_match = mapped_genes[mapped_genes["to"] == receptor_raw]["from"]
        g_match = mapped_genes[mapped_genes["to"] == gprotein_raw]["from"]

        if r_match.empty or g_match.empty:
            abort(404)

        receptor_mapped = r_match.iloc[0]
        gprotein_mapped = g_match.iloc[0]

        filename = f"AF-{receptor_mapped}-{gprotein_mapped}.cif"

        return send_from_directory('static/complex', filename)

    except Exception as e:
        print(f"Error: {e}")
        abort(500)


precog3d = pd.read_csv(path + "/static/merged_coupling_values.tsv", sep='\t')


@app.route('/api/precog3d/<receptor_id>')
def get_precog3d_api_data(
        receptor_id):
    try:
        receptor = receptor_id.split(" ")[0].strip()
        result = precog3d[precog3d["GPCR_name"] == receptor]

        if not result.empty:
            return jsonify(result.iloc[0].to_dict())
        else:
            return jsonify({"error": "Receptor not found"}), 404

    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == '__main__':
    app.run(debug=True)
