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
            if '.' in position:
                base = position.split('.')[0]
                value = int(position.split('.')[1])
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

    #print (data)
    df = pd.DataFrame(data, columns = ['position', 'order', 'base', 'value'])
    df = df.sort_values(['order','value'],ascending=[True, True])
    positions = []
    for row in df.to_numpy():
        positions.append(row[0])

    return (np.array(positions))

def extract_contacts(gprotein_given, cutoff):
    assay = ''
    for line in open(path+'/data/contacts/gprotein_best_layer.txt', 'r'):
        if gprotein_given == line.split('\t')[0]:
            assay = line.split('\t')[1]
    #print (gprotein_given)
    dic = {}; positions = []; pair_positions = []; scores = [];
    for line in open(path+'/data/contacts/all_positions_count_'+assay+'_scaled_web.txt', 'r'):
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

#
@app.route('/fetchContactsHeatmap', methods=['GET', 'POST'])
def fetchContactsHeatmap():
    if request.method == 'POST':
        data = request.get_json(force=True)
        #print (data['gpcr'])
        gprotein_given = data['gprotein']
        cutoff = float(data['cutoff'])
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff)
        for i in range(0, len(positions)):
            positions[i] = str(positions[i]).replace('.', 'x')
        return jsonify({'fetch_contactsMin': scoresMin, 'fetch_contactsMax': scoresMax, 'fetch_contacts': scores, 'positions': positions.tolist()})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

def extract_pca(gprotein, assay, pca_type):
    if pca_type != 'GPCRome':
        Xs_train_pca = np.load(path+'/static/best_PCA/'+gprotein+'.npy', allow_pickle=True)
    else:
        #Xs_train_pca = np.load(path+'/static/33layer_PCA/33layer.npy', allow_pickle=True)
        Xs_train_pca = np.load(path+'/static/best_PCA/GNAZ.npy', allow_pickle=True)
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
    if assay == 'Shedding':
        num = -1
        for line in open(path+'/data/shedding.tsv', 'r'):
            if line[0] != '#':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                score = float(line.split('\t')[num+2]) + 1
                #color = cm.get_cmap('RdYlGn', 100)
                #r,g,b,a = color(score)
                #print (score,r,g,b)
                if float(line.split('\t')[num+2]) >= -1.0:
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
                header = line.replace('\n', '').split('\t')[2:]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        flag = 1
                        break
                if flag == 0:
                    num = -1

    elif assay == 'ebBRET':
        num = -1
        for line in open(path+'/data/ebbret.tsv', 'r', encoding="utf-8"):
            if line[0] != '#':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                score = float(line.split('\t')[num+2])

                if float(line.split('\t')[num+2]) > 0.0:
                    genes_to_consider_coupling.append(gene+'|'+acc)
                    color = cm.get_cmap('Greens', 100)
                    r,g,b,a = color(score)
                    score_coupling.append('rgb('+str(r)+','+str(g)+','+str(b)+')')
                else:
                    genes_to_consider_uncoupling.append(gene+'|'+acc)
                    #score_uncoupling.append('rgb('+str(r)+','+str(g)+','+str(b)+')')
                    score_uncoupling.append('grey')
            else:
                header = line.replace('\n', '').split('\t')[2:]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        #print (gprot)
                        break

    elif assay == 'IUPHAR':
        iuphar_map = {
                      'GNAS': 'Gs', 'GNAL': 'Gs',
                      'GNAI1': 'Gi/Go', 'GNAI2': 'Gi/Go', 'GNAI3': 'Gi/Go', 'GNAO1': 'Gi/Go', 'GNAZ': 'Gi/Go', 'GoA': 'Gi/Go', 'GoB': 'Gi/Go',
                      'GNA12': 'G12/G13', 'GNA13': 'G12/G13',
                      'GNAQ': 'Gq/G11', 'GNA11': 'Gq/G11', 'GNA14': 'Gq/G11', 'GNA15': 'Gq/G11'
                      }
        gprotein_fam = iuphar_map[gprotein]
        for line in open(path+'/data/iuphar.tsv', 'r'):
            if line[0] != '#' and line.split('\t')[1] != '':
                gene = line.split('\t')[0]
                acc = line.split('\t')[1]
                if gprotein_fam in line:
                    genes_to_consider_coupling.append(gene+'|'+acc)
                    if gprotein_fam in line.split('\t')[2]:
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
            if barr in line:
                genes_to_consider_coupling.append(gene+'|'+acc)
                score_coupling.append('green')
            else:
                genes_to_consider_uncoupling.append(gene+'|'+acc)
                score_uncoupling.append('green')
    #print (genes_to_consider_coupling)
    #print (genes_to_consider_uncoupling)

    gpcr_list = []
    for line in open(path+'/static/gpcr_list_new_GN.txt', 'r'):
        gene = line.replace('\n', '').split('\t')[1]
        acc = line.replace('\n', '').split('\t')[0]
        gpcr_list.append(gene + '|' + acc)
        #gpcr_list.append(line.replace('\n', '').split('\t')[1])

    X_pos = []
    X_neg = []
    X_grey = []
    genes_to_consider_grey = []
    for gene, row in zip(gpcr_list, X):
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
        ebBRET = ['GNAS', 'GNAI1', 'GNAI2', 'GoA', 'GoB', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA11', 'GNA14', 'GNA15', 'Barr1-GRK2', 'Barr2', 'Barr2-GRK2']
        shedding = ['GNAS', 'GNAL', 'GNAI1', 'GNAI3', 'GNAO1', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15']
        both = ['GNAS', 'GNAI1', 'GNAZ', 'GNA12', 'GNA13', 'GNAQ', 'GNA14', 'GNA15']

        if 'Barr' in gprotein_given:
            assay = 'ebBRET'
            assayList = ['ebBRET', 'STRING']
        elif gprotein_given in both:
            assay = 'Shedding'
            assayList = ['Shedding', 'ebBRET', 'IUPHAR']
        elif gprotein_given in shedding:
            assay = 'Shedding'
            assayList = ['Shedding', 'IUPHAR']
        else:
            assay = 'ebBRET'
            assayList = ['ebBRET', 'IUPHAR']

        if assay_given in assayList:
            assay = assay_given

        ### MUT
        Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+gpcr_given+'.npy', allow_pickle=True)
        #print ('test',Xs_test_pca)
        x_test = Xs_test_pca[0].tolist()
        y_test = Xs_test_pca[1].tolist()
        #print (x_test)
        #print (y_test)
        ### WT
        if '_WT' not in gpcr_given:
            wt = gpcr_given.split('_')[0] + '_WT'
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+wt+'.npy', allow_pickle=True)
            #print (Xs_wt_pca)
            #x_test = Xs_test_pca[:,0].tolist()
            #y_test = Xs_test_pca[:,1].tolist()
            x_wt = Xs_wt_pca[0].tolist()
            y_wt = Xs_wt_pca[1].tolist()
            #print (x_wt)
            #print (y_wt)
        else:
            x_wt = '-'
            y_wt = '-'

        score_coupling, score_uncoupling, x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = extract_pca(gprotein_given, assay, pca_type)
        #print (x_train, y_train, x_test, y_test)
        #print (genes_to_consider_coupling)
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
                        'x_wt': x_wt,
                        'y_wt': y_wt,
                        'minX': minX,
                        'maxX': maxX,
                        'minY': minY,
                        'maxY': maxY})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

def DoBLAST(uniq_id, gpcr_given):
    handle = open(path + "/static/predictor/output/"+uniq_id+"/GPCRDBblast.txt", 'r')
    blast_records = NCBIXML.parse(handle)
    #print (blast_records)

    GPCRDB2SEQ = {}
    for blast_record in blast_records:
        #print (blast_record.query)
        if gpcr_given == blast_record.query:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    bestHIT = alignment.title.split(' ')[1]
                    for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                        if q!='-' and s!='-':
                            GPCRDB2SEQ[num + hsp.sbjct_start] = num + hsp.query_start
                    break
                break
            #print (bestHIT)
            #print (GPCRDB2SEQ)
            break
    return (GPCRDB2SEQ, bestHIT)

@app.route('/fetchContactsSequence', methods=['GET', 'POST'])
def fetchContactsSequence():
    if request.method == 'POST':
        data = request.get_json(force=True)
        #print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        print (gpcr_given, gpcr_given.split('_')[1])
        path_to_fasta = data['path_to_fasta']
        uniq_id = data['uniq_id']
        cutoff = float(data['cutoff'])
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff)

        #print ('fetch_seq', positions)
        fasta_sequence = ''; flag = 0
        for line in open(path_to_fasta):
            if line[0] == '>':
                if flag == 1:
                    break
                flag = 0
                gpcr_found = line.split('>')[1].replace('\n', '').replace(' ', '')
                if gpcr_found == gpcr_given:
                    flag = 1
            elif flag == 1:
                fasta_sequence += line.replace('\n', '')

        GPCRDB2SEQ, bestHIT = DoBLAST(uniq_id, gpcr_given)
        bestHIT_ACC = bestHIT.split('|')[1]

        BW2GPCRDB = {}
        for line in open(path + '/data/GPCRDB/GPCRDB.tsv', 'r'):
            if 'Name' not in line.split('\t')[0]:
                acc = line.split('\t')[0].split('_')[1]
                if acc == bestHIT_ACC:
                    GPCRDB = int(line.split('\t')[1][1:])
                    BW = line.split('\t')[2]
                    BW2GPCRDB[BW] = GPCRDB

        #print (BW2GPCRDB)
        #print (positions)
        seq_positions = []
        bw_positions = []
        for BW in positions:
            if BW in BW2GPCRDB:
                GPCRDB = BW2GPCRDB[BW]
                SEQ = GPCRDB2SEQ[GPCRDB]
                seq_positions.append(int(SEQ))
                bw_positions.append(BW)
            else:
                print (BW)

        #print (list(set(seq_positions)))
        #seq_positions = list(set(seq_positions))
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

        new_seq_positions = []
        new_bw_positions = []
        mutation_position = gpcr_given.split('_')[1]
        if mutation_position != 'WT':
            position = int(mutation_position[1:-1])
            if position not in seq_positions:
                for i in range(0, len(seq_positions)):
                    if seq_positions[i] < position and seq_positions[i+1] > position:
                        #print (seq_positions[i], position, seq_positions[i+1])
                        new_seq_positions.append(position)
                        new_bw_positions.append('-')
                    else:
                        new_seq_positions.append(seq_positions[i])
                        new_bw_positions.append(bw_positions[i])

            seq_positions = new_seq_positions
            bw_positions = new_bw_positions

        return jsonify({'fetch_contacts': scores,
                        'seq_positions': seq_positions,
                        'bw_positions': bw_positions,
                        'sequence': fasta_sequence})
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
        #print (pdbID)
        #print (pair_positions)

        GPCRDB2PDB = {}
        for line in open(path + '/data/PDB/GPCRDB/'+pdbID+'.txt', 'r'):
            GPCRDB2PDB[int(line.split('\t')[3].replace('\n', ''))] = int(line.split('\t')[2])
            bestHIT_ACC = line.replace('\n', '').split('\t')[4].split('|')[1]

        BW2GPCRDB = {}
        for line in open(path + '/data/GPCRDB/GPCRDB.tsv', 'r'):
            if 'Name' not in line.split('\t')[0]:
                acc = line.split('\t')[0].split('_')[1]
                if acc == bestHIT_ACC:
                    GPCRDB = int(line.split('\t')[1][1:])
                    BW = line.split('\t')[2]
                    BW2GPCRDB[BW] = GPCRDB

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
        if '_WT' not in gpcr_given:
            mutation_sequence_position = int(gpcr_given.split('_')[1][1:-1])
            handle = open(path + "/static/predictor/output/"+uniq_id+"/GPCRDBblast.txt", 'r')
            blast_records = NCBIXML.parse(handle)
            SEQ2GPCRDB = {}
            for blast_record in blast_records:
                #print (blast_record.query)
                if gpcr_given == blast_record.query:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            bestHIT = alignment.title.split(' ')[1]
                            for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                                if q!='-' and s!='-':
                                    SEQ2GPCRDB[num + hsp.query_start] = num + hsp.sbjct_start
                            break
                        break
                    break

            if mutation_sequence_position in SEQ2GPCRDB:
                mutation_GPCRDB_position = SEQ2GPCRDB[mutation_sequence_position]
                mutation_position_label = mutation_sequence_position

            if mutation_GPCRDB_position in GPCRDB2PDB:
                mutation_position = GPCRDB2PDB[mutation_GPCRDB_position]

        #print (modified_positions)
        #print (modified_pair_positions)

        return jsonify({'modified_positions': '_'.join(modified_positions),
                        'modified_positions_labels': '_'.join(modified_positions_labels),
                        'modified_num_contacts': '_'.join(modified_num_contacts),
                        'modified_pair_positions': '_'.join(modified_pair_positions),
                        'mutation_position': mutation_position,
                        'mutation_position_label': mutation_position_label
                        })
    else:
        return ("<html><h3>It was a GET request</h3></html>")

@app.route('/fetchContactsPDBStructure', methods=['GET', 'POST'])
def fetchContactsPDBStructure():
    if request.method == 'POST':
        data = request.get_json(force=True)
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        cutoff = float(data['cutoff'])
        uniq_id = data['uniq_id']
        scoresMax, scoresMin, scores, positions, pair_positions, num_contacts = extract_contacts(gprotein_given, cutoff)
        #print ('print', gprotein_given, pair_positions)
        ordered_pdbs = reorder_pdbs(uniq_id, gpcr_given, gprotein_given) ## return list of reordered PDB IDs based on GPCR
        return jsonify({'try': positions.tolist(),
                        'ordered_pdbs': ordered_pdbs,
                        'positions': ','.join(positions.tolist()),
                        'num_contacts': ','.join(num_contacts),
                        'pair_positions': ','.join(pair_positions)})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

## Function to return list of PDB IDs based on GPCR using BLAST
def reorder_pdbs(uniq_id, gpcr, gprotein):
    path_to_fasta = path+"/static/predictor/output/"+uniq_id+"/input.fasta"
    path_to_output = path+"/static/predictor/output/"+uniq_id+"/"

    os.system('blastp -query '+path_to_fasta+' -db '+ path + '/data/PDB/blastdb/allPDB -out '+path_to_output+'/blastp_output.txt')

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
            name = line.split('Query=')[1].replace('\n', '').replace(' ', '')
            dic[name] = []
        elif line[0] == '>':
            pdbid = line.split('>')[1].split('|')[0].split('_')[0].lower()
            if pdbid in chain_info:
                row = []
                #print (chain_info[pdbid])
                row.append(pdbid)
                row.append(chain_info[pdbid]['gpcr_chain'])
                row.append(chain_info[pdbid]['gprotein_chain'])
                #dic[name].append(row)
                dic[name].append(pdbid+'_'+chain_info[pdbid]['gpcr_chain']+'_'+chain_info[pdbid]['gprotein_chain'])

    return(dic[gpcr])
    #return None

@app.route('/bestGprotein', methods=['GET', 'POST'])
def bestGprotein():
    if request.method == 'POST':
        data = request.get_json(force=True)
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        #print ('Given', gpcr_given)
        gpcr = gpcr_given.split('_')[0]
        variant = gpcr_given.split('_')[1]
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
        #print (bestGprotein)
        #print (colIndex)
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
            uniq_id = precogx.main(15, input, input_file, 'all', app.root_path)
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
                if key2[1] == 'WT' or key2[1] not in ['IUPHAR', 'LogRAi', 'Emax']:
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
        return render_template('result.html',
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

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
