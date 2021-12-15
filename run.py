from flask import Flask, render_template, request, jsonify, redirect, url_for
import os, sys, json
import numpy as np
from sklearn.preprocessing import MinMaxScaler
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

def extract_contacts(gprotein_given, cutoff):
    dic = {}; positions = []; pair_positions = []; scores = []
    for line in open(path+'/static/Gprot_contacts/specific_position_score.txt', 'r'):
        if line.split(' ')[0] != 'TM_pos1':
            pos1 = line.split(' ')[0]
            pos2 = line.split(' ')[1]
            gprotein_found = line.split(' ')[2]
            score = float(line.replace('\n','').split(' ')[3])
            if gprotein_given == gprotein_found:
                if score >= cutoff:
                    #print (score)
                    if pos1 not in dic:
                        dic[pos1] = {}
                    dic[pos1][pos2] = score

                    if pos2 not in dic:
                        dic[pos2] = {}
                    dic[pos2][pos1] = score

                    positions.append(pos1)
                    positions.append(pos2)
                    pair_positions.append(pos1+':'+pos2+':'+str(score))
                scores.append(score)

    scoresMax = max(scores)
    scoresMin = min(scores)
    positions = list(set(positions))
    positions = np.array(np.sort(positions))
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

    scaler = MinMaxScaler(feature_range=(0.25,1.0))
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
        Xs_train_pca = np.load(path+'/static/33layer_PCA/33layer.npy', allow_pickle=True)
    Xs_train_pca_coupling, Xs_train_pca_uncoupling, Xs_train_pca_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = filter_gpcr_list(Xs_train_pca, assay, gprotein)
    #print ('train', Xs_train_pca_coupling)
    x_train_coupling = Xs_train_pca_coupling[:,0].tolist()
    x_train_uncoupling = Xs_train_pca_uncoupling[:,0].tolist()
    x_train_grey = Xs_train_pca_grey[:,0].tolist()
    y_train_coupling = Xs_train_pca_coupling[:,1].tolist()
    y_train_uncoupling = Xs_train_pca_uncoupling[:,1].tolist()
    y_train_grey = Xs_train_pca_grey[:,1].tolist()
    return x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey

def filter_gpcr_list(X, assay, gprotein):
    genes_to_consider_coupling = []
    genes_to_consider_uncoupling = []
    #print (assay)
    #assay = 'ebBRET'
    if assay == 'Shedding':
        num = -1
        for line in open(path+'/static/predictor/data_precog/LogRAi_values_final.tsv', 'r'):
            if line[0] != '#':
                gene = line.split('\t')[0]
                if float(line.split('\t')[num+1]) >= -1.0:
                    genes_to_consider_coupling.append(gene)
                else:
                    genes_to_consider_uncoupling.append(gene)
            else:
                header = line.replace('\n', '').split('\t')[1:]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        break

    elif assay == 'ebBRET':
        num = -1
        for line in open(path+'/static/predictor/data_precog2/emax.tsv', 'r', encoding="utf-8"):
            if line[0] != '#':
                gene = line.split('\t')[2]
                if float(line.split('\t')[num+5]) > 0.0:
                    genes_to_consider_coupling.append(gene)
                else:
                    genes_to_consider_uncoupling.append(gene)
            else:
                header = line.replace('\n', '').split('\t')[5:-3]
                for num, gprot in enumerate(header):
                    if gprot == gprotein:
                        #print (gprot)
                        break
    elif assay == 'IUPHAR':
        iuphar_map = {
                      'GNAS': 'Gs', 'GNAL': 'Gs',
                      'GNAI1': 'Gi/Go', 'GNAI2': 'Gi/Go', 'GNAI3': 'Gi/Go', 'GNAO1': 'Gi/Go', 'GNAZ': 'Gi/Go',
                      'GNA12': 'G12/G13', 'GNA13': 'G12/G13',
                      'GNAQ': 'Gq/G11', 'GNA11': 'Gq/G11', 'GNA14': 'Gq/G11', 'GNA15': 'Gq/G11'
                      }
        gprotein_fam = iuphar_map[gprotein]
        for line in open(path+'/static/predictor/data_precog2/IUPHAR_couplings.tsv', 'r'):
            if line[0] != '#' and line.split('\t')[1] != '':
                gene = line.split('\t')[1]
                if gprotein_fam in line:
                    genes_to_consider_coupling.append(gene)
                else:
                    genes_to_consider_uncoupling.append(gene)
    #print (genes_to_consider_coupling)
    #print (genes_to_consider_uncoupling)

    gpcr_list = []
    for line in open(path+'/static/gpcr_list_new_GN.txt', 'r'):
        gpcr_list.append(line.replace('\n', '').split('\t')[1])
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

    return (np.array(X_pos), np.array(X_neg), np.array(X_grey), genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey)

@app.route('/fetchPCA', methods=['GET', 'POST'])
def fetchPCA():
    if request.method == 'POST':
        data = request.get_json(force=True)
        assay = data['assay']
        pca_type = data['pca_type']
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        #print ('pca_type', gpcr_given)
        uniq_id = data['uniq_id']
        ### MUT
        Xs_test_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+gpcr_given+'.npy', allow_pickle=True)
        #print ('test',Xs_test_pca)
        #x_test = Xs_test_pca[:,0].tolist()
        #y_test = Xs_test_pca[:,1].tolist()
        x_test = Xs_test_pca[0].tolist()
        y_test = Xs_test_pca[1].tolist()
        #print (x_test)
        #print (y_test)
        ### WT
        if '_WT' not in gpcr_given:
            wt = gpcr_given.split('_')[0] + '_WT'
            Xs_wt_pca = np.load(path+'/static/predictor/output/'+uniq_id+'/PCA/'+gprotein_given+'_'+gpcr_given+'.npy', allow_pickle=True)
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

        x_train_coupling, x_train_uncoupling, x_train_grey, y_train_coupling, y_train_uncoupling, y_train_grey, genes_to_consider_coupling, genes_to_consider_uncoupling, genes_to_consider_grey = extract_pca(gprotein_given, assay, pca_type)
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

@app.route('/fetchContactsSequence', methods=['GET', 'POST'])
def fetchContactsSequence():
    if request.method == 'POST':
        data = request.get_json(force=True)
        #print (data['gpcr'])
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
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

        flag = 0
        dic = {}
        for line in open(path+'/static/predictor/output/'+uniq_id+'/temp_hmm_file.txt', 'r'):
            if line[:2] == '>>':
                flag = 0
                gpcr_found = line.split('>>')[1].replace(' ', '').replace('\n', '')
                if gpcr_found == gpcr_given:
                    flag = 1
            elif flag == 1:
                if len(line.split()) > 0:
                    if line.split()[0] == '7tm_1':
                        startDomain = line.split()[1]
                        domain = line.split()[2]
                        endDomain = line.split()[3]
                    elif gpcr_found in line.split()[0]:
                        startSequence = line.split()[1]
                        sequence = line.split()[2]
                        endSequence = line.split()[3]

                        countDomain = int(startDomain)
                        countSequence = int(startSequence)
                        for seq, dom in zip(sequence, domain):
                            if seq != '-' and dom != '.':
                                dic[countDomain] = countSequence
                                countDomain += 1
                                countSequence += 1
                            elif seq != '-':
                                countSequence += 1
                            else:
                                countDomain += 1

        #print (dic)

        '''
        BW_to_PFAM = {}
        for count, line in enumerate(open('data/7tm_1_2020_BW_matches.tsv', 'r')):
            BW_to_PFAM[line.split('\t')[1]] = count+1

        #print (BW_to_PFAM)
        seq_positions = []
        for bw_position in positions:
            #print (bw_position)
            if bw_position in BW_to_PFAM:
                pfam_position = BW_to_PFAM[bw_position]
                if pfam_position in dic:
                    #print (dic[pfam_position])
                    seq_positions.append(int(dic[pfam_position]))
        '''
        gpcr = gpcr_given.split('_')[0]
        dic_BW2SEQ = {}
        for line in open(path+'/data/GPCRDB_numbering.tsv', 'r'):
            if gpcr == line.split('\t')[0] or gpcr == line.split('\t')[1]:
                BW = line.split('\t')[-2]
                SEQ = line.split('\t')[3]
                if BW != '-' and SEQ != '-':
                    dic_BW2SEQ[BW] = SEQ

        #print (positions)
        #print (dic_BW2SEQ)
        seq_positions = []
        for bw_position in positions:
            if bw_position in dic_BW2SEQ:
                seq_positions.append(int(dic_BW2SEQ[bw_position]))

        return jsonify({'fetch_contacts': scores, 'seq_positions': seq_positions, 'sequence': fasta_sequence})
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
        mutation_pfam_position = None
        if '_WT' not in gpcr_given:
            mutation_sequence_position = int(gpcr_given.split('_')[1][1:-1])
            flag = 0
            dic = {}
            for line in open(path+'/static/predictor/output/'+uniq_id+'/temp_hmm_file.txt', 'r'):
                if line[:2] == '>>':
                    flag = 0
                    gpcr_found = line.split('>>')[1].replace(' ', '').replace('\n', '')
                    if gpcr_given == gpcr_found:
                        flag = 1
                elif flag == 1:
                    if len(line.split()) > 0:
                        if line.split()[0] == '7tm_1':
                            startDomain = line.split()[1]
                            domain = line.split()[2]
                            endDomain = line.split()[3]
                        elif gpcr_found in line.split()[0]:
                            startSequence = line.split()[1]
                            sequence = line.split()[2]
                            endSequence = line.split()[3]

                            countDomain = int(startDomain)
                            countSequence = int(startSequence)
                            for seq, dom in zip(sequence, domain):
                                if seq != '-' and dom != '.':
                                    dic[countSequence] = countDomain
                                    countDomain += 1
                                    countSequence += 1
                                elif seq != '-':
                                    countSequence += 1
                                else:
                                    countDomain += 1

            #print (dic)
            if mutation_sequence_position in dic:
                mutation_pfam_position = dic[mutation_sequence_position]
                #print (mutation_pfam_position)

        BW_to_PFAM = {}
        for count, line in enumerate(open(path+'/data/7tm_1_2020_BW_matches.tsv', 'r')):
            BW_to_PFAM[line.split('\t')[1]] = count+1
        PFAM_to_PDB = {}
        for count, line in enumerate(open(path+'/data/hmmsearchPDB/'+pdbID+'.txt', 'r')):
            PFAM_to_PDB[int(line.split('\t')[3].replace('\n', ''))] = line.split('\t')[2]

        #print (BW_to_PFAM)
        #print (PFAM_to_PDB)
        #print (num_contacts)
        modified_positions = []
        modified_num_contacts = []
        for num, position in enumerate(positions.split(',')):
            if position in BW_to_PFAM:
                #print (BW_to_PFAM[position], end=' ')
                if BW_to_PFAM[position] in PFAM_to_PDB:
                    pdbPosition = PFAM_to_PDB[BW_to_PFAM[position]]
                    modified_positions.append(pdbPosition)
                    modified_num_contacts.append(str(num_contacts.split(',')[num]))
        #print (modified_num_contacts)

        modified_pair_positions = []
        #print (pair_positions)
        for position in pair_positions.split(','):
            pos1 = position.split(':')[0]
            pos2 = position.split(':')[1]
            score = position.split(':')[2]
            if pos1 in BW_to_PFAM and pos2 in BW_to_PFAM:
                #print (BW_to_PFAM[position], end=' ')
                if BW_to_PFAM[pos1] in PFAM_to_PDB and BW_to_PFAM[pos2] in PFAM_to_PDB:
                    pdbPosition1 = PFAM_to_PDB[BW_to_PFAM[pos1]]
                    pdbPosition2 = PFAM_to_PDB[BW_to_PFAM[pos2]]
                    modified_pair_positions.append(pdbPosition1+':'+pdbPosition2+':'+score)

        #print (modified_positions)
        mutation_position = '-'
        if mutation_pfam_position != None:
            for pfamPosition in PFAM_to_PDB:
                if pfamPosition == mutation_pfam_position:
                    #print (PFAM_to_PDB[pfamPosition])
                    mutation_position = PFAM_to_PDB[pfamPosition]
                    break

        return jsonify({'modified_positions': '_'.join(modified_positions),
                        'modified_num_contacts': '_'.join(modified_num_contacts),
                        'modified_pair_positions': '_'.join(modified_pair_positions),
                        'mutation_position': mutation_position
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
        ordered_pdbs = reorder_pdbs(gpcr_given, uniq_id) ## return list of reordered PDB IDs based on GPCR
        return jsonify({'try': positions.tolist(),
                        'ordered_pdbs': ordered_pdbs,
                        'positions': ','.join(positions.tolist()),
                        'num_contacts': ','.join(num_contacts),
                        'pair_positions': ','.join(pair_positions)})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

## Function to return list of PDB IDs based on GPCR using BLAST
def reorder_pdbs(gpcr, uniq_id):
    path_to_fasta = path+"/static/predictor/output/"+uniq_id+"/input.fasta"
    path_to_output = path+"/static/predictor/output/"+uniq_id+"/"

    os.system('blastp -query '+path_to_fasta+' -db '+ path + '/data/fasta/blastdb/all_pdbs -out '+path_to_output+'/blastp_output.txt')

    chain_info = {}
    for line in open(path+'/data/pdblist.txt', 'r'):
        pdbid = line.split(' ')[0]
        gpcr_chain = line.split(' ')[1]
        gprotein_chain = line.split(' ')[2].replace('\n', '')
        chain_info[pdbid] = {}
        chain_info[pdbid]['gpcr_chain'] = gpcr_chain
        chain_info[pdbid]['gprotein_chain'] = gprotein_chain

    dic = {}
    for line in open(path_to_output+'/blastp_output.txt', 'r'):
        if 'Query=' in line:
            name = line.split('Query=')[1].replace('\n', '').replace(' ', '')
            dic[name] = []
        elif line[0] == '>':
            pdbid = line.split('>')[1].split('|')[0].split('_')[0].lower()
            row = []
            #print (chain_info[pdbid])
            row.append(pdbid)
            row.append(chain_info[pdbid]['gpcr_chain'])
            row.append(chain_info[pdbid]['gprotein_chain'])
            #dic[name].append(row)
            dic[name].append(pdbid+'_'+chain_info[pdbid]['gpcr_chain']+'_'+chain_info[pdbid]['gprotein_chain'])

    return(dic[gpcr])
    #return None

## Route to output page
@app.route('/input', methods=['GET', 'POST'])
def input():
    if request.method == 'POST':
        input = request.form['input']
        input_file = None ## Upload FASTA file
        ## Run the predictor
        uniq_id = precogx.main(input, input_file, 'all', app.root_path)
        #uniq_id = 'OXDUB'
        return redirect('/output/'+uniq_id)
    else:
        return ("<html><h3>It was a GET request</h3></html>")

@app.route('/output/<uniq_id>', methods=['GET', 'POST'])
def output(uniq_id):
    if request.method == 'GET' or request.method == 'POST':
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
        gpcr_list = []
        for key1 in d:
            for num, key2 in enumerate(d[key1]):
                gpcr = key2[0]
                if key2[1] == 'WT' or key2[1] not in ['IUPHAR', 'LogRAi', 'Emax']:
                    variant = '_' + key2[1]
                else:
                    variant = '_WT'
                if num == 0:
                    first_entry = gpcr+variant
                    #print (key1, d[key1])
                gpcr_list.append(gpcr+variant)
                #break

        #print (first_entry)
        #path_to_json_output = "/static/predictor/output/"+uniq_id+"/out.json"
        #path_to_fasta = "/static/predictor/output/"+uniq_id+"/input.fasta"
        return render_template('result.html',
                                path_to_json_output=json.dumps(path_to_json_output),
                                path_to_fasta=json.dumps(path_to_fasta),
                                first_entry=json.dumps(first_entry),
                                gpcr_list=json.dumps(gpcr_list),
                                uniq_id=json.dumps(uniq_id))
    else:
        return ("<html><h3>It was a GET request</h3></html>")

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
