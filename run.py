from flask import Flask, render_template, request, jsonify, redirect, url_for
import os, sys, json
import numpy as np
sys.path.insert(1, 'static/predictor/')
import precogx

app = Flask(__name__)

## Route to home page
@app.route('/', methods=['GET', 'POST'])
@app.route('/home', methods=['GET', 'POST'])
def home():
    return render_template('home.html')

def extract_contacts(gprotein_given, cutoff):
    dic = {}; positions = []; pair_positions = []
    for line in open('static/Gprot_contacts/specific_position_score.txt', 'r'):
        if line.split(' ')[0] != 'TM_pos1':
            pos1 = line.split(' ')[0]
            pos2 = line.split(' ')[1]
            gprotein_found = line.split(' ')[2]
            score = float(line.replace('\n','').split(' ')[3])
            if score >= cutoff and gprotein_given == gprotein_found:
                #print (score)
                if pos1 not in dic:
                    dic[pos1] = {}
                dic[pos1][pos2] = score

                if pos2 not in dic:
                    dic[pos2] = {}
                dic[pos2][pos1] = score

                positions.append(pos1)
                positions.append(pos2)
                pair_positions.append(pos1+'-'+pos2)

    positions = list(set(positions))
    positions = np.array(np.sort(positions))
    data = []
    for pos1 in positions:
        row = []
        for pos2 in positions:
            if pos2 in dic[pos1]:
                row.append(dic[pos1][pos2])
            else:
                #row.append(0)
                row.append(None)
        data.append(row)
    #print (cutoff, gprotein_given)
    #print ('positions', positions)
    return data, positions, pair_positions

#
@app.route('/fetchContactsHeatmap', methods=['GET', 'POST'])
def fetchContactsHeatmap():
    if request.method == 'POST':
        data = request.get_json(force=True)
        #print (data['gpcr'])
        gprotein_given = data['gprotein']
        cutoff = float(data['cutoff'])
        scores, positions, pair_positions = extract_contacts(gprotein_given, cutoff)
        for i in range(0, len(positions)):
            positions[i] = str(positions[i]).replace('.', 'x')
        return jsonify({'fetch_contacts': scores, 'positions': positions.tolist()})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

def extract_pca(gprotein, label):
    Xs_train_pca = np.load('static/Xs_train_pca/'+gprotein+'.npy', allow_pickle=True)
    Xs_train_pca_coupling, Xs_train_pca_uncoupling = filter_gpcr_list(Xs_train_pca, label, gprotein)
    x_train_coupling = Xs_train_pca_coupling[:,0].tolist()
    x_train_uncoupling = Xs_train_pca_uncoupling[:,0].tolist()
    y_train_coupling = Xs_train_pca_coupling[:,1].tolist()
    y_train_uncoupling = Xs_train_pca_uncoupling[:,1].tolist()
    return x_train_coupling, x_train_uncoupling, y_train_coupling, y_train_uncoupling

def filter_gpcr_list(X, label, gprotein):
    genes_to_consider_coupling = []
    genes_to_consider_uncoupling = []
    #print (label)
    #label = 'ebBRET'
    if label == 'Shedding':
        num = -1
        for line in open('static/predictor/data_precog/LogRAi_values_final.tsv', 'r'):
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

    elif label == 'ebBRET':
        num = -1
        for line in open('static/predictor/data_precog2/emax.tsv', 'r'):
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
        #print (genes_to_consider_coupling)
        #print (genes_to_consider_uncoupling)

    gpcr_list = []
    for line in open('static/predictor/gpcr_list.txt', 'r'):
        gpcr_list.append(line.replace('\n', ''))

    X_pos = []
    X_neg = []
    for gene, row in zip(gpcr_list, X):
        if gene in genes_to_consider_coupling:
            X_pos.append(row)
        elif gene in genes_to_consider_uncoupling:
            X_neg.append(row)

    return (np.array(X_pos), np.array(X_neg))

@app.route('/fetchPCA', methods=['GET', 'POST'])
def fetchPCA():
    if request.method == 'POST':
        data = request.get_json(force=True)
        label = data['label']
        gprotein_given = data['gprotein']
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        Xs_test_pca = np.load('static/predictor/output/'+uniq_id+'/'+gprotein_given+'.npy', allow_pickle=True)
        x_test = Xs_test_pca[:,0].tolist()
        y_test = Xs_test_pca[:,1].tolist()
        x_train_coupling, x_train_uncoupling, y_train_coupling, y_train_uncoupling = extract_pca(gprotein_given, label)
        #print (x_train, y_train, x_test, y_test)
        minimum = min(x_train_coupling + x_train_uncoupling)
        maximum = max(x_train_coupling + x_train_uncoupling)
        return jsonify({'x_train_coupling': x_train_coupling,
                        'x_train_uncoupling': x_train_uncoupling,
                        'y_train_coupling': y_train_coupling,
                        'y_train_uncoupling': y_train_uncoupling,
                        'x_test': x_test,
                        'y_test': y_test,
                        'min': minimum,
                        'max': maximum})
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
        scores, positions, pair_positions = extract_contacts(gprotein_given, cutoff)
        positions = positions

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
        for line in open('static/predictor/output/'+uniq_id+'/temp_hmm_file.txt', 'r'):
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
        gpcr_given = data['gpcr']
        uniq_id = data['uniq_id']
        #print (pdbID)
        #print (pair_positions)
        mutation_pfam_position = None
        if '_WT' not in gpcr_given:
            mutation_sequence_position = int(gpcr_given.split('_')[1][1:-1])
            flag = 0
            dic = {}
            for line in open('static/predictor/output/'+uniq_id+'/temp_hmm_file.txt', 'r'):
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
        for count, line in enumerate(open('data/7tm_1_2020_BW_matches.tsv', 'r')):
            BW_to_PFAM[line.split('\t')[1]] = count+1
        PFAM_to_PDB = {}
        for count, line in enumerate(open('data/hmmsearchPDB/'+pdbID+'.txt', 'r')):
            PFAM_to_PDB[int(line.split('\t')[3].replace('\n', ''))] = line.split('\t')[2]

        #print (BW_to_PFAM)
        #print (PFAM_to_PDB)
        modified_positions = []
        for position in positions.split(','):
            if position in BW_to_PFAM:
                #print (BW_to_PFAM[position], end=' ')
                if BW_to_PFAM[position] in PFAM_to_PDB:
                    pdbPosition = PFAM_to_PDB[BW_to_PFAM[position]]
                    modified_positions.append(pdbPosition)

        modified_pair_positions = []
        #print (pair_positions)
        for position in pair_positions.split(','):
            pos1 = position.split('-')[0]
            pos2 = position.split('-')[1]
            if pos1 in BW_to_PFAM and pos2 in BW_to_PFAM:
                #print (BW_to_PFAM[position], end=' ')
                if BW_to_PFAM[pos1] in PFAM_to_PDB and BW_to_PFAM[pos2] in PFAM_to_PDB:
                    pdbPosition1 = PFAM_to_PDB[BW_to_PFAM[pos1]]
                    pdbPosition2 = PFAM_to_PDB[BW_to_PFAM[pos2]]
                    modified_pair_positions.append(pdbPosition1+'-'+pdbPosition2)

        #print (modified_positions)
        mutation_position = '-'
        if mutation_pfam_position != None:
            for pfamPosition in PFAM_to_PDB:
                if pfamPosition == mutation_pfam_position:
                    #print (PFAM_to_PDB[pfamPosition])
                    mutation_position = PFAM_to_PDB[pfamPosition]
                    break

        return jsonify({'modified_positions': '_'.join(modified_positions),
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
        scores, positions, pair_positions = extract_contacts(gprotein_given, cutoff)
        ordered_pdbs = reorder_pdbs(gpcr_given, uniq_id) ## return list of reordered PDB IDs based on GPCR
        return jsonify({'ordered_pdbs': ordered_pdbs, 'positions': ','.join(positions.tolist()), 'pair_positions': ','.join(pair_positions)})
    else:
        return ("<html><h3>It was a GET request</h3></html>")

## Function to return list of PDB IDs based on GPCR using BLAST
def reorder_pdbs(gpcr, uniq_id):
    path_to_fasta = "static/predictor/output/"+uniq_id+"/input.fasta"
    path_to_output = "static/predictor/output/"+uniq_id+"/"

    os.system('blastp -query '+path_to_fasta+' -db data/fasta/blastdb/all_pdbs -out '+path_to_output+'/blastp_output.txt')

    chain_info = {}
    for line in open('data/pdblist.txt', 'r'):
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
        uniq_id = precogx.main(input, input_file, 'all')
        #uniq_id = 'OXDUB'
        return redirect('/output/'+uniq_id)
    else:
        return ("<html><h3>It was a GET request</h3></html>")

@app.route('/output/<uniq_id>', methods=['GET', 'POST'])
def output(uniq_id):
    if request.method == 'GET' or request.method == 'POST':
        #print (os.getcwd())
        path_to_json_output = "/static/predictor/output/"+uniq_id+"/out.json"
        path_to_fasta = os.getcwd()+ "/static/predictor/output/"+uniq_id+"/input.fasta"

        ## extract first entry
        with open(os.getcwd() + path_to_json_output) as f:
            d = json.load(f)
        first_entry = 'Hello'
        gpcr_list = []
        for key1 in d:
            for num, key2 in enumerate(d[key1]):
                if num == 0:
                    first_entry = key2[0]
                    #print (key1, d[key1])
                gpcr_list.append(key2[0])
                #break

        #print (first_entry)
        #path_to_json_output = "/static/predictor/output/"+uniq_id+"/out.json"
        #path_to_fasta = "/static/predictor/output/"+uniq_id+"/input.fasta"
        return render_template('embedded.html',
                                path_to_json_output=json.dumps(path_to_json_output),
                                path_to_fasta=json.dumps(path_to_fasta),
                                first_entry=json.dumps(first_entry),
                                gpcr_list=json.dumps(gpcr_list),
                                uniq_id=json.dumps(uniq_id))
    else:
        return ("<html><h3>It was a GET request</h3></html>")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
