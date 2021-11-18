from flask import Flask, render_template, request, jsonify
import os, json

app = Flask(__name__)

## Route to home page
@app.route('/', methods=['GET', 'POST'])
@app.route('/home', methods=['GET', 'POST'])
def home():
    return render_template('home.html')

@app.route('/order_pdbs', methods=['GET', 'POST'])
def order_pdbs():
    if request.method == 'POST':
        data = request.get_json(force=True)
        print (data['gpcr'])
        ordered_pdbs = reorder_pdbs(data['gpcr'])
        return jsonify({'ordered_pdbs': ordered_pdbs})
    else:
        return ("<html><h3>It was a GET request</h3></html>")


def reorder_pdbs(gpcr):
    path_to_fasta = "static/OL820/temp/new_fasta_file.txt"
    path_to_output = "static/OL820/temp/"

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
            print (chain_info[pdbid])
            row.append(pdbid)
            row.append(chain_info[pdbid]['gpcr_chain'])
            row.append(chain_info[pdbid]['gprotein_chain'])
            #dic[name].append(row)
            dic[name].append(pdbid+'_'+chain_info[pdbid]['gpcr_chain']+'_'+chain_info[pdbid]['gprotein_chain'])

    #print(dic[gpcr])

    return(dic[gpcr])
    #return None

## Route to output page
@app.route('/output', methods=['GET', 'POST'])
def output():
    if request.method == 'POST':
        #Path to sample output
        path_to_json_output = "static/OL820/out.json"
        path_to_fasta = "static/OL820/temp/new_fasta_file.txt"
        ## extract first entry
        with open(path_to_json_output) as f:
            d = json.load(f)
        first_entry = 'Hello'
        for key1 in d:
            for key2 in d[key1]:
                first_entry = key2[0]
                break

        return render_template('embedded.html', path_to_json_output=json.dumps(path_to_json_output), path_to_fasta=json.dumps(path_to_fasta), first_entry=json.dumps(first_entry))
    return ("<html><h3>It was a GET request</h3></html>")

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')
