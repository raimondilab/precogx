import os
import argparse, requests, urllib

for files in os.listdir('.'):
    '''
    if '_' in files and files.endswith('.pdb') == True:
        print (files)
        os.system('mv '+files+' '+files.replace('_',':'))
    '''
    if ':' in files and files.endswith('.pdb') == True:
        print (files)
        os.system('mv '+files+' AF:'+files.replace('_',':'))
    '''
    if '-' in files and files.endswith('.pdb') == True:
        print (files)
        os.system('mv '+files+' '+files.split('.')[0]+'.pdb')
    '''
    '''
    if '-' in files and files.endswith('.pdb') == True:
        url = 'https://www.uniprot.org/uploadlists/'
        #print (name)
        params = {
        'to': 'GENENAME',
        'from': 'ACC',
        'format': 'tab',
        'taxon': '9606',
        'query': files.split('-')[0]
        }

        data = urllib.parse.urlencode(params)
        data = data.encode('utf-8')
        req = urllib.request.Request(url, data)
        with urllib.request.urlopen(req) as f:
           response = f.read()
        acc = None
        seq = None
        for line in response.decode('utf-8').split('\n'):
            print (line)
            continue
            if 'From' not in line:
                if len(line.split())>0:
                    gene = str(line.split('\t')[1].replace('\n', ''))
                    break
        break
    '''
