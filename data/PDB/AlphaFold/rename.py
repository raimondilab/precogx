import os
import argparse, requests, urllib

dic = {}
for line in open('/home/gurdeep/projects/DB/uniprot/uniprot_sprot.fasta', 'r'):
    if line[0] == '>':
        acc = line.split('|')[1]
        if 'GN=' in line:
            gene = line.split('GN=')[1].split()[0]
            dic[acc] = gene

dic['P09471-2']='GoB'

l = []
for files in os.listdir('.'):
    '''
    if '_' in files and files.endswith('.pdb') == True:
        print (files)
        os.system('mv '+files+' '+files.replace('_',':'))
    '''
    '''
    if ':' in files and files.endswith('.pdb') == True:
        print (files)
        os.system('mv '+files+' AF:'+files.replace('_',':'))
    '''
    '''
    if ':' in files and files.endswith('.pdb') == True:
        acc1 = files.split(':')[1]
        acc2 = files.split(':')[2].split('.')[0]
        if acc1 not in dic:
            print ('not found', acc1)
        if acc2 not in dic:
            print ('not found', acc2)
        else:
            l.append(dic[acc2])

        newName = 'AF:'+dic[acc1]+':'+dic[acc2]+'.pdb'
        os.system('mv '+files+' '+newName)
    '''
    if ':' in files and files.endswith('.pdb') == True:
        gene1 = files.split(':')[1]
        gene2 = files.split(':')[2].split('.')[0]
        if gene2 == 'GoA':
            newName = 'AF:'+gene1+':GoB.pdb'
            os.system('mv '+files+' '+newName)
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

print (list(set(l)))
