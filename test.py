import os, sys, gzip

os.chdir('data/')
os.system('rm -rf pdb_chain_pfam.tsv.gz')
os.system('wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz')

dic = {}
for line in gzip.open('pdb_chain_pfam.tsv.gz', 'rt'):
    if line[0] != '#' and line.split()[0] != 'PDB':
        #print (line)
        pdbid = line.split('\t')[0]
        pfam = line.split('\t')[3]
        if pfam == 'PF00001' or pfam == 'PF00503':
            if pdbid not in dic:
                dic[pdbid] = {}
                dic[pdbid]['GPCR'] = '-'
                dic[pdbid]['GPROT'] = '-'
            if pfam == 'PF00001':
                dic[pdbid]['GPCR'] = line.split('\t')[1]
            elif pfam == 'PF00503':
                dic[pdbid]['GPROT'] = line.split('\t')[1]

l = ''
for pdbid in dic:
    if dic[pdbid]['GPCR'] != '-' and dic[pdbid]['GPROT'] != '-':
        l += pdbid + ' ' + dic[pdbid]['GPCR'] + ' ' + dic[pdbid]['GPROT'] + '\n'
        os.system('wget https://www.rcsb.org/fasta/chain/'+pdbid.upper()+'.'+dic[pdbid]['GPCR'] + ' -O ../static/fasta/'+ pdbid.upper()+'.'+dic[pdbid]['GPCR']+'.fasta')

open('../static/pdblist.txt', 'w').write(l)
