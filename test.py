import os, sys, gzip

os.chdir('data/')
os.system('wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz')

l = ''
for line in gzip.open('pdb_chain_pfam.tsv.gz', 'rt'):
    if line[0] != '#' and line.split()[0] != 'PDB':
        #print (line)
        pdbid = line.split('\t')[0]
        pfam = line.split('\t')[3]
        if pfam == 'PF00001':
            l += pdbid + '\n'

open('../static/help.txt', 'w').write(l)
