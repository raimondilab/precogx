## to run on workflow trigger

import os, sys, gzip

os.chdir('data/')

## Remove the old version of SIFT data
os.system('rm -rf pdb_chain_pfam.tsv.gz')

## Download recent version of SIFT data
os.system('wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz')

## Extract PDB ID, GPCR and G-prot chain information from SIFT data
dic = {}
for line in gzip.open('pdb_chain_pfam.tsv.gz', 'rt'):
    if line[0] != '#' and line.split()[0] != 'PDB':
        #print (line)
        pdbid = line.split('\t')[0]
        pfam = line.split('\t')[3]
        if pfam in ['PF00001', 'PF00503', 'PF00339', 'PF02752']:
            if pdbid not in dic:
                dic[pdbid] = {}
                dic[pdbid]['GPCR'] = '-'
                dic[pdbid]['GPROT'] = '-'
                dic[pdbid]['BARR'] = '-'
            if pfam == 'PF00001':
                dic[pdbid]['GPCR'] = line.split('\t')[1]
            elif pfam == 'PF00503':
                dic[pdbid]['GPROT'] = line.split('\t')[1]
            elif pfam in ['PF00339', 'PF02752']:
                dic[pdbid]['BARR'] = line.split('\t')[1]
                #print (pfam)

## Save PDB ID and chain information as well as
## fetch FASTA files of PDB IDs
l = ''
for pdbid in dic:
    if (dic[pdbid]['GPCR'] != '-' and dic[pdbid]['GPROT'] != '-') or (dic[pdbid]['GPCR'] != '-' and dic[pdbid]['BARR'] != '-'):
        l += pdbid + ' ' + dic[pdbid]['GPCR'] + ' ' + dic[pdbid]['GPROT'] + ' ' + dic[pdbid]['BARR'] + '\n'
        os.system('wget https://www.rcsb.org/fasta/entry/'+pdbid.lower()+'/download' + ' -O ../data/fasta/'+ pdbid.lower()+'.fasta')

#print (l)
#sys.exit()
open('../data/pdblist.txt', 'w').write(l)

## Concatenate all FASTA files into one and make a blastdb in blastdb/
os.chdir('../data/fasta/')
# Delete existing all_pdbs.fasta
os.system("rm -rf all_pdbs.fasta")
l = ''
for files in os.listdir('.'):
    if files.endswith('.fasta'):
        for line in open(files, 'r'):
            l += line

open ('all_pdbs.fasta', 'w').write(l)
if os.path.isfile('blastdb/') == False:
    os.system("mkdir blastdb/")
os.system("rm -rf blastdb/all_pdbs*")
os.system("makeblastdb -in all_pdbs.fasta -dbtype 'prot' -out blastdb/all_pdbs")
