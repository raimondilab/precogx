#!/usr/bin/env python3
from Bio import SearchIO
from Bio.Blast import NCBIXML
import os, sys, gzip

class mapping:
    def __init__(self, pdbID):
        self.pdbID = ''
        self.map = {}
        self.bestHIT = {}

PDB2CHAIN = {}
for line in gzip.open('data/pdb_chain_pfam.tsv.gz', 'rt'):
    if line[0] != '#':
        pdbID = line.split('\t')[0]
        chainID = line.split('\t')[1]
        pfamID = line.split('\t')[3]
        if pfamID == 'PF00001':
            PDB2CHAIN[pdbID] = chainID

#os.system('blastp -query data/fasta/all_pdbs.fasta -outfmt 5 -out ' + 'data/PDB/GPCRDBblast.txt -db data/GPCRDB/blastdb/GPCRDB')

handle = open('data/PDB/GPCRDBblast.txt', 'r')
blast_records = NCBIXML.parse(handle)
#print (blast_records)

dic = {}
for blast_record in blast_records:
    #print (blast_record.query)
    query = blast_record.query
    if '[auth' in query:
        chainID = query.split('|')[1].split('[auth')[1].split(']')[0].replace(' ', '')
    else:
        chainID = query.split('|')[1].split('Chain')[1].replace(' ', '')

    pdbID = query.split('|')[0].split('_')[0].lower()

    if chainID == PDB2CHAIN[pdbID]:
        if pdbID not in dic:
            dic[pdbID] = mapping(pdbID)
        dic[pdbID].map[chainID] = {}

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                bestHIT = alignment.title.split(' ')[1]
                dic[pdbID].bestHIT[chainID] = bestHIT
                for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                    if q!='-' and s!='-':
                        dic[pdbID].map[chainID][num + hsp.sbjct_start] = num + hsp.query_start
                break
            break

for pdbID in dic:
    l = ''
    for chainID in dic[pdbID].map:
        for s in dic[pdbID].map[chainID]:
            q = dic[pdbID].map[chainID][s]
            bestHIT = dic[pdbID].bestHIT[chainID]
            l += str(pdbID) + '\t' + chainID + '\t' + str(s) + '\t' + str(q) + '\t' + bestHIT + '\n'
    #print (l)
    open('data/PDB/GPCRDB/'+pdbID+'.txt', 'w').write(l)
