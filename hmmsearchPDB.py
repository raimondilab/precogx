#!/usr/bin/env python3

import os, sys

#os.system('hmmsearch -o data/fasta/all_pdbs_hmmsearch.txt data/7tm_1.hmm data/fasta/all_pdbs.fasta')
os.system('hmmsearch -o data/fasta/all_pdbs_hmmsearch.txt data/SCOP_7TM_348.hmm data/fasta/all_pdbs.fasta')

class mapping:
    def __init__(self, pdbID):
        self.pdbID = ''
        self.map = {}

flag = 0
dic = {}
for line in open('data/fasta/all_pdbs_hmmsearch.txt', 'r'):
    if line[:2] == '>>':
        if '[auth' in line:
            chainID = line.split('|')[1].split('[auth')[1].split(']')[0].replace(' ', '')
        else:
            chainID = line.split('|')[1].split('Chain')[1].replace(' ', '')
        pdbID = line.split('|')[0].split('>>')[1].replace(' ', '').split('_')[0].lower()
        entityID = line.split('|')[0].split('>>')[1].replace(' ', '')
        print (chainID, pdbID)
        if pdbID not in dic:
            dic[pdbID] = mapping(pdbID)
        dic[pdbID].map[chainID] = {}
        flag = 1
    elif flag == 1:
        if len(line.split()) > 0:
            if line.split()[0] == '0037432':
                startDomain = line.split()[1]
                domain = line.split()[2]
                endDomain = line.split()[3]
            elif entityID in line.split()[0]:
                startSequence = line.split()[1]
                sequence = line.split()[2]
                endSequence = line.split()[3]

                countDomain = int(startDomain)
                countSequence = int(startSequence)
                for seq, dom in zip(sequence, domain):
                    if seq != '-' and dom != '.':
                        dic[pdbID].map[chainID][countDomain] = countSequence
                        countDomain += 1
                        countSequence += 1
                    elif seq != '-':
                        countSequence += 1
                    else:
                        countDomain += 1

for pdbID in dic:
    l = ''
    for chainID in dic[pdbID].map:
        for countDomain in dic[pdbID].map[chainID]:
            countSequence = dic[pdbID].map[chainID][countDomain]
            l += str(pdbID) + '\t' + chainID + '\t' + str(countSequence) + '\t' + str(countDomain) + '\n'
    print (l)
    open('data/hmmsearchPDB/'+pdbID+'.txt', 'w').write(l)
