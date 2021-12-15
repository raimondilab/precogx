## to run on workflow trigger

import os, sys, gzip
from Bio.PDB import *
from Bio.SeqUtils import seq1
from Bio import SearchIO
from Bio.Blast import NCBIXML

## Remove the old version of PFAM clans
os.system('wget http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.clans.tsv.gz -O data/Pfam/Pfam-A.clans.tsv.gz')
gpcrPFAM = []

for line in gzip.open('data/Pfam/Pfam-A.clans.tsv.gz', 'rt'):
    if line.split('\t')[2] == 'GPCR_A':
        gpcrPFAM.append(line.split('\t')[0])

#print (gpcrPFAM)
#sys.exit()
## Remove the old version of SIFT data
os.system('rm -rf data/SIFTS/pdb_chain_pfam.tsv.gz')

## Download recent version of SIFT data
os.system('wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/pdb_chain_pfam.tsv.gz -O data/SIFTS/pdb_chain_pfam.tsv.gz')

## Extract PDB ID, GPCR and G-prot chain information from SIFT data
dic = {}
for line in gzip.open('data/SIFTS/pdb_chain_pfam.tsv.gz', 'rt'):
    if line[0] != '#' and line.split()[0] != 'PDB':
        #print (line)
        pdbID = line.split('\t')[0]
        pfam = line.split('\t')[3]
        if pfam in ['PF00503', 'PF00339', 'PF02752'] or pfam in gpcrPFAM:
            if pdbID not in dic:
                dic[pdbID] = {}
                dic[pdbID]['GPCR'] = '-'
                dic[pdbID]['GPROT'] = '-'
                dic[pdbID]['BARR'] = '-'
            if pfam in gpcrPFAM:
                dic[pdbID]['GPCR'] = line.split('\t')[1]
            elif pfam == 'PF00503':
                dic[pdbID]['GPROT'] = line.split('\t')[1]
            elif pfam in ['PF00339', 'PF02752']:
                dic[pdbID]['BARR'] = line.split('\t')[1]
                #print (pfam)

def makeMapFASTA(pdbID, dic):
    print (pdbID)
    SEQ2PDB = {};
    fasta = ''
    x = ''
    pdbl = PDBList()
    try:
        pdbl.retrieve_pdb_file(pdbID.upper(), pdir='data/PDB/pdir', obsolete=False)
        parser = MMCIFParser()
        structure = parser.get_structure(x, 'data/PDB/pdir/'+pdbID+'.cif')
    except:
        pdbl.retrieve_pdb_file(pdbID.upper(), pdir='data/PDB/pdir', obsolete=False, file_format="pdb")
        parser = PDBParser()
        structure = parser.get_structure(x, 'data/PDB/pdir/pdb'+pdbID+'.ent')
    print ('pass')
    for model in structure:
        for chain in model:
            if chain.id == dic[pdbID]['GPCR']:
                fasta = '>'+pdbID+'|'+dic[pdbID]['GPCR'] + '\n'
                num = 1
                for residue in chain:
                    for atom in residue:
                        if atom.id == 'CA':
                            _, position, _ = residue.id
                            aa = seq1(residue.resname)
                            fasta += aa
                            SEQ2PDB[num] = int(position)
                            num += 1
                            #print(chain.id, aa, position)
                break
    #print (fasta)
    open('data/PDB/fasta/'+pdbID+'.fasta', 'w').write(fasta)
    os.system('blastp -query data/PDB/fasta/'+pdbID+'.fasta'+' -outfmt 5 -out ' + 'data/PDB/fasta/'+pdbID+'.txt -db data/GPCRDB/blastdb/GPCRDB')
    #print (map)

    handle = open('data/PDB/fasta/'+pdbID+'.txt', 'r')
    blast_records = NCBIXML.parse(handle)
    #print (blast_records)

    GPCRDB2SEQ = {}
    for blast_record in blast_records:
        #print (blast_record.query)
        query = blast_record.query
        chainID = query.split('|')[1].replace('\n', '')

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                bestHIT = alignment.title.split(' ')[1]
                for num, (q, s) in enumerate(zip(hsp.query, hsp.sbjct)):
                    if q!='-' and s!='-':
                        GPCRDB2SEQ[num + hsp.sbjct_start] = num + hsp.query_start
                break
            break

    #print (SEQ2PDB)
    #print (GPCRDB2SEQ)
    l = ''
    for s in GPCRDB2SEQ:
        q = GPCRDB2SEQ[s]
        if q in SEQ2PDB:
            pdbPosition = SEQ2PDB[q]
            l += str(pdbID) + '\t' + chainID + '\t' + str(s) + '\t' + str(q) + '\t' + bestHIT + '\n'
            #print (l)
        open('data/PDB/GPCRDB/'+pdbID+'.txt', 'w').write(l)

    return fasta

## Save PDB ID and chain information as well as
## fetch FASTA files of PDB IDs
l = ''
pdblist = []; allPDB = ''
for pdbID in dic:
    if (dic[pdbID]['GPCR'] != '-' and dic[pdbID]['GPROT'] != '-') or (dic[pdbID]['GPCR'] != '-' and dic[pdbID]['BARR'] != '-'):
        l += pdbID + ' ' + dic[pdbID]['GPCR'] + ' ' + dic[pdbID]['GPROT'] + ' ' + dic[pdbID]['BARR'] + '\n'
        #os.system('wget https://www.rcsb.org/fasta/entry/'+pdbID.lower()+'/download' + ' -O ../data/fasta/'+ pdbID.lower()+'.fasta')
        allPDB += makeMapFASTA(pdbID, dic) + '\n'
        pdblist.append(pdbID)
        #break

#print (l)
open('data/PDB/pdblist.txt', 'w').write(l)

## Concatenate all FASTA files into one and make a blastdb in blastdb/
os.system("rm -rf data/PDB/allPDB.fasta")
open ('data/PDB/allPDB.fasta', 'w').write(allPDB)
if os.path.isfile('data/PDB/blastdb/') == False:
    os.system("mkdir blastdb/")
os.system("rm -rf data/PDB/blastdb/allPDB*")
os.system("makeblastdb -in data/PDB/allPDB.fasta -dbtype 'prot' -out data/PDB/blastdb/allPDB")
print ('complete')
