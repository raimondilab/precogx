## to run on workflow trigger
import os, sys, gzip
from Bio.PDB import *
from Bio.SeqUtils import seq1
from Bio.PDB.Polypeptide import is_aa  
from Bio.Blast import NCBIXML

# Ensure necessary directories exist immediately
os.makedirs('data/PDB/fasta', exist_ok=True)
os.makedirs('data/PDB/GPCRDB', exist_ok=True)
os.makedirs('data/PDB/pdir', exist_ok=True)
os.makedirs('data/PDB/blastdb', exist_ok=True)

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

## Extract ID, GPCR and G-prot chain information from AlphaFold data
for files in os.listdir('data/PDB/AlphaFold/'):
    if files.endswith('.pdb'):
        pdbID = files.split('.pdb')[0]
        dic[pdbID] = {}
        dic[pdbID]['GPCR'] = 'A'
        dic[pdbID]['GPROT'] = 'B'
        dic[pdbID]['BARR'] = '-'

def makeMapFASTA(pdbID, dic):
    print(f"--> Processing {pdbID}...")
    SEQ2PDB = {}
    fasta = ''
    
    # 1. Define paths
    pdir = 'data/PDB/pdir'
    fasta_dir = 'data/PDB/fasta'
    gpcrdb_dir = 'data/PDB/GPCRDB'
    cif_path = f'{pdir}/{pdbID}.cif'
    
    # 2. Load Structure (Robust Method)
    structure = None
    try:
        # Handle AlphaFold structures
        if pdbID.startswith('AF:'):
             parser = PDBParser(QUIET=True)
             structure = parser.get_structure('struct', f'data/PDB/AlphaFold/{pdbID}.pdb')
        else:
            # Handle Standard PDBs (CIF format preferred)
            parser = MMCIFParser(QUIET=True)
            
            # If file doesn't exist or is empty, force download via wget
            if not os.path.isfile(cif_path) or os.path.getsize(cif_path) == 0:
                print(f"   Downloading {pdbID}.cif...")
                # Redirect output to /dev/null to keep logs clean
                os.system(f'wget https://www.ebi.ac.uk/pdbe/entry-files/{pdbID}.cif -O {cif_path} >/dev/null 2>&1')
            
            # Attempt to parse
            if os.path.isfile(cif_path) and os.path.getsize(cif_path) > 0:
                structure = parser.get_structure('struct', cif_path)
            else:
                print(f"❌ Error: Could not download/find {pdbID}.cif")
                return ''

    except Exception as e:
        print(f"❌ CRITICAL ERROR loading structure {pdbID}: {e}")
        return ''

    # 3. Extract Sequence from Chain
    target_chain = dic[pdbID]['GPCR']
    found = False
    
    # Iterate through models and chains to find the target GPCR chain
    for model in structure:
        for chain in model:
            if chain.id == target_chain:
                fasta = f'>{pdbID}|{target_chain}\n'
                num = 1
                for residue in chain:
                    # check if it is a valid Amino Acid (standard or modified)
                    if is_aa(residue, standard=False):
                        for atom in residue:
                            if atom.id == 'CA': # Use Alpha Carbon for numbering
                                _, position, _ = residue.id
                                try:
                                    # Try to convert 3-letter code to 1-letter
                                    aa = seq1(residue.resname)
                                except:
                                    # Fallback for unknown residues to prevent crash
                                    aa = 'X' 
                                fasta += aa
                                SEQ2PDB[num] = int(position)
                                num += 1
                found = True
                break
        if found: break

    # Validation: If chain not found or sequence is too short (empty)
    if not found or len(fasta.split('\n')[-1]) < 10:
        print(f"⚠️ WARNING: GPCR Chain '{target_chain}' empty or not found in {pdbID}")
        return ''

    # 4. Save FASTA and Run BLAST
    fasta_file = f'{fasta_dir}/{pdbID}.fasta'
    blast_out = f'{fasta_dir}/{pdbID}.txt'
    
    with open(fasta_file, 'w') as f:
        f.write(fasta)
    
    # Execute BLASTP
    # Ensure 'data/PDB/blastdb/GPCRDB' exists before running this
    os.system(f'blastp -query {fasta_file} -outfmt 5 -out {blast_out} -db data/PDB/blastdb/GPCRDB')

    # 5. Process BLAST Results (Map PDB to GPCRdb numbers)
    if not os.path.exists(blast_out) or os.path.getsize(blast_out) == 0:
        print(f"⚠️ BLAST failed for {pdbID} (No output file)")
        return fasta

    try:
        handle = open(blast_out, 'r')
        blast_records = NCBIXML.parse(handle)
        
        GPCRDB2SEQ = {}
        bestHIT = 'None'
        chainID = target_chain

        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    bestHIT = alignment.title.split(' ')[1]
                    num_s = hsp.sbjct_start
                    num_q = hsp.query_start
                    
                    # Map residue by residue
                    for q, s in zip(hsp.query, hsp.sbjct):
                        if q != '-' and s != '-':
                            GPCRDB2SEQ[num_s] = num_q
                            num_q += 1
                            num_s += 1
                        elif q != '-':
                            num_q += 1
                        else:
                            num_s += 1
                    break # Take only the first HSP
                break # Take only the best alignment

        # Write the final mapping file
        l = ''
        for s in GPCRDB2SEQ:
            q = GPCRDB2SEQ[s]
            if q in SEQ2PDB:
                pdbPosition = SEQ2PDB[q]
                l += f"{pdbID}\t{chainID}\t{pdbPosition}\t{s}\t{bestHIT}\n"
        
        with open(f'{gpcrdb_dir}/{pdbID}.txt', 'w') as f:
            f.write(l)
            
        print(f"✅ Success: {pdbID} mapped to {bestHIT}")

    except Exception as e:
        print(f"❌ Error parsing BLAST XML for {pdbID}: {e}")

    return fasta

## Save PDB ID and chain information as well as
## fetch FASTA files of PDB IDs
l = ''
pdblist = []; allPDB = ''
for pdbID in dic:
    #if pdbID == '7f58' or pdbID == '7F58':
    if True:
        if (dic[pdbID]['GPCR'] != '-' and dic[pdbID]['GPROT'] != '-') or (dic[pdbID]['GPCR'] != '-' and dic[pdbID]['BARR'] != '-'):
            l += pdbID + ' ' + dic[pdbID]['GPCR'] + ' ' + dic[pdbID]['GPROT'] + ' ' + dic[pdbID]['BARR'] + '\n'
            #os.system('wget https://www.rcsb.org/fasta/entry/'+pdbID.lower()+'/download' + ' -O ../data/fasta/'+ pdbID.lower()+'.fasta')
            allPDB += makeMapFASTA(pdbID, dic) + '\n'
            pdblist.append(pdbID)
            #break

#sys.exit()
#print (l)
open('data/PDB/pdblist.txt', 'w').write(l)

## Concatenate all FASTA files into one and make a blastdb in blastdb/
os.system("rm -rf data/PDB/allPDB.fasta")
open ('data/PDB/allPDB.fasta', 'w').write(allPDB)
if os.path.isfile('data/PDB/blastdb/') == False:
    os.system("mkdir data/PDB/blastdb/")
os.system("rm -rf data/PDB/blastdb/allPDB*")
os.system("makeblastdb -in data/PDB/allPDB.fasta -dbtype 'prot' -out data/PDB/blastdb/allPDB")
print ('complete')
