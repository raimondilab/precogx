import os, sys
from Bio.PDB import *

#pdbl = PDBList()
#   pdbl.retrieve_pdb_file("3SN6")

parser = MMCIFParser()
x = ''
structure = parser.get_structure(x, "data/PDB/pdir/6cmo.cif")
for model in structure:
    for chain in model:
        for residue in chain:
            for atom in residue:
                if atom.id == 'CA':
                    print(chain.id, residue.resname, residue.id(2), atom.id)
