#!/usr/bin/env python3

import sys
from Bio.PDB import PDBParser, MMCIFIO

name=sys.argv[1]
p = PDBParser()
struc = p.get_structure("", name+".pdb")
io = MMCIFIO()
io.set_structure(struc)
io.save(name+".cif")

