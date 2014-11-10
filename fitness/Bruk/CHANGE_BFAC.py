#! /usr/bin/python
import numpy as np
import sys
from Bio.PDB import PDBIO
from Bio.PDB.PDBParser import PDBParser
parser=PDBParser()
parser.__init__(QUIET=True)


## USAGE:
# ./CHANGE_BFAC.py [fitvals.txt] [pdb file] [ouput name]

fit_vals = open(sys.argv[1], 'r').readline().split()
print len(fit_vals)
fit_vals = np.array(fit_vals, dtype='float')/10.


# load PDB file
structure = parser.get_structure(sys.argv[2], sys.argv[2])
model = structure[0]
output_name = sys.argv[3]

# change beta
for chain in model:
  for residue in chain:
     for atom in residue:
        atom.bfactor = 0.0

for i, fitness in enumerate(fit_vals):
  for chain in model:
     for residue in chain:
        if residue.id[1]==i+1:
           for atom in residue:
              atom.bfactor=fitness

w = PDBIO()
w.set_structure(structure)
w.save("{0}".format(output_name))
