#! /usr/local/bin/python



import __main__
__main__.pymol_argv = [ 'pymol', '-qc'] # Quiet and no GUI

import sys, time, os
import pymol
from pymol import cmd

pymol.finish_launching()

##python [structure.py] [your input structure (in our case,fintness_scoreB_factor pdb file)]
# Read User Input
spath = os.path.abspath(sys.argv[1])
sname = spath.split('/')[-1].split('.')[0]

# Load Structures in pymol

pymol.cmd.load(spath, sname)
pymol.cmd.disable("all")
pymol.cmd.enable(sname)

pymol.cmd.bg_color("white") #change color to white
pymol.cmd.remove("solvent") #remove solvent
pymol.cmd.show_as("cartoon")
pymol.cmd.hide("sticks")
pymol.cmd.color("gray50")

#change color by b factor
pymol.cmd.spectrum("b", blue_red, minimum=-0.5, maximum=0)

#Change residues with most significant b-factors to a colored sphere
cutoff = ("residue 1-10")
pymol.cmd.show("spheres",cutoff)
pymol.cmd.color("marine",cutoff)
pymol.cmd.set_view((\
	"0.215765774,    0.370535314,    0.903409004,\
	-0.027041638,    0.927114367,   -0.373799652,\
    -0.976070225,    0.056223251,    0.210059240,\
    0.000000000,    0.000000000, -119.776687622,\
    30.295255661,   28.771011353,   15.131137848,\
    94.432853699,  145.120529175,  -20.000000000"))


pymol.cmd.png("test_image.png") #save image to a new file

# Close pymol 
pymol.cmd.quit()