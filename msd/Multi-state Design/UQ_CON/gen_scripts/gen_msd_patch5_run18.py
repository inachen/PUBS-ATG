from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)
from Bio.PDB import NeighborSearch
from Bio.PDB import Selection
#You can only run this script while ssh'ed into iqint because biopython is not installed on the normal submission login.


radius = 10 # size of repack shell

patches = {}
# patchname = "patch1"
# residue_list = [1,19,20,57,60,62,63] #list of integers, assumes chain A
# patches[patchname] = residue_list

# patchname = "patch2"
# residue_list = [42,44,45,46,48,49] #list of integers, assumes chain A
# patches[patchname] = residue_list

# patchname = "patch3"
# residue_list = [6,8,42,44,66,68,70,73] #list of integers, assumes chain A
# patches[patchname] = residue_list

# patchname = "patch4"
# residue_list = [44, 45, 57, 58, 59, 60, 62, 65] #list of integers, assumes chain A
# patches[patchname] = residue_list

patchname = "patch5"
residue_list = [8, 42, 44, 68, 70, 72] #list of integers, assumes chain A
patches[patchname] = residue_list

structures = {}
#The keys for this dictionary are strings and the string cannot begin with a number (due to parsing issues in Rosetta).
structures["UBQ"] = "/netapp/home/thompson/UBQ/PDBs/1UBQ_yeast_relax.pdb"
structures["CUE"] = "/netapp/home/thompson/UBQ/PDBs/CUE_yeast_relax.pdb"
structures["OTU"] = "/netapp/home/thompson/UBQ/PDBs/OTU_yeast_relax.pdb"
structures["RPN13"] = "/netapp/home/thompson/UBQ/PDBs/RPN13_yeast_relax.pdb"
structures["SH3"] = "/netapp/home/thompson/UBQ/PDBs/SH3_yeast_relax.pdb"
structures["UQ_CON"] = "/netapp/home/thompson/UBQ/PDBs/UQ_con_yeast_relax.pdb"

startenergies = {}
#Note that if you use your own patch, these should to be optimized.
#The current values are the energies from the minimized structures, which are too low for these simulations.
#Talk to Samuel for help with the optimization.
#Make sure that the keys for the startenergies dictionary match those for the structures dictionary 
startenergies["patch1"] = {}
startenergies["patch1"]["UBQ"] = "-127.264"
startenergies["patch1"]["CUE"] = "-193.851"
startenergies["patch1"]["OTU"] = "-436.815"
startenergies["patch1"]["RPN13"] = "-293.743"
startenergies["patch1"]["SH3"] = "-205.083"
startenergies["patch1"]["UQ_CON"] = "-580.214"

startenergies["patch2"] = {}
startenergies["patch2"]["UBQ"] = "-128.23"
startenergies["patch2"]["CUE"] = "-191.393"
startenergies["patch2"]["OTU"] = "-427.361"
startenergies["patch2"]["RPN13"] = "-289.546"
startenergies["patch2"]["SH3"] = "-196.747"
startenergies["patch2"]["UQ_CON"] = "-589.052"

startenergies["patch3"] = {}
startenergies["patch3"]["UBQ"] = "-126.755"
startenergies["patch3"]["CUE"] = "-190.438"
startenergies["patch3"]["OTU"] = "-413.491"
startenergies["patch3"]["RPN13"] = "-289.394"
startenergies["patch3"]["SH3"] = "-216.911"
startenergies["patch3"]["UQ_CON"] = "-588.108"

startenergies["patch4"] = {}
startenergies["patch4"]["UBQ"] = "-125.813"
startenergies["patch4"]["CUE"] = "-187.539"
startenergies["patch4"]["OTU"] = "-429.016"
startenergies["patch4"]["RPN13"] = "-286.033"
startenergies["patch4"]["SH3"] = "-196.728"
startenergies["patch4"]["UQ_CON"] = "-577.719"

startenergies["patch5"] = {}
startenergies["patch5"]["UBQ"] = "-127.03"
startenergies["patch5"]["CUE"] = "-190.665"
startenergies["patch5"]["OTU"] = "-418.467"
startenergies["patch5"]["RPN13"] = "-290.882"
startenergies["patch5"]["SH3"] = "-217.574"
startenergies["patch5"]["UQ_CON"] = "-591.748"

k_values = {}
#These are the primary values that you should adjust in your experiments.
k_values["UBQ"] = "0"
k_values["CUE"] = "0"
k_values["OTU"] = "0"
k_values["RPN13"] = "0"
k_values["SH3"] = "1"
# OUR COMPOUND
k_values["UQ_CON"] = "1.0"

s_values = {}
#You should experiment with these values only if you have more time after screening k-values.
s_values["UBQ"] = "0.6"
s_values["CUE"] = "0.6"
s_values["OTU"] = "0.6"
s_values["RPN13"] = "0.6"
s_values["SH3"] = "0.6"
s_values["UQ_CON"] = "0.6"

#========== You should not need to edit below here ==========

for patch in patches.keys():
	residue_list = patches[patch]
	
	entityfilename = str("entity_file_%s.txt") % patch
	entityfile = open(entityfilename, 'w')
	patch_size = len(residue_list)
	outstring = str("%d\nALLAA EX 1 EX 2 EX_CUTOFF 0\nstart") % patch_size
	entityfile.write(outstring)
	entityfile.close()
	
	correspondencefilename = str("correspondence_file_%s.txt") % patch
	correspondencefile = open(correspondencefilename, 'w')
	for num in range(patch_size):
		entity_pos = num + 1
		PDB_pos = residue_list[num]
		outstring = str("%d %d A\n") % (entity_pos, PDB_pos)
		correspondencefile.write(outstring)
	correspondencefile.close()
	
	fitnessfilename = str("fitness_file_%s.txt") % patch
	fitnessfile = open(fitnessfilename, 'w')
	fitnessstring = "FITNESS -1"
	
	submissionfilename = str("MSD_MPI_%s.sh") % patch
	submissionfile = open(submissionfilename, 'w')
	outstring = str("""#
/bin/bash
#
#$ -S /bin/bash
#$ -l arch=linux-x64    # Specify architecture, required
#$ -l mem_free=4G       # Memory usage, required.  Note that this is per slot
#$ -pe ompi 4           # Specify parallel environment and number of slots, required
#$ -R yes               # SGE host reservation, highly recommended
#$ -V                   # Pass current environment to exec node, required
#$ -cwd                 # Current working directory
#$ -t 1					# Number of repeats
#$ -l h_rt=168:00:00    # Time limit

# Load OpenMPI-1.5 environment
module load openmpi-x86_64

# Run application
mpirun -np $NSLOTS /netapp/home/thompson/Rosetta/rosetta_2014.22/source/bin/mpi_msd.mpi.linuxgccrelease -database /netapp/home/thompson/Rosetta/rosetta_2014.22/database -msd:fitness_file %s -entity_resfile %s -ms:generations 200 -ms:pop_size 250 -ms:fraction_by_recombination 0.1 -msd::double_lazy_ig_mem_limit 1024 -ignore_unrecognized_res -ignore_zero_occupancy false -mute core.pack.annealer.MultiCoolAnnealer""") % (fitnessfilename, entityfilename)
	submissionfile.write(outstring)
	submissionfile.close()

	
	for structure_id in structures.keys():
	
		structurefilename = structures[structure_id]
		structure = parser.get_structure(structure_id, structurefilename)
		
		resfilename = str("sec_resfile_%s_%s.txt") % (structure_id, patch)
		resfile = open(resfilename, 'w')
		resfile.write("NATRO\nstart\n")
		
		startenergy = startenergies[patch][structure_id]
		s_value = s_values[structure_id]
		k_value = k_values[structure_id]
		
		outstring = str("""STATE %s_currentenergy %s %s %s 
SCALAR_EXPRESSION %s_startenergy = %s
SCALAR_EXPRESSION %s_o = %s_startenergy + 4.91
SCALAR_EXPRESSION %s_s = %s
SCALAR_EXPRESSION %s_exp = exp(%s_s*(%s_currentenergy-%s_o))
SCALAR_EXPRESSION %s_k = %s
SCALAR_EXPRESSION %s_sig = (1-%s_k)+(%s_k/(1+%s_exp))

""") %(structure_id, structurefilename, correspondencefilename, resfilename, structure_id, startenergy, structure_id, structure_id, structure_id, s_value, structure_id, structure_id, structure_id, structure_id, structure_id, k_value, structure_id, structure_id, structure_id, structure_id)
		fitnessfile.write(outstring)

		fitnessstring = fitnessstring + str("*%s_sig") % (structure_id)
		

		nucleosome = structure[0]

		atom_list = Selection.unfold_entities(nucleosome, 'A') # A for atoms
		neighbor_search = NeighborSearch(atom_list)

		contacts_list = neighbor_search.search_all(radius, level = 'R')
		
		repack_residues = []
		for contact in contacts_list:
			res1 = contact[0]
			res2 = contact[1]
			res1id = int(res1.get_id()[1])
			chain1 = res1.get_parent()
			chain1id = chain1.get_id()
			res2id = int(res2.get_id()[1])
			chain2 = res2.get_parent()
			chain2id = chain2.get_id()
			
			res1_in_patch = False
			res2_in_patch = False
			
			if (res1id in residue_list) and chain1id in "aA":
				res1_in_patch = True
			
			if (res2id in residue_list) and chain2id in "aA":
				res2_in_patch = True

			if res1_in_patch and not res2_in_patch:
				resstring = str("%d_%s") %(res2id, chain2id)
				repack_residues.append(resstring)
				
			if res2_in_patch and not res1_in_patch:	
				resstring = str("%d_%s") %(res1id, chain1id)
				repack_residues.append(resstring)
			
			repack_set = set(repack_residues)
		repack_set = sorted(repack_set)
		
		for residue in repack_set:
			pos, chain = residue.split("_")
			outstring = str("%s %s NATAA EX 1 EX 2 EX_CUTOFF 0\n") % (pos, chain)
			resfile.write(outstring)
		resfile.close()
	
	fitnessfile.write(fitnessstring)	
	fitnessfile.close()
	
	message = str("Now, run '$qsub %s'") % submissionfilename
	print message
				