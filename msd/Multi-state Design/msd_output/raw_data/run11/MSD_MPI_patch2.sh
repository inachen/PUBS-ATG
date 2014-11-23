#
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
mpirun -np $NSLOTS /netapp/home/thompson/Rosetta/rosetta_2014.22/source/bin/mpi_msd.mpi.linuxgccrelease -database /netapp/home/thompson/Rosetta/rosetta_2014.22/database -msd:fitness_file fitness_file_patch2.txt -entity_resfile entity_file_patch2.txt -ms:generations 200 -ms:pop_size 250 -ms:fraction_by_recombination 0.1 -msd::double_lazy_ig_mem_limit 1024 -ignore_unrecognized_res -ignore_zero_occupancy false -mute core.pack.annealer.MultiCoolAnnealer