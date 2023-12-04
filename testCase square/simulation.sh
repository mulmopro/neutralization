#!/bin/bash
#SBATCH --mail-user=antonello.raponi@polito.it
#SBATCH --mail-type=FAIL    	                # NONE, BEGIN, END, FAIL, REQUEUE, ALL, TIME_LIMIT, TIME_LIMIT_90
#SBATCH --error=simulation.err                 # standard error file
#SBATCH --output=simulation.out                # standard output file
#SBATCH --account=araponi                      # account name
#SBATCH --mem=5                   	        # Memory Size (maximum of 128G per node)
#SBATCH --time=24:00:00   	                # Time-limit
######################################

#SBATCH --job-name=mesh4TMixerDH3mm       # Name of job
#SBATCH -p global                         # Queue
#SBATCH -N1 -n16                           # <n> cores on <N> nodes

## Please note that Casper has 17 nodes each with 36 cores (AMD Bulldozer)
##              and Hactar has 24 nodes each with 24 cores (Intel XEON v3)

######################################

module load openmpi3
module load openfoam/8.0
source /share/apps/openfoam/8/OpenFOAM-8/etc/bashrc
module load intel/python/3

/bin/echo " "
/bin/echo Running on host: $SLURM_SUBMIT_HOST
/bin/echo In directory: $SLURM_SUBMIT_DIR 
/bin/echo "Got $SLURM_CPUS_ON_NODE processors."
/bin/echo " "

echo "$SLURM_CPUS_ON_NODE"

## myJob script
# ./Allclean
./Allrun $SLURM_CPUS_ON_NODE
