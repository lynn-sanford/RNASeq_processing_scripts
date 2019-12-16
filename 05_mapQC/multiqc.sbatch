#!/bin/bash
#SBATCH --job-name=multiQC 	# Job name
#SBATCH --mail-type=END 	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 		# Number of CPU (processer cores i.e. tasks)
#SBATCH --time=00:30:00 	# Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=1gb 		# Memory limit
#SBATCH --output=<Path to OUTERR directory>/multiqc_%j.out
#SBATCH --error=<Path to OUTERR directory>/multiqc_%j.err

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors


##### SPECIFY JOB SPECIFIC VARIABLES ##############

module load python/2.7.14

###################################################
# This is only useful if you've run the preseq and RSeQC scripts first - it accumulates all of that data

/Users/lysa8537/.local/bin/multiqc -o /Users/lysa8537/2018_RNASeq/accumulated_data/05_mapping_qc/ /Users/lysa8537/2018_RNASeq/accumulated_data/

echo date
