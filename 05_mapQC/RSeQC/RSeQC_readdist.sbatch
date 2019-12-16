#!/bin/bash
#SBATCH --job-name=RSeQC_readdist 	# Job name
#SBATCH --mail-type=END 	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 		# Number of CPU (processer cores i.e. tasks)
#SBATCH --time=08:00:00 	# Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=8gb 		# Memory limit
#SBATCH --output=<Path to OUTERR directory>/RSeQC_readdist_%j.out
#SBATCH --error=<Path to OUTERR directory>/RSeQC_readdist_%j.err

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors


##### SPECIFY JOB SPECIFIC VARIABLES ##############

SCRATCH='/scratch/Users/lysa8537/2018_RNASeq/accumulated_data/HISAT_mapping'
OUTDIR='/scratch/Users/lysa8537/2018_RNASeq/accumulated_data/map_qc'

##### LOAD NECESSARY MODULES ######################

module load python/2.7.14/rseqc

###################################################

#analyze read distributions using RSeQC
# Makes a separate output file for each BAM based on part of the string of the input filename that designates sample

for sortedbam in $SCRATCH/*.sorted.bam; do
read_distribution.py \
-i $sortedbam \
-r /scratch/Users/lysa8537/Genomes/mm10/Annotations/Genes/mm10.bed \
> $OUTDIR/${sortedbam:73:5}_map1_RSeQC_readdist.txt
done

echo `date`
