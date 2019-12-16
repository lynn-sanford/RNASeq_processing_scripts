#!/bin/bash
#SBATCH --job-name=preseq_lcextrap 	# Job name
#SBATCH --mail-type=END 	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 		# Number of CPU (processer cores i.e. tasks)
#SBATCH --time=12:00:00 	# Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=8gb 		# Memory limit
#SBATCH --output=<Path to OUTERR directory>/preseq_lcextrap_%j.out
#SBATCH --error=<Path to OUTERR directory>/preseq_lcextrap_%j.err

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

module load preseq/2.0.3

###################################################

#lc_extrap : can use inputs from from tyles [bed/bam] and will estimate future yields for complexity if the sample is sequenced at higher read depths. Specify -B for sorted bam files, -P for paired end reads

for sortedbam in $SCRATCH/*.sorted.bam; do
preseq lc_extrap -B \
-o $OUTDIR/${sortedbam:73:5}_map1_preseq_lc_extrap.txt \
$sortedbam
done

echo `date`
