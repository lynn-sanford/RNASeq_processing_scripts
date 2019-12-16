#!/bin/bash
#SBATCH --job-name=Preseq_cc 	# Job name
#SBATCH --mail-type=END 	# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1 		# Number of CPU (processer cores i.e. tasks)
#SBATCH --time=12:00:00 	# Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=8gb 		# Memory limit
#SBATCH --output=<Path to OUTERR directory>/preseq_ccurve_%j.out
#SBATCH --error=<Path to OUTERR directory>/preseq_ccurve_%j.err

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

#c_curve : can use inputs from file types [bed/bam] and will plot the estimated complexity of a sample. Specify -B for sorted bam files, -P for paired end read files

for sortedbam in $SCRATCH/*.sorted.bam; do
preseq c_curve -B \
-o $OUTDIR/${sortedbam:73:5}_map1_preseq_ccurve.txt $sortedbam
done

echo `date`
