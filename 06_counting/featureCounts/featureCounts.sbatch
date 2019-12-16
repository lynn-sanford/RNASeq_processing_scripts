#!/bin/bash
#SBATCH --job-name=featureCounts      # Job name
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<InsertEmail> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=04:00:00         # Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=4gb               # Memory limit
#SBATCH --output=<PathToOutDirectory>/featureCounts_%j.out
#SBATCH --error=<PathToErrorDirectory>/featureCounts_%j.err


##### OUTPUT RUN PARAMETERS #######################
# These go into output file and may help in troubleshooting later

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors

##### LOAD MODULES AND RLIB LOCATION ##############
# See /Shares/palmer/Lynn/R_packages_README.txt for general explanation in getting R libraries into your home directory
# Change RLIB path to the location of those libraries in your home directory

module load R/3.3.0
export RLIB='/Users/lysa8537/R/x86_64-pc-linux-gnu-library/3.3/'

##### SPECIFY VARIABLES ###########################

COUNTDIR='/Users/lysa8537/2016_07_RNASeq/2018_analysis/counting'
BAMDIR='/scratch/Users/lysa8537/2018_RNASeq/accumulated_data/HISAT_mapping'

###################################################
# Run this to create file lists and sample lists - the numbers in the ${filename:73:5} correspond to the string locations for sample IDs in my BAMfile names - they will differ for yours

# Create file with all filenames to analyze - done previously
# for filename in $BAMDIR/*.map2.hisat2.bam; do
# echo $filename >> $COUNTDIR/map2_file_list.txt
# echo ${filename:73:5} >> $COUNTDIR/map2_sample_list.txt
# done

###################################################
# Between the defined COUNTDIR variable and this path below, must be the location of the corresponding R script to run featureCounts

### Run featureCounts by running rscript
Rscript $COUNTDIR/SCRIPTS/featureCounts.r

#Final time stamp
echo Time is `date`
