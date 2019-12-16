#!/bin/bash
#SBATCH --job-name=DESeq2      # Job name
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email> # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=4gb               # Memory limit
#SBATCH --output=<Path to out directory>/DESeq2_%j.out
#SBATCH --error=<Path to error directory>/DESeq2_%j.err


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

SCRIPTDIR='/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex'

###################################################
# Between the defined SCRIPTDIR variable and this path below, must be the location of the corresponding R script to run DESeq2

### Run featureCounts by running rscript
Rscript $SCRIPTDIR/SCRIPTS/DESeq2.r

#Final time stamp
echo Time is `date`
