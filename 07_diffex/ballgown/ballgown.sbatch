#!/bin/bash
#SBATCH --job-name=Ballgown      # Job name
#SBATCH --mail-type=END         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=lynn.sanford@colorado.edu # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1              # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=01:00:00         # Time limit hrs:min:sec
#SBATCH --partition=short
#SBATCH --mem=4gb               # Memory limit
#SBATCH --output=/Users/lysa8537/2018_RNASeq/accumulated_data/07_diff_ex/OUTERR/ballgown_exgtf_accum_map2_%j.out
#SBATCH --error=/Users/lysa8537/2018_RNASeq/accumulated_data/07_diff_ex/OUTERR/ballgown_exgtf_accum_map2_%j.err


##### OUTPUT RUN PARAMETERS #######################

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors

##### LOAD MODULES AND RLIB LOCATION ##############

module load R/3.3.0
export RLIB='/Users/lysa8537/R/x86_64-pc-linux-gnu-library/3.3/'

##### SPECIFY VARIABLES ###########################

SCRIPTDIR='/Users/lysa8537/2018_RNASeq/accumulated_data/07_diff_ex'

###################################################


### Run ballgown by running rscript
Rscript $SCRIPTDIR/SCRIPTS/ballgown_exgtf_accum_map2.r

#Final time stamp
echo Time is `date`
