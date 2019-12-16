#!/bin/bash
#SBATCH --job-name=HISAT2_index_build                                  # Job name
#SBATCH --mail-type=END                               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email>                           # Where to send mail
#SBATCH --nodes=1                                          # Number of nodes -- almost always one
#SBATCH --ntasks=8                                         # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=06:00:00                                    # Time limit hrs:min:sec
#SBATCH --partition=short                                # Queue you want to submit to
#SBATCH --mem=32gb                                         # Memory limit
#SBATCH --output=<Path to OUTERR directory>/HISAT2_index_build_%j.out
#SBATCH --error=<Path to OUTERR directory>/HISAT2_index_build_%j.err

######################################################
# Include info in output file to help with troubleshooting

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors


##### SPECIFY JOB SPECIFIC VARIABLES ##############

GENOME='/scratch/Users/lysa8537/Genomes/mm10'

##### LOAD NECESSARY MODULES ######################

module load samtools/1.3.1
module load hisat2/2.1.0

###################################################

#Sync over genome fa to my scratch

rsync /scratch/Shares/public/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa $GENOME/

#Index reference genome

hisat2-build -p 8 -f $GENOME/genome.fa mm10

#Final time stamp
echo Time is `date`
