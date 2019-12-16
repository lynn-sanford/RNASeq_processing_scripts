#!/bin/bash
#SBATCH --partition=short                       # Partition or queue. short, long, or highmem
#SBATCH --job-name=fastQC                       # Job name
#SBATCH --mail-type=END                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert Email>
#SBATCH --nodes=1                               # Only use a single node (should always be 1)
#SBATCH --ntasks=1                              # Run on a single CPU (only more than 1 for parallel computing)
#SBATCH --mem=4gb                             # Memory limit
#SBATCH --time=04:00:00                         # Time limit hrs:min:sec
#SBATCH --output=<Path to out directory>/fastQC_%j.out   # Standard output and error log
#SBATCH --error=<Path to error directory>/fastQC_%j.err    # %j inserts job number

##################################################
# Put information into output that may be useful for troubleshooting

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors

##################################################

### Make shortcuts
SCRATCH='/scratch/Users/lysa8537/2018_RNASeq/2018_08_13'

##################################################
### Load fastqc module
module load fastqc/0.11.5

##################################################
### Move all fastq inputs to scratch
rsync /Users/lysa8537/2018_RNASeq/2018_08_13/01_fastq/*.fastq.gz $SCRATCH/fastq/

### Run fastqc on fastq files already in scratch, output to QC folder
for FASTQ in $SCRATCH/fastq/*.fastq.gz; do
 fastqc -o $SCRATCH/fastQC/ $FASTQ
done

### Move outputs back to QC repository
rsync -a $SCRATCH/fastQC/* /Users/lysa8537/2018_RNASeq/2018_08_13/02_fastQC/

#Final time stamp
echo Time is `date`
