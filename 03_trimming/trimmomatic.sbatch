#!/bin/bash
#SBATCH --partition=short                       # Partition or queue. short, long, or highmem
#SBATCH --job-name=trimmomatic                  # Job name
#SBATCH --mail-type=END                         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email>
#SBATCH --nodes=1                               # Only use a single node (should always be 1)
#SBATCH --ntasks=8                              # Run on a single CPU (only more than 1 for parallel computing)
#SBATCH --mem=16gb                               # Memory limit
#SBATCH --time=08:00:00                         # Time limit hrs:min:sec
#SBATCH --output=<Path to OUTERR directory>/trimmomatic_%j.out   # Standard output and error log
#SBATCH --error=<Path to OUTERR directory>/trimmomatic_%j.err    # %j inserts job number

### Load modules
module load trimmomatic/0.36 
module load java/1.8.0_101

### Define locations
SCRATCH='/scratch/Users/lysa8537/2018_RNASeq/accumulated_data'

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors

### Run trimmomatic on fastq files already in scratch
# I had to run a separate script for all of the fastq files in each condition because of the naming scheme from the core facility
# If you're trimming one file at a time, remove the loop and all of the string manipulations and just use the filename of that file (or files if paired end)

# This loop looks at all of the fastqs that I had and filters only the "Control" condition
for FASTQ in $SCRATCH/fastq/*R1_001.fastq.gz; do

# The numbers in the brackets refer to the portion of the filename string that contains the condition information
if [ ${FASTQ:65:7} = "Control" ]
then
java -jar /opt/trimmomatic/0.36/trimmomatic-0.36.jar PE -phred33 -threads 8 \
# This designates the R1 and R2 files for the same sample
$FASTQ ${FASTQ:0:79}R2_001.fastq.gz \
# These bracketed numbers refer to the portion of the filename string that I wanted to use to name my output files
$SCRATCH/trimmed_fastq/${FASTQ:59:20}R1_trimmed.fastq.gz \
$SCRATCH/trimmed_fastq/${FASTQ:59:20}R1_trimmed_unpaired.fastq.gz \
$SCRATCH/trimmed_fastq/${FASTQ:59:20}R2_trimmed.fastq.gz \
$SCRATCH/trimmed_fastq/${FASTQ:59:20}R2_trimmed_unpaired.fastq.gz \
# This command trims the first 10 bases off of every read
HEADCROP:10 \
# This command trimmed adapters off the end of reads - this file contains the adapters used for my library
ILLUMINACLIP:/opt/trimmomatic/0.36/adapters/TruSeq2-PE.fa:2:30:10 MINLEN:30
fi

done

#Final time stamp
echo Time is `date`
