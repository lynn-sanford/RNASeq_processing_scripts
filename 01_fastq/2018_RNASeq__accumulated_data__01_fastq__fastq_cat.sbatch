#!/bin/bash
#SBATCH --partition=short              # Partition or queue. In this case, short!
#SBATCH --job-name=fastq_cat    # Job name
#SBATCH --mail-type=END               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email>
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=1gb                   # Memory limit
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=<Path to OUTERR>/fastq_cat_%j.out   # Standard output and error log
#SBATCH --error=<Path to OUTERR>/fastq_cat_%j.err   # %j inserts job number

######################################

# This script was to concatenate fastq files from two separate rounds of sequencing into "accumulated files" - one for each sample that is basically one fastq file right after the other
# If you don't do any other data manipulation first, paired read order should still be fine upon concatenating


# For this I had to sync one set of my two fastq files to scratch - others were already there
rsync /Users/lysa8537/2018_RNASeq/2018_07_29/01_fastq/*.fastq.gz /scratch/Users/lysa8537/2018_RNASeq/2018_07_29/fastq/

# Takes each file in one dataset and concatenates with the filename that matches in the other dataset
for FASTQ in /scratch/Users/lysa8537/2018_RNASeq/2018_07_29/fastq/*.fastq.gz; do
cat ${FASTQ} /scratch/Users/lysa8537/2018_RNASeq/2018_08_13/fastq/${FASTQ##*/} > /scratch/Users/lysa8537/2018_RNASeq/accumulated_data/fastq/accum_${FASTQ##*/}
done

rsync /scratch/Users/lysa8537/2018_RNASeq/accumulated_data/fastq/accum_* /Users/lysa8537/2018_RNASeq/accumulated_data/01_fastq/
