#!/bin/bash
#SBATCH --partition=short              # Partition or queue. In this case, short!
#SBATCH --job-name=rsync_md5sum    # Job name
#SBATCH --mail-type=END               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email>
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=8gb                   # Memory limit
#SBATCH --time=01:00:00               # Time limit hrs:min:sec
#SBATCH --output=<Path to OUTERR directory>/rsync_md5sum_%j.out   # Standard output and error log
#SBATCH --error=<Path to OUTERR directory>/rsync_md5sum_%j.err   # %j inserts job number

rsync /Users/lysa8537/2018_RNASeq/accumulated_data/01_fastq/*.fastq.gz /Users/lysa8537/GEO_mouse_Sanford_upload/

for FASTQ in /Users/lysa8537/GEO_mouse_Sanford_upload/*.fastq.gz; do
md5sum ${FASTQ} >> /Users/lysa8537/GEO_mouse_Sanford_upload/md5sum_list.txt
done
