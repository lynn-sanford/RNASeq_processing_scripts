#!/bin/bash
#SBATCH --partition=short              # Partition or queue. In this case, short!
#SBATCH --job-name=picard    # Job name
#SBATCH --mail-type=END               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email>
#SBATCH --nodes=1                    # Only use a single node
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=16gb                   # Memory limit
#SBATCH --time=00:10:00               # Time limit hrs:min:sec
#SBATCH --output=<Path to output directory>/insertsize_%j.out   # Standard output and error log
#SBATCH --error=<Path to error directory>/insertsize_%j.err   # %j inserts job number

### Run Picard tools CollectInsertSizeMetrics
# This generally describes the distribution of reads across transcripts
# See Picardtools documentation for more information

java -jar /opt/picard/2.6.0/picard-2.6.0.jar CollectInsertSizeMetrics I=/Users/lysa8537/2016_07_RNASeq/Samples/KCl_only/Rep1/04_mapping/KCl_Rep1_trim_beg10.tophat.igenome.bam O=/Users/lysa8537/2016_07_RNASeq/Samples/KCl_only/Rep1/04_mapping/KCl_Rep1_trim_beg10.igenome.insert_size_metrics.txt H=/Users/lysa8537/2016_07_RNASeq/Samples/KCl_only/Rep1/04_mapping/KCl_Rep1_trim_beg10.igenome.insert_size_hist.pdf

#Final time stamp
echo Time is `date`
