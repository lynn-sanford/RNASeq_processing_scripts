#!/bin/bash
#SBATCH --job-name=HISAT2_acc_cont1                                  # Job name
#SBATCH --mail-type=END                               # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email>                           # Where to send mail
#SBATCH --nodes=1                                          # Number of nodes -- almost always one
#SBATCH --ntasks=8                                         # Number of CPU (processer cores i.e. tasks)
#SBATCH --time=02:00:00                                    # Time limit hrs:min:sec
#SBATCH --partition=short                                # Queue you want to submit to
#SBATCH --mem=32gb                                         # Memory limit
#SBATCH --output=<Path to output directory>/HISAT2_%j.out
#SBATCH --error=<Path to error directory>/HISAT2_%j.err


####################################################
# Include information in output file that may help troubleshooting
# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors


##### SPECIFY JOB SPECIFIC VARIABLES ##############

SCRATCH='/scratch/Users/lysa8537/2018_RNASeq/accumulated_data'
MAPPING='/Users/lysa8537/2018_RNASeq/accumulated_data/04_mapping/map2'
GENOME='/scratch/Users/lysa8537/Genomes/mm10/HISAT2Index'
BASENAME='accum_cont1.map2'

##### LOAD NECESSARY MODULES ######################

module load samtools/1.3.1
module load hisat2/2.1.0

###################################################

# Map trimmed fastq files to reference genome --> sam files include all reads
# This assumes a HISAT2 index has already been created - see separate indexing script

hisat2 -p 8 -q --phred33 \
 --very-sensitive --max-intronlen 50000 --rna-strandness FR \
 -x $GENOME/mm10 \
 -1 $SCRATCH/trimmed_fastq/accum_Control_1_S20_R1_trimmed.fastq.gz \
 -2 $SCRATCH/trimmed_fastq/accum_Control_1_S20_R2_trimmed.fastq.gz \
 > $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.sam \
 2> $MAPPING/${BASENAME}.hisat2.stderr

# Convert SAM to BAM and sort BAM by position
samtools view -b -o $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.bam \
 $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.sam

samtools sort -T $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.sorted \
 -o $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.sorted.bam \
 $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.bam

# Index BAM

samtools index $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.sorted.bam

# Transfer BAM back to mapping directory

rsync $SCRATCH/HISAT_mapping/${BASENAME}.hisat2.sorted.* $MAPPING/

#Final time stamp
echo Time is `date`
