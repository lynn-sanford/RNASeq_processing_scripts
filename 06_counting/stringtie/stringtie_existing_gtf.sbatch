#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=stringtie_batch				# Job name
#SBATCH --mail-type=END						# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=<Insert email>			# Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=8	     					# Run on <number> CPUs
#SBATCH --mem=8gb						# Memory limit
#SBATCH --time=10:00:00						# Time limit hrs:min:sec
#SBATCH --output=<Path to OUTERR>/stringtie_exgtf_%j.out
#SBATCH --error=<Path to OUTERR>/stringtie_exgtf_%j.err

### This script runs stringtie and tries a bunch of outputs...

##### OUTPUT RUN PARAMETERS #######################

# Calculate the number of processors allocated to this run.
NPROCS=$SLURM_CPUS_ON_NODE

### Display the job context
echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`
echo Using ${NPROCS} processors

##### LOAD MODULES AND RLIB LOCATION ##############

module load stringtie/1.3.3b 

##### SPECIFY VARIABLES ###########################

COUNTDIR='/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting'
BAMDIR='/scratch/Users/lysa8537/2018_RNASeq/accumulated_data/HISAT_mapping'
GTF='/Users/lysa8537/Genomes_and_annotations/Mus_musculus/genes.gtf'

###################################################
### STRINGTIE OPTIONS ############################################
# -v 			Verbose mode
# -o <fullpath.gtf>	Sets name of output gtF file in which StringTie will write the assembled transcripts
# -A <gene_abund.tab>	Report gene abundances in output file with given name (tab delimited)
# -G <ref_ann.gtf>	Use reference annotation file (in GTF or GFF3 format) to guide assembly. Output includes expressed reference transcripts as well as novel assembled transcripts. Requires -B, -b, -e, -C
# -C			Output a file with all transcripts in the provided gtf reference that are covered by reads
# -B			Enables output of Ballgown input table containing coverage data (DIRECT REPLACEMENT FOR HAVING TO CALL TABLEMAKER)
# -e			limits processing of read alignments to only estimate and output assembled transcripts matching reference transctipts given by -G option. This will only consider read bundles that map to annotated transcript (MIGHT WANT TO SKIP THIS)
# -b			Specifies the output of the Ballgown file (-B flag), if different from -o location

#stringtie ${SORTEDBAM} -vB -G ${HG38GTF} -o ${OUTGTF} -C ${COVGTF} -A ${GENEABUNDANCE}

#################################################
### CALCULTE COVERAGE WITH RESPECT TO GTF

# This loops based on a file list given to it after the "done" command
while read -r line; do

  # Define a BAM filename
  BAM=${line:0:-3}'sorted.bam'

  ### Define output files based on sample name in the filename string
  NEWOUTPATH=${COUNTDIR}/stringtie_exgtf/${line:73:10}'/'
  mkdir ${NEWOUTPATH}

  OUTGTF=${NEWOUTPATH}${line:73:10}'_stringtie_FINAL.gtf'
		
  # Recalculate expression levels using new merged GTF
  stringtie -p 8 -e -B -G ${GTF} -o ${OUTGTF} ${BAM}

done < $COUNTDIR/map2_file_list.txt

echo Time is `date`
