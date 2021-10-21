These scripts were used in RNA-seq processing in Amy Palmer's lab prior to 2020



In these scripts, anything in <brackets> should be replaced without brackets
-Any paths should almost certainly be replaced - these were all for mine and will almost certainly not work with yours

Some scripts are PBS scripts from a previous cluster scheduler - would need to be converted to sbatch scripts


I generally did RNA-Seq analysis in these steps:

01_fastq - processing fastq files
02_fastQC - doing QC on fastq files
03_trimming - trimming reads on files to allow better mapping
04_mapping - mapping reads to genome
05_mapQC - doing QC on mapping to see performance
06_counting - counting reads over annotated areas of the genomes (genes)
07_diffex - doing differential expression analysis between conditions
08_downstream - prepping files for GSEA and DAVIDtools

00_general - any other types of scripts for syncing or processing files
