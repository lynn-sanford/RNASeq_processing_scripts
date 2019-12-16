library("Rsubread")

# featurecounts_file_list.txt is a file that contains all of the filenames of the mapped BAM files to count (see examples in scripts folder)
# This can be generated from featureCounts.sbatch script
inputData <- read.table("/Users/lysa8537/2016_07_RNASeq/2018_analysis/counting/featurecounts_file_list.txt")

# Transpose this file list and put it in a variable
fileList <- t(inputData$V1)

# featurecounts_sample_list.txt is a file that contains all of the headings to make it readable - must be in same order as mapped file list
# This can be generated from featureCounts.sbatch script
columnNames <- read.table("/Users/lysa8537/2016_07_RNASeq/2018_analysis/counting/featurecounts_sample_list.txt")

# Transpose this sample list and put it in a variable
sampleList <- t(columnNames$V1)

# Run featureCounts on these BAM files - look up featureCounts manual to understand parameters
coverage <- featureCounts(fileList,
# annotation file is the GTF file containing the genes over which to count
annot.ext="/Users/lysa8537/Genomes_and_annotations/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf",
isGTFAnnotationFile=TRUE,
GTF.featureType="exon",
GTF.attrType="gene_id",
# level of summarization
useMetaFeatures=TRUE,
# overlap between reads and features
allowMultiOverlap=FALSE,
largestOverlap=FALSE,
# strandness
strandSpecific=2,
# parameters specific to paired end reads
isPairedEnd=TRUE,
requireBothEndsMapped=TRUE,
checkFragLength=TRUE,
minFragLength=50,
maxFragLength=60000,
countChimericFragments=FALSE,
autosort=TRUE,
# number of CPU threads
nthreads=1)

# Write out annotation and count data as tab delimited txt files
# Count data is just counts
# Annotation data has the gene length information for calculating TPM

colnames(coverage$counts) <- sampleList

write.table(coverage$counts,file="/Users/lysa8537/2016_07_RNASeq/2018_analysis/counting/rat_gene_count_s2.txt",append=FALSE,sep='\t',quote=FALSE)
write.table(coverage$annotation,file="/Users/lysa8537/2016_07_RNASeq/2018_analysis/counting/rat_gene_count_annotation_s2.txt",append=FALSE,sep='\t',quote=FALSE)
