library("DESeq2")

# Generate data frames with treatment and replicate information, and normalization factor info
# Replace Rep data with your replicates and Treatment with your sample names

colData_all <- data.frame(Rep=c(1,2,3,1,2,3,1,2,3),Treatment=c("kclo","kclo","kclo","kclt","kclt","kclt","kclz","kclz","kclz"))
rownames(colData_all) <- c("kclo1","kclo2","kclo3","kclt1","kclt2","kclt3","kclz1","kclz2","kclz3")

# Load in count and annotation data from files
# These are the files output from featureCounts

countTable <- read.table("/Users/lysa8537/2016_07_RNASeq/2018_analysis/counting/rat_gene_count_s2.txt",header=TRUE,sep='\t')
annotData <- read.table("/Users/lysa8537/2016_07_RNASeq/2018_analysis/counting/rat_gene_count_annotation_s2.txt",header=TRUE,sep='\t')

# Generate DESeq count matrix

dds <- DESeqDataSetFromMatrix(countData = countTable, colData = colData_all, design = ~ Treatment)

# Calculate TPM values
# This can be done independently of DESeq2 with the count file, as well
gene_counts <- counts(dds)
RPK <- (gene_counts)/(annotData$Length/1000)
per_mil_scaling_factor <- (colSums(RPK)/1000000)
TPM <- t(t(RPK)/per_mil_scaling_factor) 
rownames(TPM) <- annotData$GeneID

# Write out TPM file
write.table(as.data.frame(TPM),file="/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/all_samples_TPM.txt",sep="\t",quote=FALSE)

# Calculate version of the dds table for cluster analysis between samples
# This is for a rough analysis to visualize a dendrogram and PCA of samples - more robust methods are recommended for further analysis

rld <- rlog(dds, blind = FALSE)

# Calculate Euclidean distances

sampleDists <- dist(t(assay(rld)))

# Print out dendrogram

png("/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_all_samples_dendrogram.png",type='cairo')
plot(hclust(sampleDists))
dev.off()
graphics.off()

# Print out PCA

png("/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_all_samples_PCA.png",type='cairo')
plotPCA(rld,intgroup='Treatment')
dev.off()
graphics.off()

# Do differential expression analysis and split up pairwise comparisons
# You may have more rows depending on how many pairwise comparisons you have
# Ensures that regardless of the annotation labeling, you'll get the gene identifier to label rows
# Can compare replicates instead of treatment conditions to see if any replicates are way off

dds <- DESeq(dds)

res_kclo_kclt <- results(dds, contrast=c("Treatment","kclt","kclo"))
res_kclo_kclt_labeled <- data.frame(res_kclo_kclt,gene_symbol = annotData$GeneID)
res_kclo_kclz <- results(dds, contrast=c("Treatment","kclz","kclo"))
res_kclo_kclz_labeled <- data.frame(res_kclo_kclz,gene_symbol = annotData$GeneID)
res_kclt_kclz <- results(dds, contrast=c("Treatment","kclz","kclt"))
res_kclt_kclz_labeled <- data.frame(res_kclt_kclz,gene_symbol = annotData$GeneID)

# Write results to files
# You may have more files depending on comparisons

write.table(as.data.frame(res_kclo_kclt_labeled),file="/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_gene_kclo_kclt.txt",sep="\t",quote=FALSE)
write.table(as.data.frame(res_kclo_kclz_labeled),file="/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_gene_kclo_kclz.txt",sep="\t",quote=FALSE)
write.table(as.data.frame(res_kclt_kclz_labeled),file="/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_gene_kclt_kclz.txt",sep="\t",quote=FALSE)

# Make and write out several MA plots with different p-values
# These are a little rough and not very customizable - can use ggplot2 to make better ones

png("/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_gene_kclo_kclt.p05.png",type='cairo')
plotMA(res_kclo_kclt, alpha=0.05, main="DESeq2", ylim=c(-1,1))
dev.off()
graphics.off()

png("/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_gene_kclo_kclz.p05.png",type='cairo')
plotMA(res_kclo_kclz, alpha=0.05, main="DESeq2", ylim=c(-1,1))
dev.off()
graphics.off()

png("/Users/lysa8537/2016_07_RNASeq/2018_analysis/diffex/DESeq2_gene_kclt_kclz.p05.png",type='cairo')
plotMA(res_kclt_kclz, alpha=0.05, main="DESeq2", ylim=c(-1,1))
dev.off()
graphics.off()
