# Libraries
library("ballgown")

# This is path where stringtie output is located
inpath <- "/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/"

# These sample ids must match stringtie outputs
filetable <- data.frame(ids=c('cont1.map2','cont2.map2','cont3.map2','cont4.map2','cont5.map2','glut1.map2','glut2.map2','glut3.map2','glut4.map2','glut5.map2','glzn1.map2','glzn2.map2','glzn3.map2','glzn4.map2','glzn5.map2','gtpa1.map2','gtpa2.map2','gtpa3.map2','gtpa4.map2','gtpa5.map2'),treatment=c('cont','cont','cont','cont','cont','glut','glut','glut','glut','glut','glzn','glzn','glzn','glzn','glzn','gtpa','gtpa','gtpa','gtpa','gtpa'),replicate=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5))

samples <- c('/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/cont1.map2/', 
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/cont2.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/cont3.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/cont4.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/cont5.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glut1.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glut2.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glut3.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glut4.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glut5.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glzn1.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glzn2.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glzn3.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glzn4.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/glzn5.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/gtpa1.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/gtpa2.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/gtpa3.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/gtpa4.map2/',
	'/Users/lysa8537/2018_RNASeq/accumulated_data/06_read_counting/stringtie_exgtf/gtpa5.map2/')

# select samples to compare from file table ##########

bg <- ballgown(samples, pData = filetable)

# To view the samples contained in the gown object (bg)
#sampleNames(bg)

# To access the exon, intron, transcript, or gene expression values use *expr(gown_object, <measure>)
# * = {e, i, t, g}
# <measure> = {'FPKM', 'cov', 'all', 'mcov', ...}
transcript_fpkm <- texpr(bg, 'FPKM')	#FPKM values for all transcripts
junction_rcount <- iexpr(bg)
gene_expression <- gexpr(bg)

# Visualizing the coverage for various isoforms of a gene
#library("ggplot2")
#dir <- '/scratch/Users/jast1849/sread2018/'
#gene <- 'MSTRG.6453'
#main <- gene
#samp <- sampleNames(bg)
#pdf(paste0(dir,"./transcripts_", gene, ".pdf"))
#plotTranscripts(gene, bg, samples=samp, meas='FPKM', colorby='transcript', main=main)
#dev.off()

# Visualizing mean abundances for given gene's isoforms between groups
#pdf(paste0(dir,"./means_", gene, ".pdf"))
#plotMeans(gene, bg, groupvar='condition', meas='FPKM', colorby='transcript')
#dev.off()


# DEA
DEres_gene <- stattest(bg, feature='gene', meas='FPKM', covariate='treatment')
DEres_transcript <- stattest(bg, feature='transcript', meas='FPKM', covariate='treatment')

# Filtering
#sigres_gene_0.1 <- subset(DEres_gene, qval <= 0.1)
#sigres_transcript_0.1 <- subset(DEres_transcript,qval <= 0.1)

#geneID <- as.vector(sigres_gene_0.1, mode="integer")
#transcriptID <- as.vector(sigres_transcript_0.1, mode = "integer")

#sigres_gene_names <- geneNames(bg)[c(geneID)]
#sigres_trans_names <- transcriptNames(bg)[c(transcriptID)]

# Outputting results files
write.table(DEres_gene,file="/Users/lysa8537/2018_RNASeq/accumulated_data/07_diff_ex/ballgown/all_samples_exgtf_gene_DE.txt",append=FALSE,sep='\t',quote=FALSE)
write.table(DEres_transcript,file="/Users/lysa8537/2018_RNASeq/accumulated_data/07_diff_ex/ballgown/all_samples_exgtf_transcript_DE.txt",append=FALSE,sep='\t',quote=FALSE)
