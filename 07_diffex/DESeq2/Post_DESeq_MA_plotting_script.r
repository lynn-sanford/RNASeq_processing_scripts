library(ggplot2)
library(ggrepel)

# This was based on variables already loaded into memory from the DESeq2 script
# Generally take the results variables from that script and put into plotting variables
# This is one plot example - just copy paste for each of the pairwise comparisons you need

glut_gtpa_to_plot <- data.frame(mean=log(res_glut_gtpa$baseMean,base=10),fc=res_glut_gtpa$log2FoldChange,p=res_glut_gtpa$padj,sig=(res_glut_gtpa$padj<0.05))

glut_gtpa_sig <- (res_glut_gtpa$padj<0.05 & complete.cases(res_glut_gtpa$padj))
glut_gtpa_sig_to_plot <- data.frame(mean = log(res_glut_gtpa[glut_gtpa_sig,]$baseMean,base=10), fc = res_glut_gtpa[glut_gtpa_sig,]$log2FoldChange, p = res_glut_gtpa[glut_gtpa_sig,]$padj)
row.names(glut_gtpa_sig_to_plot) <- row.names(res_glut_gtpa[glut_gtpa_sig,])
glut_gtpa_notsig_to_plot <- data.frame(mean = log(res_glut_gtpa[!glut_gtpa_sig,]$baseMean,base=10), fc = res_glut_gtpa[!glut_gtpa_sig,]$log2FoldChange, p = res_glut_gtpa[!glut_gtpa_sig,]$padj)
row.names(glut_gtpa_notsig_to_plot) <- row.names(res_glut_gtpa[!glut_gtpa_sig,])

glut_gtpa_05 <- ggplot() +
  geom_point(data = glut_gtpa_notsig_to_plot, aes(x = mean, y = fc), color = "gray60", size = .5) +
  geom_point(data = glut_gtpa_sig_to_plot, aes(x = mean, y = fc), color = "darkorange3", size = .7) +
  #geom_text_repel(data = glut_gtpa_sig_to_plot,aes(x = mean, y = fc,label = row.names(glut_gtpa_sig_to_plot)), point.padding = .1,show.legend=FALSE) +
  labs(x = "Log10(Mean read count)", y = "Log2(Fold change)", title = "Differential Expression Glutamate vs. Glutamate/TPA") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  
png(file="D:/Sequencing_data/2018_RNASeq/accum/diffex/map2_repfiltered/Glu_GTPA_p05_MAplot.png",width=1500,height=1000,units="px",res = 300)
glut_gtpa_05
dev.off()

