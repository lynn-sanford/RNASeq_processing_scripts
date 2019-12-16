library("ggplot2")
library("reshape2")
library("dplyr")
library("purrr")
library("cowplot")
library("extrafont")
library("ggpubr")
library("ggrepel")

# I usually ran this after DESeq2, so the TPM variable from that script was still in memory
# If running from an external text TPM file, you will need to load that into R

TPM_to_PCA <- t(TPM)

# Filter all NA and TPM 0 genes
TPM_to_PCA <- TPM_to_PCA[,apply(TPM_to_PCA, 2, var, na.rm=TRUE) != 0]

# Do principle component analysis
project.pca <- prcomp(TPM_to_PCA)

# These were also based on the variables already defined by the DESeq2 script
PCA_to_plot <- data.frame(project.pca$x,sample=colData_all$Treatment,rep=c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5,1,2,3,4,5))

# Separately plot all of the first 10 PCs

PC1_PC2_plot <- ggplot(PCA_to_plot,aes(x=PC1,y=PC2,col=sample)) +
  geom_point(size=2) +
  scale_color_manual(labels = c("Control","Glutamate","Glutamate/Zn","Glutamate/TPA"), values = c("cont"="turquoise4","glut"="darkorchid4","glzn"="orange3","gtpa"="red2")) +
  labs(x = "PC1", y = "PC2", title = "PC2 vs PC1") +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = rep), point.padding = .1,show.legend=FALSE)

# Output plot to file - only did for the PC1 vs. PC2 - all others just visualized in R
png(file="Z:/Courses_and_program_reqs/2018_11_Committee_meeting/PCA_glutamate_all_reps_PC1_PC2.png",width=1200,height=900,units="px",res = 260)
PC1_PC2_plot
dev.off()

PC3_PC4_plot <- ggplot(PCA_to_plot,aes(x=PC3,y=PC4,col=sample)) +
  geom_point(size=2) +
  scale_color_manual(labels = c("Control","Glutamate","Glutamate/Zn","Glutamate/TPA"), values = c("cont"="turquoise4","glut"="darkorchid4","glzn"="orange3","gtpa"="red2")) +
  labs(x = "PC3", y = "PC4", title = "PC4 vs PC3") +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = rep), point.padding = .1,show.legend=FALSE)

PC3_PC4_plot

PC5_PC6_plot <- ggplot(PCA_to_plot,aes(x=PC5,y=PC6,col=sample)) +
  geom_point(size=2) +
  scale_color_manual(labels = c("Control","Glutamate","Glutamate/Zn","Glutamate/TPA"), values = c("cont"="turquoise4","glut"="darkorchid4","glzn"="orange3","gtpa"="red2")) +
  labs(x = "PC5", y = "PC6", title = "PC5 vs PC6") +
  theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = rep), point.padding = .1,show.legend=FALSE)

PC5_PC6_plot

# To look at individual gene contributions - can output this to text file
loadings <- project.pca$rotation
