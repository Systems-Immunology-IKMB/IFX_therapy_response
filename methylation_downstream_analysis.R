library(ggplot2)
library(reshape2)
library(plyr)

#Code to get number of differentially methylated sites (DMPs) and differentially methylated regions (DMRs) per timepoint and disease status defined by the automatically selected rank cutoff--------

#Extract number of differentially methylated sites for remission samples
dmp_6w_R <- read.csv("Analysis/IFX/comparison_0h_6w/R/reports_differential/differential_methylation_data/diffMethTable_site_cmp1.csv", 
                   header = TRUE, sep = ',')
dmp_2w_R <- read.csv("Analysis/IFX/comparison_0h_2w/R/reports_differential/differential_methylation_data/diffMethTable_site_cmp1.csv", 
                   header = TRUE, sep = ',')
dmp_6w_R <- dmp_6w_R[dmp_6w_R$combinedRank <= 197485, ]
dmp_2w_R <- dmp_2w_R[dmp_2w_R$combinedRank <= 149541, ]
uplen_6w <- nrow(subset(dmp_6w_R, dmp_6w_R$mean.X0h <= dmp_6w_R$mean.X6w))
uplen_2w <- nrow(subset(dmp_2w_R, dmp_2w_R$mean.X0h <= dmp_2w_R$mean.X2w))
downlen_6w <- nrow(subset(dmp_6w_R, dmp_6w_R$mean.X0h > dmp_6w_R$mean.X6w))
downlen_2w <- nrow(subset(dmp_2w_R, dmp_2w_R$mean.X0h > dmp_2w_R$mean.X2w))

dmp_data_R <- data.frame(tp=c("2w", "6w", "2w", "6w"), 
                         direction=c("up", "up", "down", "down"),
                         number=c(uplen_2w, uplen_6w, downlen_2w, downlen_6w), 
                         status=rep('R', 4))

#Extract number of differentially methylated sites for non-remission samples
dmp_6w_NR <- read.csv("Analysis/IFX/comparison_0h_6w/NR/reports_differential/differential_methylation_data/diffMethTable_site_cmp1.csv", 
                     header = TRUE, sep = ',')
dmp_2w_NR <- read.csv("Analysis/IFX/comparison_0h_2w/NR/reports_differential/differential_methylation_data/diffMethTable_site_cmp1.csv", 
                     header = TRUE, sep = ',')
dmp_6w_NR <- dmp_6w_NR[dmp_6w_NR$combinedRank <= 12202, ]
dmp_2w_NR <- dmp_2w_NR[dmp_2w_NR$combinedRank <= 190429, ]
uplen_6w <- nrow(subset(dmp_6w_NR, dmp_6w_NR$mean.X0h <= dmp_6w_NR$mean.X6w))
uplen_2w <- nrow(subset(dmp_2w_NR, dmp_2w_NR$mean.X0h <= dmp_2w_NR$mean.X2w))
downlen_6w <- nrow(subset(dmp_6w_NR, dmp_6w_NR$mean.X0h > dmp_6w_NR$mean.X6w))
downlen_2w <- nrow(subset(dmp_2w_NR, dmp_2w_NR$mean.X0h > dmp_2w_NR$mean.X2w))

dmp_data_NR <- data.frame(tp=c("2w", "6w", "2w", "6w"), 
                         direction=c("up", "up", "down", "down"),
                         number=c(uplen_2w, uplen_6w, downlen_2w, downlen_6w), 
                         status=rep('NR', 4))

dmp_data <- rbind(dmp_data_R, dmp_data_NR)
dmp_data$number <- ifelse(dmp_data$direction == "up", dmp_data$number, -1*dmp_data$number)

#Plot bar plot as in Figure 4B
pdf(file = "Downstream_analysis/IFX/dmp_number_comparison_sites.pdf", width = 4, height = 4)
p <- ggplot(data = dmp_data, mapping = aes(x=tp, y=number, color=status, fill=status, alpha=direction, group=status)) + 
  geom_bar(stat = "identity", position = "dodge", color="black")
p <- p + scale_alpha_manual(values=c(0.5, 1)) + scale_fill_manual(values=c("#00CC00", "#3399FF"))
p <- p + xlab("Timepoints") + ylab("Number of sites") 
p <- p + theme_minimal() + theme(axis.text=element_text(size=14, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=12), 
                                 legend.title = element_text(size=0), plot.title = element_text(size=20, face="bold", hjust = 0.5))
p
dev.off()

#Get overlapping DMPs between remission and non-remission samples

R_dmps <- union(dmp_2w_R$cgid, dmp_6w_R$cgid)
NR_dmps <- union(dmp_2w_NR$cgid, dmp_6w_NR$cgid)

library(VennDiagram)

#Plot venn diagram as in Figure 4C
pdf("Downstream_analysis/IFX/DMP_sites_R_NR_venn.pdf")
grid.newpage()
draw.pairwise.venn(area1 = length(R_dmps), area2 = length(NR_dmps), cross.area = length(intersect(R_dmps, NR_dmps)), 
                   category = c("R", "NR"), col = c("#00CC00", "#3399FF"), 
                   fill = c("white", "white"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.025, 0.04), 
                   cat.cex = rep(1.7, 2), 
                   cat.fontfamily = rep("Helvetica", 2), cex = rep(1.5, 3), fontfamily = rep("Helvetica", 3))
dev.off()


#Extract number of differentially methylated regions for remission/non-remission samples
dmpr_6w <- read.csv("Analysis/IFX/comparison_0h_6w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv", 
                    header = TRUE, sep = ',')
dmpr_2w <- read.csv("Analysis/IFX/comparison_0h_2w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_promoters.csv", 
                    header = TRUE, sep = ',')
dmpr_6w <- dmpr_6w[dmpr_6w$combinedRank <= 161, ]
dmpr_2w <- dmpr_2w[dmpr_2w$combinedRank <= 100, ]
num_dmpr_6w <- nrow(dmpr_6w)
num_dmpr_2w <- nrow(dmpr_2w)


dmg_6w <- read.csv("Analysis/IFX/comparison_0h_6w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_genes.csv", 
                   header = TRUE, sep = ',')
dmg_2w <- read.csv("Analysis/IFX/comparison_0h_2w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_genes.csv", 
                   header = TRUE, sep = ',')
dmg_6w <- dmg_6w[dmg_6w$combinedRank <= 98, ]
dmg_2w <- dmg_2w[dmg_2w$combinedRank <= 80, ]
num_dmg_6w <- nrow(dmg_6w)
num_dmg_2w <- nrow(dmg_2w)


dmcpg_6w <- read.csv("Analysis/IFX/comparison_0h_6w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_cpgislands.csv", 
                     header = TRUE, sep = ',')
dmcpg_2w <- read.csv("Analysis/IFX/comparison_0h_2w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_cpgislands.csv", 
                     header = TRUE, sep = ',')
dmcpg_6w <- dmcpg_6w[dmcpg_6w$combinedRank <= 76, ]
dmcpg_2w <- dmcpg_2w[dmcpg_2w$combinedRank <= 0, ]
num_dmcpg_6w <- nrow(dmcpg_6w)
num_dmcpg_2w <- nrow(dmcpg_2w)


dmde_6w <- read.csv("Analysis/IFX/comparison_0h_6w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_ensembleRegBuildBPdistal.csv", 
                    header = TRUE, sep = ',')
dmde_2w <- read.csv("Analysis/IFX/comparison_0h_2w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_ensembleRegBuildBPdistal.csv", 
                    header = TRUE, sep = ',')
dmde_6w <- subset(dmde_6w, dmde_6w$combinedRank <= 623)
dmde_2w <- subset(dmde_2w, dmde_2w$combinedRank <= 1659)
num_dmde_6w <- nrow(dmde_6w)
num_dmde_2w <- nrow(dmde_2w)


dmpe_6w <- read.csv("Analysis/IFX/comparison_0h_6w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_ensembleRegBuildBPproximal.csv", 
                    header = TRUE, sep = ',')
dmpe_2w <- read.csv("Analysis/IFX/comparison_0h_2w/R/reports_differential/differential_methylation_data/diffMethTable_region_cmp1_ensembleRegBuildBPproximal.csv", 
                    header = TRUE, sep = ',')
dmpe_6w <- subset(dmpe_6w, dmpe_6w$combinedRank <= 935)
dmpe_2w <- subset(dmpe_2w, dmpe_2w$combinedRank <= 1648)
num_dmpe_6w <- nrow(dmpe_6w)
num_dmpe_2w <- nrow(dmpe_2w)


freq_data <- data.frame(Annotation=c("Promoters", "Genes", "CpG", "Distal Enhancer", "Proximal Enhancer"), 
                        Week2=c(num_dmpr_2w, num_dmg_2w, num_dmcpg_2w, num_dmde_2w, num_dmpe_2w), 
                        Week6=c(num_dmpr_6w, num_dmg_6w, num_dmcpg_6w, num_dmde_6w, num_dmpe_6w))
library(reshape2)
freq_data <- melt(freq_data)
freq_data[freq_data$value == 0,]$value = 1

#Plot bar plot as in Figure S3D and S3E 
pdf("Analysis/IFX/R_DMR_number_comparison.pdf")
p <- ggplot(data = freq_data, mapping = aes(x=Annotation, y=value, fill=variable, group=variable)) + 
  geom_bar(stat = "identity", position = "dodge", color="black")
p <- p + scale_fill_manual(values=c("#99CCFF", "#000099"))
p <- p + xlab("Annotation") + ylab("Number of regions") 
p <- p + scale_x_discrete(limits = c("Promoters", "Genes", "CpG", "Proximal Enhancer", "Distal Enhancer"))
p <- p + theme_bw() + 
  theme(axis.text=element_text(size=14, color="black"), axis.title=element_text(size=16), legend.text = element_text(size=12), 
        legend.title = element_text(size=0), plot.title = element_text(size=20, face="bold", hjust = 0.5), 
        axis.text.x = element_text(angle=45, hjust=1))
p
dev.off()

