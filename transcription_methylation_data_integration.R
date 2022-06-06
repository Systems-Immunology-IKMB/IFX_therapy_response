#Code for correlation analysis between transcription and methylation--------

#Read and filter methylation data for sites
methylation_data <- read.csv("Downstream_analysis/IFX/site_meth_data_hg38.txt", header = TRUE, 
                             sep = '\t', row.names = 1)
meth_sample_info <- read.csv("EpicArray_IFX.csv", header = TRUE, sep = ',')
meth_sample_info$name <- paste("IFX", meth_sample_info$Patient, meth_sample_info$Timepoint, 
                               meth_sample_info$R.NR.C, sep = '_')
R_meth_sample_info <- subset(meth_sample_info, meth_sample_info$R.NR.C == 'R')
R_meth_sample_info$Sample_ID <- paste("X", R_meth_sample_info$Sample_ID, sep = '') 

R_meth_data <- methylation_data[, as.character(R_meth_sample_info$Sample_ID)]
colnames(R_meth_data) <- R_meth_sample_info$name

#Read and filter transcriptome data for genes
transcriptome_data <- read.csv("count_files/IFX_R_NR_sf_normalized_counts.txt", 
                               header = TRUE, sep = '\t', row.names = 1)
common_samples <- intersect(colnames(transcriptome_data), colnames(R_meth_data))
transcriptome_data <- transcriptome_data[, common_samples]
R_meth_data <- R_meth_data[, common_samples]

#Read list of differentially expressed genes 
deg_gene_list <- read.table("DESeq2_results/Timepoint_comparison/IFX/All_rem_only_genes.txt")
deg_gene_list <- deg_gene_list$V1

#Read and filter gene and linked methylation sites info, remove sites on sex chromosomes
gene_meth_sites <- read.csv("Downstream_analysis/IFX/gene_meth_sites_5000bp.txt", header = TRUE, 
                            sep = '\t')
gene_meth_sites <- subset(gene_meth_sites, gene_meth_sites$Chr != 'chrX' & gene_meth_sites$Chr != 'chrY' & 
                            gene_meth_sites$Chr != 'chrM')
rownames(gene_meth_sites) <- gene_meth_sites$Gene_id

deg_gene_list <- intersect(deg_gene_list, rownames(gene_meth_sites))

#Calculate correlation coefficient between the gene expression of each gene and the methylation intensity of each of its link methylated sites
#Calculate FDR using a permutation approach
deg_gene_meth_sites <- gene_meth_sites[as.character(deg_gene_list), ]
deg_gene_meth_sites <- subset(deg_gene_meth_sites, deg_gene_meth_sites$no_of_meth_sites > 0)

deg_gene_meth_sites$meth_sites <- as.character(deg_gene_meth_sites$meth_sites)
gene_meth_site_correlation <- matrix(nrow = sum(deg_gene_meth_sites$no_of_meth_sites), ncol = 5)
rownum = 1
for (i in 1:nrow(deg_gene_meth_sites)) {
  meth_sites <- strsplit(deg_gene_meth_sites[i,8], split = ';')
  for (j in 1:length(meth_sites[[1]])) {
    site_id_info <- strsplit(meth_sites[[1]][j], split = ':')
    site_id <- site_id_info[[1]][1]
    distance <- site_id_info[[1]][3]
    corr_data_frame <- data.frame(meth=t(R_meth_data[site_id,]), 
                                  expr=t(transcriptome_data[as.character(deg_gene_meth_sites[i,1]),]))
    corr_data_frame <- corr_data_frame[complete.cases(corr_data_frame),]
    rho <- cor(corr_data_frame[,1], corr_data_frame[,2], method="spearman")
    fdr <- 0
    for (k in 1:10000) {
      x <- sample(corr_data_frame[,1])
      y <- sample(corr_data_frame[,2])
      rand_rho <- cor(x, y, method="spearman")
      if(abs(rand_rho) >= abs(rho)){
        fdr <- fdr + 1
      }
    }
    fdr <- fdr/10000
    gene_meth_site_correlation[rownum, ] <- c(as.character(deg_gene_meth_sites[i,1]), site_id, distance, rho, fdr)
    rownum <- rownum + 1
  }
}
gene_meth_site_correlation <- as.data.frame(gene_meth_site_correlation)
colnames(gene_meth_site_correlation) <- c("Gene", "Site", "Distance from TSS", "Rho", "FDR")

write.table(gene_meth_site_correlation, file = "Downstream_analysis/IFX/R_DEG_meth_site_correlation_5000bp_10000rep.txt", quote = FALSE,
            sep = '\t')



#Code for analysis of DMPs correlated with DEGs------

#Get DMPs in remission samples
diffmeth_6w <- read.csv("Analysis/IFX/comparison_0h_6w/R/reports_differential/differential_methylation_data/diffMethTable_site_cmp1.csv", 
                        header = TRUE, sep = ',')
diffmeth_2w <- read.csv("Analysis/IFX/comparison_0h_2w/R/reports_differential/differential_methylation_data/diffMethTable_site_cmp1.csv", 
                        header = TRUE, sep = ',')

diffmeth_sites <- union(diffmeth_6w[diffmeth_6w$combinedRank <= 197485, ]$cgid, 
                        diffmeth_2w[diffmeth_2w$combinedRank <= 149541, ]$cgid)

#Filter for significant DMP-DEG correlations
sig_correlations <- subset(gene_meth_site_correlation, gene_meth_site_correlation$FDR < 0.01)


#Extract methylation data for DMPs significantly correlated with DEGs
meth_table2 <- R_meth_data[as.character(unique(sig_correlations$Site)), ]
meth_table2 <- meth_table2[rownames(meth_table2) %in% diffmeth_sites, ]

#Generate heatmap of methylation data for DMPs significantly correlated with DEGs for each sample as in Figure S4C and figure S4D
condition_colors <- unlist(lapply(colnames(meth_table2),function(x){
  if(grepl("12", x)) '#990000' #darkred
  else if(grepl("13", x)) '#999900' #redgreen
  else if(grepl('14',x)) '#009999' #greenblue
  else if(grepl('15',x)) '#4C0099' #purple
  else if(grepl('40',x)) '#990099' #magenta
  else if(grepl('42',x)) '#FFFF66' #yellow
  else if(grepl('43',x)) '#663300' #brown
}))

hypo_table <- meth_table2[rownames(hypo_meth_sites),]
hypo_distance <- as.dist(1-cor(t(hypo_table), method = "spearman"))
hypo_hclust <- hclust(hypo_distance)
hypo_table <- hypo_table[c(hypo_hclust$order), ]
hyper_table <- meth_table2[rownames(hyper_meth_sites),]
hyper_distance <- as.dist(1-cor(t(hyper_table), method = "spearman"))
hyper_hclust <- hclust(hyper_distance)
hyper_table <- hyper_table[c(hyper_hclust$order), ]

full_table <- rbind(hypo_table, hyper_table)

pdf("Downstream_analysis/IFX/R_DEG_correlated_DMPs_FDR_01_heatmap.pdf", height = 8, width = 10)
heatmap.2(as.matrix(full_table), trace="none", ColSideColors = condition_colors, Colv=NA,Rowv = NA,
          dendrogram = "none", density="none", scale="row", col=rev(brewer.pal(11, "RdYlBu")), lhei = c(1,5),
          colsep = c(3, 6, 9, 12, 15, 18), labCol = rep(c("0h", "2w", "6w"), 7), labRow = NA, margins = c(5,9))
dev.off()

