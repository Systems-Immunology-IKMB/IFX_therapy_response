library(DESeq2)
library(ggplot2)
library(magrittr)
library(WGCNA)
library(ComplexHeatmap)
library(reshape2)
library(RColorBrewer)
library(stringr)

#Code for constructing co-expression modules------

treatment <- "Infliximab"

#Read raw counts and sample info file
count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
col_data <- read.table("info_files/Sequenced_sample_information.csv", header = TRUE, sep = ',')

#Filter samples - Samples can be filtered for different timepoints and disease status groups
col_data <- subset(col_data, col_data$Treatment == treatment)
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))
col_data <- subset(col_data, col_data$Time.Point %in% c('0h'))

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point)
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
count_data <- count_data[, as.character(col_data$file_name)]

#Create DESeq object and normalize counts
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Remission.NR.C)
dds_counts <- estimateSizeFactors(dds_counts)
normalized_expression_measures <- counts(dds_counts, normalized=TRUE)

#Filter genes
sig_genes <- read.table("DESeq2_results/Timepoint_comparison/IFX/All_rem_only_genes.txt")
sig_genes <- sig_genes$V1

sig_genes_expression_data <- normalized_expression_measures[as.character(sig_genes),]
col_data$colnames <- paste("IFX", col_data$Patient.Nr., col_data$Time.Point, col_data$Remission.NR.C, sep = '_')
colnames(sig_genes_expression_data) <- col_data$colnames

#Log transform count data
log_counts <- log2(sig_genes_expression_data + 1)

#Define wgcna_matrix
wgcna_matrix <- t(log_counts)
s <- bicor(t(log_counts))

#Select parameters for scale free topology
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(t(log_counts), powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

#Create adjacency matrix
adj_mat <- adjacency.fromSimilarity(s, power=16, type = 'signed')
TOM = TOMsimilarity(adj_mat)
dissTOM = 1-TOM

#Make gene tree and define modules
gene_tree <- hclust(as.dist(1 - adj_mat), method = "average")

module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15, 
                                   deepSplit=TRUE)
module_colors <- labels2colors(module_labels)

plotDendroAndColors(gene_tree, module_colors, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')

mergedColor <- mergeCloseModules(t(log_counts),module_colors,cutHeight=.2)$color

#Output modules
gene_module <- data.frame(Gene = colnames(wgcna_matrix), ModuleLabel = module_labels, ModuleColor = mergedColor)
write.table(gene_module, "Coexpression_analysis/gene_module.txt", 
            quote = FALSE, sep = '\t', row.names = FALSE)

#Calculate module eigengenes
MEsO <- moduleEigengenes(wgcna_matrix, mergedColor)$eigengenes
MEs <- orderMEs(MEsO)
rownames(MEs) <- rownames(wgcna_matrix)

write.table(MEs, "Coexpression_analysis/module_eigengenes.txt", 
            quote = FALSE, sep = '\t')

MEs[, c("Treatment", "Patient", "Timepoint", "Status")] <- str_split_fixed(rownames(MEs), '_', 4)

#Code to create barplot as in Figure 2E-------

#Calculate number of genes per module
module_gene_counts <- as.data.frame(table(gene_module$ModuleColor))
module_gene_counts <- subset(module_gene_counts, module_gene_counts$Var1 != "grey")

module_color_names <- read.csv("Coexpression_analysis/module_color_names.txt", row.names = 1, 
                               sep = '\t')
module_gene_counts$Module <- module_color_names[as.character(module_gene_counts$Var1), "Module.name"]

pdf("Coexpression_analysis/module_gene_counts.pdf", 
    width = 6.5, height = 4)
p <- ggplot(module_gene_counts, aes(x=Module, y=Freq)) + geom_bar(stat = "identity", fill="grey", color="black")
p <- p + scale_x_discrete(limits=c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", 
                                   "M14", "M15", "M16", "M17", "M18", "M19", "M20", "M21", "M22", "M23", "M24"))
p <- p + ylab("#Genes")
p <- p + theme(axis.text=element_text(size=12, color = "black"), axis.title=element_text(size=16), legend.text = element_text(size=14), 
               legend.title = element_text(size=0), plot.title = element_text(size=16, face="bold", hjust = 0.5), axis.text.x = element_text(angle=45, hjust = 1)) +
  theme(axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), panel.background = element_blank())
p
dev.off()

#Code to create correlation heatmap as in Figure 2E and Figure S2C-----

#Read clinical data table
clinical_data <- read.csv("clinical_data/baseline_clinical_parameters.csv")
clinical_data <- clinical_data[rownames(MEs), ]

#Calculate correlation between module eigengenes and clinical parameters
module_trait_cor <- cor(MEs, clinical_data, use='p', method = "spearman")
module_trait_cor_pvalue <- corPvalueStudent(module_trait_cor, nrow(ME))
module_trait_cor_adj_pvalue <- p.adjust(module_trait_cor_pvalue, method = 'BH')
write.table(module_trait_cor, 
            file = "Coexpression_analysis/module_trait_correlation.txt", 
            sep = '\t', quote = FALSE)
write.table(module_trait_cor_pvalue, 
            file = "Coexpression_analysis/module_trait_correlation_P.txt", 
            sep = '\t', quote = FALSE)

module_trait_cor <- module_trait_cor[-25,]
module_trait_cor_pvalue <- module_trait_cor_pvalue[-25,]

rownames(module_trait_cor) <- gsub('ME', '', rownames(module_trait_cor))
rownames(module_trait_cor_pvalue) <- gsub('ME', '', rownames(module_trait_cor_pvalue))
module_trait_cor <- module_trait_cor[rownames(module_color_names), ]
module_trait_cor_pvalue <- module_trait_cor_pvalue[rownames(module_color_names), ]
rownames(module_trait_cor) <- as.character(module_color_names$Module.name)
rownames(module_trait_cor_pvalue) <- as.character(module_color_names$Module.name)

#Plot correlation heatmap

pdf(file = file.path("Coexpression_analysis/module_trait_correlation4.pdf"), 
    width=6, height = 5)
corrplot(t(module_trait_cor), method="color", type="full", col=viridis::viridis(30), p.mat = t(module_trait_cor_pvalue), sig.level = c(.001, .01, .05), pch.cex=1.5, pch.col="white", 
         tl.col = "black", insig="label_sig", tl.srt = 90, mar = c(0,0,1,0), cl.cex=1, cl.align.text = "l", cl.ratio = 0.3)
dev.off()

#Read cell fraction data obtained from CIBERSORTx
cell_fraction_data <- read.csv("CIBERSORT_analysis/IFX_blood/data/IFX_Blood_CIBERSORTx_Job25_Results.csv")
cell_fraction_data <- inner_join(col_data, cell_fraction_data, by = c("file_name"="Mixture"))
cell_fraction_data$Sample <- paste(cell_fraction_data$Patient.Nr., cell_fraction_data$Time.Point, sep = '_')
rownames(cell_fraction_data) <- cell_fraction_data$Sample

Cell_Name_vector <- c("B.cells.naive", "B.cells.memory",
                      "Plasma.cells", "T.cells.CD8", 
                      "T.cells.CD4.naive", 
                      "T.cells.CD4.memory.resting", # none measured 
                      "T.cells.CD4.memory.activated", 
                      "T.cells.follicular.helper", 
                      "T.cells.regulatory..Tregs.", 
                      "T.cells.gamma.delta", 
                      "NK.cells.resting",
                      "Monocytes", 
                      "Macrophages.M0", 
                      "Dendritic.cells.activated",
                      "Mast.cells.resting",
                      "Neutrophils")

cell_fraction_data <- cell_fraction_data[, Cell_Name_vector]
baseline_cell_fraction_data <- cell_fraction_data[rownames(MEs), ]

#Calculate correlation between module eigengenes and cell fraction data
module_cell_cor <- cor(MEs, baseline_cell_fraction_data, use='p', method = "spearman")
module_cell_cor_pvalue <- corPvalueStudent(module_cell_cor, nrow(ME))

module_cell_cor <- module_cell_cor[-25,]
module_cell_cor_pvalue <- module_cell_cor_pvalue[-25,]

module_cell_cor_adj_pvalue <- p.adjust(as.vector(module_cell_cor_pvalue), method = 'BH')
module_cell_cor_adj_pvalue <- matrix(module_cell_cor_adj_pvalue, nrow = nrow(module_cell_cor_pvalue), ncol = ncol(module_cell_cor_pvalue))
rownames(module_cell_cor_adj_pvalue) <- rownames(module_cell_cor_pvalue)
colnames(module_cell_cor_adj_pvalue) <- colnames(module_cell_cor_pvalue)

rownames(module_cell_cor) <- gsub('ME', '', rownames(module_cell_cor))
rownames(module_cell_cor_adj_pvalue) <- gsub('ME', '', rownames(module_cell_cor_adj_pvalue))
module_cell_cor <- module_cell_cor[rownames(module_color_names), ]
module_cell_cor_pvalue <- module_cell_cor_pvalue[rownames(module_color_names), ]
module_cell_cor_adj_pvalue <- module_cell_cor_adj_pvalue[rownames(module_color_names), ]
rownames(module_cell_cor) <- as.character(module_color_names$Module.name)
rownames(module_cell_cor_pvalue) <- as.character(module_color_names$Module.name)
rownames(module_cell_cor_adj_pvalue) <- as.character(module_color_names$Module.name)


#Plot correlation heatmap

pdf(file = file.path("DESeq2_results/Timepoint_comparison/IFX/Coexpression_analysis_R_only_genes_separate_timepoints/Coexpression_analysis_0h/module_cell_fraction_correlation.pdf"), 
    width=10, height = 12)
corrplot(module_cell_cor, method="color", type="full", col=viridis::viridis(30), p.mat = module_cell_cor_adj_pvalue, sig.level = c(.001, .01, .05), pch.cex=1.5, pch.col="white", 
         tl.col = "black", insig="label_sig", tl.srt = 45, mar = c(0,0,1,0), cl.cex=1.5, cl.align.text = "l", cl.ratio = 0.3, tl.cex = 1.5)
dev.off()




#Code for module preservation analysis-----

#all_colors is a merged module data for baseline, week 2 and week 6 - this can be created for remission and non-remission samples
all_colors <- data.frame(baseline=baseline_coexpression_modules$ModuleColor, early=early_coexpression_modules$ModuleColor, 
                         late=late_coexpression_modules$ModuleColor)

#Read raw counts and sample info file

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
col_data <- read.table("info_files/Sequenced_sample_information.csv", header = TRUE, sep = ',')

col_data <- subset(col_data, col_data$Treatment == "Infliximab")
col_data <- subset(col_data, col_data$Remission.NR.C %in% c('R', 'NR'))
col_data <- subset(col_data, col_data$Time.Point %in% c('0h', '2w', '6w'))

col_data$Patient.Nr. <- factor(col_data$Patient.Nr.)
col_data$Time.Point <- factor(col_data$Time.Point)
col_data$Remission.NR.C <- factor(col_data$Remission.NR.C)
col_data$colnames <- paste("IFX", col_data$Patient.Nr., col_data$Time.Point, col_data$Remission.NR.C, sep = '_')
count_data <- count_data[, as.character(col_data$file_name)]

#Create DESeq object and normalize counts

dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Remission.NR.C + Time.Point)
dds_counts <- estimateSizeFactors(dds_counts)
normalized_expression_measures <- counts(dds_counts, normalized=TRUE)
expr_dat <- normalized_expression_measures[as.character(all_colors$Gene), ]
colnames(expr_dat) <- col_data$colnames

#Module preservation in remission/non-remission samples

expr_dat <- t(log(expr_dat+1))

multiExpr <- list()
multiExpr[[1]] = list(data = expr_dat[grep("0h", rownames(expr_dat)), ])
multiExpr[[2]] = list(data = expr_dat[grep("2w_R", rownames(expr_dat)), ])
multiExpr[[3]] = list(data = expr_dat[grep("6w_R", rownames(expr_dat)), ])

setLabels <- c("baseline", "early", "late")
names(multiExpr) <- setLabels

colorList <- list(all_colors$baseline[, drop=TRUE], all_colors$early[, drop=TRUE], all_colors$late[, drop=TRUE])

names(colorList) = setLabels

nSets = 3

system.time( {
  mp = modulePreservation(multiExpr, colorList,
                          referenceNetworks = c(1),
                          loadPermutedStatistics = FALSE,
                          nPermutations = 200,
                          verbose = 3)
} )
# Save the results
save(mp, file = "Coexpression_analysis/Module_preservation/baseline_early_late_module_R_Preservation.RData")




#Code to plot Zsummary heatmap as in Figure 2F and Figure 5D-----

#Load preservation data for remission and non-remission samples
load(file = "Coexpression_analysis/Module_preservation/baseline_early_late_R_modulePreservation.RData")
R_mp <- mp
load(file = "Coexpression_analysis/Module_preservation/baseline_early_late_NR_modulePreservation.RData")
NR_mp <- mp

#Extract Zsummary scores from the preservation objects
early_R_statsZ = cbind(R_mp$quality$Z[[1]][[2]][, -1], R_mp$preservation$Z[[1]][[2]][, -1])
late_R_statsZ = cbind(R_mp$quality$Z[[1]][[3]][, -1], R_mp$preservation$Z[[1]][[3]][, -1])
early_NR_statsZ = cbind(NR_mp$quality$Z[[1]][[2]][, -1], NR_mp$preservation$Z[[1]][[2]][, -1])
late_NR_statsZ = cbind(NR_mp$quality$Z[[1]][[3]][, -1], NR_mp$preservation$Z[[1]][[3]][, -1])

Zsummary_data_frame <- data.frame(module_color=rownames(early_R_statsZ), early_R_statsZ=early_R_statsZ$Zsummary.pres, 
                                  late_R_statsZ=late_R_statsZ$Zsummary.pres, early_NR_statsZ=early_NR_statsZ$Zsummary.pres, 
                                  late_NR_statsZ=late_NR_statsZ$Zsummary.pres)
Zsummary_data_frame <- subset(Zsummary_data_frame, !(Zsummary_data_frame$module_color %in% c("gold", "grey")))

Zsummary_data_frame <- inner_join(module_color_name, Zsummary_data_frame, by=c("Module.color"="module_color"))

rownames(Zsummary_data_frame) <- Zsummary_data_frame$Module.name
module_colors <- as.character(Zsummary_data_frame$Module.color)
Zsummary_data_frame[, c("Module.color", "Module.name")] <- NULL
condition_colors <- c(rep("#00CC00", 2), rep("#0066CC", 2))

#Plot Zsummary heatmap

Zsummary_data_frame <- Zsummary_data_frame[, c("early_NR_statsZ", "early_R_statsZ", "late_NR_statsZ", "late_R_statsZ")]

col_ha = rowAnnotation(df = data.frame(Status=c("NR", "R", "NR", "R")),
                       col = list(Status=c('R'="#00CC00", 'NR'="#3399FF")),
                       gp = gpar(col = "black"))

pdf("Coexpression_analysis/Module_preservation/baseline_early_late_R_NR_Zsummary_heatmap.pdf", height=3, width=6.5)
Heatmap(scale(t(Zsummary_data_frame)), col=brewer.pal(9, "Blues"), left_annotation = col_ha, cluster_rows = FALSE, cluster_columns = FALSE, 
        name = "Scaled\nZsummary\nscore", row_labels = c("2w NR", "2w R", "6w NR", "6w R"))
dev.off()
#Code to plot Zsummary dotplot as in Figure S2D and Figure S6D------

Zsummary_all_data <- melt(Zsummary_data_frame, id.vars = c("Module.color", "Module.name"))
Zsummary_all_data[, c("Timepoint", "Status")] <- str_split_fixed(Zsummary_all_data$variable, '_', 3)[, 1:2]

pdf("Coexpression_analysis/Module_preservation/baseline_early_late_R_NR_differentially_preserved_modules_Zsummary.pdf")
p <- ggplot(data=Zsummary_all_data, aes(x=Module.name, y=value, color=Status, shape=Timepoint))
p <- p + geom_point(size=3)
p <- p + geom_hline(yintercept = 2, color="blue", lty=2) + geom_hline(yintercept = 10, color="green", lty=2)
p <- p + xlab("Module") + ylab("Zsummary")
p <- p + scale_color_manual(values = c("#0066CC", "#00CC00"))
p <- p + scale_x_discrete(limits = unique(Zsummary_all_data$Module.name))
p <- p + theme_bw() + theme(axis.text=element_text(color="black", size=14), axis.title=element_text(size=16), 
                            legend.text = element_text(size=12), legend.title = element_text(size=0), axis.text.x = element_text(angle = 45, hjust = 1),
                            plot.title = element_text(size=16, face="bold", hjust = 0.5), 
                            panel.border = element_blank(), axis.line = element_line(size=0.4, color = "black"))
p
dev.off()

#Plot kME heatmaps for differentially preserved modules as in Figures S6E and S6F---------

multiExpr <- list()
multiExpr[[1]] = list(data = expr_dat[grep("0h", rownames(expr_dat)), ])
multiExpr[[2]] = list(data = expr_dat[grep("2w_R", rownames(expr_dat)), ])
multiExpr[[3]] = list(data = expr_dat[grep("6w_R", rownames(expr_dat)), ])
multiExpr[[4]] = list(data = expr_dat[grep("2w_NR", rownames(expr_dat)), ])
multiExpr[[5]] = list(data = expr_dat[grep("6w_NR", rownames(expr_dat)), ])

setLabels <- c("baseline", "early_R", "late_R", "early_NR", "late_NR")
names(multiExpr) <- setLabels

nSets = 5

# Get eigengenes
mes = list()
for (set in 1:nSets)
{
  mes[[set]] = moduleEigengenes(multiExpr[[set]]$data, baseline_coexpression_modules$ModuleColor)$eigengenes
}

kME <- list()
for (set in 1:nSets)
{
  kME[[set]] = signedKME(multiExpr[[set]]$data, mes[[set]])
}

#Extract kME for each module and plot heatmap

all_mes <- moduleEigengenes(expr_dat, baseline_coexpression_modules$ModuleColor)$eigengenes

plum_genes <- as.character(subset(baseline_coexpression_modules, baseline_coexpression_modules$ModuleColor == "plum")$Gene)

plum_kme <- data.frame(baseline=kME[[1]]$kMEplum, early_NR=kME[[4]]$kMEplum, early_R=kME[[2]]$kMEplum, late_NR=kME[[5]]$kMEplum, late_R=kME[[3]]$kMEplum)
rownames(plum_kme) <- rownames(kME[[1]])
plum_kme <- plum_kme[plum_genes, ]
rownames(plum_kme) <- gene_info_file[rownames(plum_kme), "Gene_name"]

pdf("DESeq2_results/Timepoint_comparison/IFX/Coexpression_analysis_R_only_genes_separate_timepoints/Module_preservation/plum_kME2.pdf", width = 4, height = 12)
Heatmap(t(scale(t(plum_kme))), col = brewer.pal(9, "Greens"), cluster_columns = FALSE, row_names_gp = gpar(fontsize=7), name = "kME")
dev.off()


