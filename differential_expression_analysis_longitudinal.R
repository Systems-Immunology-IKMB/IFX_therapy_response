library(ImpulseDE2)
library(DESeq2)
library(magrittr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(compiler)
#Code for pairwise differential expression analysis---------

#create result directory if it does not already exists
if (!dir.exists("DESeq2_results/")){
  dir.create("DESeq2_results/")
}

#Intialize parameters
#Parameters can be modified for different timepoints, disease status groups and treatments
condition1 <- '0h'
condition2 <- '4h'
condition_column_name <- 'Time.Point'
patient_column_name <- 'Patient.Nr.'
file_column_name <- 'file_name'
treatment <- 'Infliximab'
filtering_options <- data.frame(column_name=c('Treatment', 'Remission.NR.C'), values=c('Infliximab', 'R'))

#Read raw counts and sample info file

count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
col_data <- read.table("info_files/Sequenced_sample_information.csv", header = TRUE, sep = ',')

#Create output folders

analysis_folder <- file.path("DESeq2_results/", paste("0-", condition2, "_comparison/", sep = ''))
if (!dir.exists(analysis_folder)) {
  dir.create(analysis_folder)
}
output_folder <-file.path(analysis_folder, ifelse(filtering_options$values[1] == '', '',
                             (ifelse(filtering_options$values[1] == 'Infliximab', 'IFX/', 'VDZ/'))))

if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

#Get patients for conditions
pat_use1 <- subset(col_data, col_data[condition_column_name] == condition1)[patient_column_name]
pat_use2 <- subset(col_data, col_data[condition_column_name] == condition2)[patient_column_name]
pat_use <- intersect(pat_use1[,1], pat_use2[,1])

#Create sample names
for(i in 1:length(col_data[,condition_column_name])){
  trt <- ifelse(col_data[i, as.vector(filtering_options$column_name[1])]=='Infliximab', 'IFX', 
                ifelse(col_data[i, as.vector(filtering_options$column_name[1])]=='Vedolizumab', 'VDZ', 'TZZ'))
  col_data$colnames[i] <- paste(col_data[i, condition_column_name], trt, 
                                col_data[i, patient_column_name], sep = "_")
}

#Filter data
col_data <- subset(col_data, col_data[,patient_column_name] %in% pat_use)
col_data <- subset(col_data, col_data[condition_column_name] == condition1 | col_data[condition_column_name] == condition2)
for(i in 1:length(filtering_options$column_name)){
  filter_column <- as.vector(filtering_options$column_name[i])
  value <- filtering_options$values[i]
  if(value == '' & filter_column == 'Treatment'){
    value <- treatment
  }
  col_data <- subset(col_data, col_data[,filter_column] %in% value)
}

col_data[,patient_column_name] <- factor(as.character(col_data[,patient_column_name]))
col_data[,condition_column_name] <- factor(col_data[,condition_column_name])
count_data <- count_data[, as.character(col_data[,file_column_name])]

#Create DESeq object

dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Patient.Nr. + Time.Point)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

#Run DESeq
dds <- DESeq(dds_counts, betaPrior = FALSE)
res <- results(dds, independentFiltering = TRUE, alpha = 0.05)

#Output results to file
sink(file = paste(output_folder, "DESeq2result_summary_", condition1, "vs", condition2, "_", filtering_options$values[2], ".txt", sep = ""))
cat(summary(res))
sink()
png(file=paste(output_folder, "MAplot_", condition1, "_vs_", condition2, "_", filtering_options$values[2], ".png", sep = ""), width = 500, height = 500)
plotMA(res, alpha=0.05, main=paste(condition1, "vs.", condition2, filtering_options$values[1], filtering_options$values[2], sep = " "), ylim=c(-6,6))
dev.off()
res_sorted <- res[order(res$padj), ]
write.table(res_sorted, file = paste(output_folder, "DESeq2result_", condition1, "vs", condition2, "_", filtering_options$values[2], ".txt", sep = ""), 
            sep="\t", quote=FALSE)

#Add gene symbols
mart_export <- read.table("../reference/human/GRCh38_mart_export.txt", header = T, sep="\t")
unique_mart <- subset(mart_export, duplicated(mart_export$Gene_stable_ID) == FALSE)
rownames(unique_mart) <- unique_mart$Gene_stable_ID
res_sorted$gene <- unique_mart[rownames(res_sorted), ]$Gene_name
write.table(res_sorted, file = paste(output_folder, "DESeq2result_genenames_", condition1, "vs", condition2, "_", filtering_options$values[2], ".txt", sep = ""), 
            sep="\t", quote=FALSE)


#Code for longitudinal differential expression analysis---------


#Read raw counts and sample info file
count_data <- read.table("count_files/all_merged_gene_counts.txt", header = TRUE, sep = '\t', row.names = 1)
col_data <- read.table("info_files/Sequenced_sample_information.csv", header = TRUE, sep = ',')


#Filter data - data can be filtered for each disease status group
col_data$Patient.Nr. <- factor(as.character(col_data$Patient.Nr.))
col_data$Time.Point <- factor(col_data$Time.Point)
col_data$sampleName <- paste(col_data$Time.Point, col_data$Patient.Nr., sep = "_")
col_data <- subset(col_data, col_data$Treatment == "Infliximab")
col_data <- subset(col_data, col_data$Remission.NR.C == 'R')

count_data <- count_data[, as.character(col_data$file_name)]
colnames(count_data) <- col_data$sampleName

#Create DESeq object
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Patient.Nr. + Time.Point)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]

#Create annotation data for ImpulseDE2
matCountData <- counts(dds_counts)
dfAnnotation <- data.frame(Sample=col_data$sampleName, Condition=rep('case', nrow(col_data)), 
                           TimeCateg=col_data$Time.Point, Batch=col_data$Patient.Nr.)
dfAnnotation$Time <- ifelse(dfAnnotation$TimeCateg=='0h', 1, ifelse(dfAnnotation$TimeCateg=='4h', 2, 
                                                                    ifelse(dfAnnotation$TimeCateg=='24h', 3, 
                                                                           ifelse(dfAnnotation$TimeCateg=='72h', 4, 
                                                                                  ifelse(dfAnnotation$TimeCateg=='2w', 5,
                                                                                         ifelse(dfAnnotation$TimeCateg=='6w', 6, 7))))))

dfAnnotation <- dfAnnotation[order(dfAnnotation$Time),]
matCountData <- matCountData[, as.character(dfAnnotation$Sample)]

#Run ImpulseDE2
objectImpulseDE2 <- runImpulseDE2(
  matCountData           = matCountData, 
  dfAnnotation           = dfAnnotation,
  boolCaseCtrl           = FALSE,
  vecConfounders         = c("Batch"),
  boolIdentifyTransients = TRUE,
  scaNProc               = 2 )

res_sorted2 <- objectImpulseDE2$dfImpulseDE2Results[order(objectImpulseDE2$dfImpulseDE2Results$padj),]

save.image("Longitudinal_analysis/IFX/impulseDE2_R.RData")
write.table(res_sorted2, file = "Longitudinal_analysis/IFX/ImpulseDE2_result_transient_batch_case_R.txt", sep="\t", quote=FALSE)




#Code to generate barplot as in Figure 2B--------

#Get number of DEGs at each timepoint for remission patients
working_folder <- "DESeq2_results/"
timepoint <- c("4h", "24h", "72h", "2w", "6w", "14w")
R_pairwise_degs <- c()
deg_data_R <- data.frame(timepoint=as.character(), upgenes=as.integer(), downgenes=as.integer())
for (tp in timepoint) {
  file <- paste(working_folder, "0-", tp, "_comparison/IFX/DESeq2result_0hvs", tp, "_R.txt", sep = "")
  data <- read.table(file = file, header = TRUE, sep = "\t")
  R_pairwise_degs <- union(R_pairwise_degs, rownames(subset(data, data$padj < 0.05)))
  uplen <- nrow(subset(data, data$padj < 0.05 & data$log2FoldChange > 0))
  downlen <- nrow(subset(data, data$padj < 0.05 & data$log2FoldChange < 0))
  deg_data_R <- rbind(deg_data_R, data.frame(tp, uplen, downlen))
}


#Get number DEGs at each timepoint for non remission patients
working_folder <- "DESeq2_results/"
timepoint <- c("4h", "24h", "72h", "2w", "6w", "14w")
NR_pairwise_degs <- c()
deg_data_NR <- data.frame(timepoint=as.character(), upgenes=as.integer(), downgenes=as.integer())
for (tp in timepoint) {
  file <- paste(working_folder, "0-", tp, "_comparison/IFX/DESeq2result_0hvs", tp, "_NR.txt", sep = "")
  data <- read.table(file = file, header = TRUE, sep = "\t")
  NR_pairwise_degs <- union(NR_pairwise_degs, rownames(subset(data, data$padj < 0.05)))
  uplen <- nrow(subset(data, data$padj < 0.05 & data$log2FoldChange > 0))
  downlen <- nrow(subset(data, data$padj < 0.05 & data$log2FoldChange < 0))
  deg_data_NR <- rbind(deg_data_NR, data.frame(tp, uplen, downlen))
}

library(reshape2)
#Reformat data

deg_data_R <- melt(deg_data_R)
deg_data_R$status <- "R"
deg_data_NR <- melt(deg_data_NR)
deg_data_NR$status <- "NR"
deg_data <- rbind(deg_data_R, deg_data_NR)

index <- deg_data$variable == "downlen"
deg_data$value[index] <- -1*deg_data$value[index]

#Get number of longitudinal DEGs for remission and non-remission patients
impulse_data_R <- read.table("Longitudinal_analysis/IFX/ImpulseDE2_result_transient_batch_case_R.txt", 
                             header = TRUE, sep = '\t')
impulse_data_R <- subset(impulse_data_R, impulse_data_R$padj < 0.05)
impulse_data_NR <- read.table("Longitudinal_analysis/IFX/ImpulseDE2_result_transient_batch_case_NR.txt", 
                              header = TRUE, sep = '\t')
impulse_data_NR <- subset(impulse_data_NR, impulse_data_NR$padj < 0.05)
impulse_plot_data <- data.frame(Status=c("R", "NR"), 
                                Number=c(nrow(impulse_data_R), nrow(impulse_data_NR)))

deg_data <- rbind(deg_data, data.frame(tp="Longitudinal", variable="impulse", value=impulse_plot_data$Number, status=impulse_plot_data$Status))

#Plot bar plot
library(ggplot2)

pdf(file = "DESeq2_results/Timepoint_comparison/IFX/deg_impulse_number_comparison.pdf", width = 9, height = 7)
p <- ggplot(data = deg_data, mapping = aes(x=tp, y=value, color=status, fill=status, alpha=variable, group=status)) + 
  geom_bar(stat = "identity", position = "dodge", color="black")
p <- p + scale_alpha_manual(values=c(1, 0.5, 1)) + scale_fill_manual(values=c("#3399FF", "#00CC00"))
p <- p + xlab("Time points") + ylab("Number of genes") + 
  scale_y_continuous(breaks = round_any(seq(min(deg_data$value), max(deg_data$value), by=1000), 100))
p <- p + theme_minimal() + theme(axis.text=element_text(size=14), axis.title=element_text(size=16), legend.text = element_text(size=12), 
                                 legend.title = element_text(size=0), plot.title = element_text(size=20, face="bold", hjust = 0.5))
p
dev.off()
#Code to generate venn diagram as in Figure 2C-----

R_degs <- union(R_pairwise_degs, rownames(impulse_data_R))
NR_degs <- union(NR_pairwise_degs, rownames(impulse_data_NR))

library(VennDiagram)
pdf(paste("Timepoint_comparison/IFX/pw_and_impulse_R_NR.pdf", sep = ""))
grid.newpage()
draw.pairwise.venn(area1 = length(R_degs), area2 = length(NR_degs), cross.area = length(intersect(R_degs, NR_degs)), 
                   category = c("R", "NR"), col = c("#00CC00", "#0066CC"), 
                   fill = c("white", "white"), alpha = rep(0.5, 2), cat.pos = c(0, 0), cat.dist = c(0.025, 0.04), 
                   cat.cex = rep(1.7, 2), 
                   cat.fontfamily = rep("Arial", 2), cex = rep(1.5, 3), fontfamily = rep("Arial", 3))
dev.off()




#Code to generate heatmap as in Figure 2D-----

#Create DESeq object and normalize counts
dds_counts <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Patient.Nr.+ Time.Point)
dds_counts <- dds_counts[ rowSums(counts(dds_counts)) > 1, ]
dds_counts <- estimateSizeFactors(dds_counts)

#Select top 10 DEGs for each timepoint in remission patients
upgenes <- c()
downgenes <- c()

working_folder <- "DESeq2_results/"
timepoint <- c("4h", "24h", "72h", "2w", "6w", "14w")
for (tp in timepoint) {
  file <- paste(working_folder, "0-", tp, "_comparison/IFX/DESeq2result_0hvs", tp, "_R.txt", sep = "")
  data <- read.table(file = file, header = TRUE, sep = "\t")
  up <- rownames(subset(data, data$padj < 0.05 & data$log2FoldChange > 0))
  down <- rownames(subset(data, data$padj < 0.05 & data$log2FoldChange < 0))
  if (length(up) > 10){
    up <- up[1:10]
  }
  if (length(down) > 10){
    down <- down[1:10]
  }
  upgenes <- c(upgenes, up)
  downgenes <- c(downgenes, down)
}

#Select top 30 longitudinal DEGs in remission patients
file <- "Longitudinal_analysis/IFX/ImpulseDE2_result_transient_batch_case_R.txt"
data <- read.table(file = file, header = TRUE, sep = "\t")
impulse_deg <- rownames(subset(data, data$padj < 0.05))
if (length(impulse_deg) > 30){
  impulse_deg <- impulse_deg[1:30]
}

#Combine all DEGs
upgenes <- upgenes[!duplicated(upgenes)]
downgenes <- downgenes[!duplicated(downgenes)]
impulse_deg <- impulse_deg[!(impulse_deg %in% upgenes)]
impulse_deg <- impulse_deg[!(impulse_deg %in% downgenes)]

#Remove DEGs overlapping with non-remission patients
overlap_genes <- read.table("DESeq2_results/Timepoint_comparison/IFX/All_overlap_rem_nrem_genes.txt")
overlap_genes <- overlap_genes$V1

deg_list <- c(upgenes, downgenes, impulse_deg)
deg_list_unique <- deg_list[!(deg_list %in% overlap_genes)]
timepoint <- c("0h", "4h", "24h", "72h", "2w", "6w", "14w")

#Calculate median expression for each gene at each timepoint in remission and non-remission samples
R_expression_counts <- data.frame()
NR_expression_counts <- data.frame()

for (gene in deg_list_unique) {
  data <- plotCounts(dds_counts, gene=as.character(gene), intgroup = "Time.Point", returnData=TRUE)
  R_data <- data[grep("_R", rownames(data)),]
  NR_data <- data[grep("_NR", rownames(data)),]
  
  R_gene_median <- c()
  NR_gene_median <- c()
  for (tp in timepoint) {
    tp_data <- subset(R_data, R_data$Time.Point == tp)
    tp_median <- median(tp_data$count)
    R_gene_median <- c(R_gene_median, tp_median)
    
    tp_data <- subset(NR_data, NR_data$Time.Point == tp)
    tp_median <- median(tp_data$count)
    NR_gene_median <- c(NR_gene_median, tp_median)
    
  }
  
  R_expression_counts <- rbind(R_expression_counts, R_gene_median)
  NR_expression_counts <- rbind(NR_expression_counts, NR_gene_median)
  
}
colnames(R_expression_counts) <- timepoint
colnames(NR_expression_counts) <- timepoint
rownames(R_expression_counts) <- deg_list_unique
rownames(NR_expression_counts) <- deg_list_unique

#Cluster genes based on median expression in remission patients
data_distance <- as.dist(1-cor(t(R_expression_counts), method = "spearman"))
data_hclust <- hclust(data_distance)
R_expression_counts <- R_expression_counts[c(data_hclust$order), ]
NR_expression_counts <- NR_expression_counts[c(data_hclust$order), ]

#Combine expression counts and add gene names
heatmap_expression_counts <- cbind(R_module_expression_counts, NR_module_expression_counts)

genelist <- rownames(heatmap_expression_counts)
genenames <- unique_mart[genelist, ]$Gene_name

rownames(heatmap_expression_counts) <- genenames

#Plot heatmap

condition_colors <- c(rep("#00CC00", 7), rep("#0066CC", 7))

pdf(file="Timepoint_comparison/IFX/all_timepoint_heatmap_IFX_R_pw_and_impulse_genes.pdf", height = 9, width = 7)  
par(cex.main=1.5,mar=c(1,1,1,1))

heatmap.2(as.matrix(heatmap_expression_counts), trace="none", Colv=NA, Rowv=NA,
          dendrogram = "none", density="none", scale="row", col=rev(brewer.pal(11, "RdYlBu")), 
          ColSideColors = condition_colors, cexRow=0.7, cexCol=1.5, margins = c(5,9), lhei = c(1,5.5))

legend(.7,0.98,legend=c("R", "NR") ,
       fill=c('#00CC00', '#0066CC'),cex=1.5, bty = "n")

dev.off()
