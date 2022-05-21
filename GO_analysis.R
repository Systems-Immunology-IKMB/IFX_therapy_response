library(org.Hs.eg.db)
library(topGO)
library(GO.db)
library(stringr)
library(ggplot2)
library(reshape2)

#Code for GO term enrichment analysis------

#This code is for the GO enrichment analysis of co-expression modules. Same code can be used for the GO analysis of DEGs by changing the gene set (foreG).

#Read gene set and background gene set for enrichment analysis
gene_module_assignment <- read.csv("Coexpression_analysis/gene_module.txt", header = TRUE, sep = '\t')
expression_counts <- read.csv("count_files/IFX_R_NR_sf_normalized_counts.txt", header = TRUE, 
                              sep = '\t', row.names = 1)

all_genes <- rownames(expression_counts)

sig_genes <- gene_module_assignment$Gene

modules <- unique(gene_module_assignment$ModuleColor)

inUniverse <- all_genes %in% all_genes

tab <- vector('list', length(modules))
for (i in 1:length(modules)) {
  foreG <- as.character(gene_module_assignment[gene_module_assignment$ModuleColor == modules[i],]$Gene)
  inSelection <- all_genes %in% foreG
  alg <- factor( as.integer( inSelection[inUniverse] ) )
  names(alg) <- all_genes[inUniverse]
  
  
  #Create topGO object
  tgd <- new( "topGOdata", ontology="BP", allGenes = alg, nodeSize=5,
              annot=annFUN.org, mapping="org.Hs.eg.db", ID = "ensembl" )
  
  #Run Fisher test
  resultTopGO.elim <- runTest(tgd, algorithm = "elim", statistic = "Fisher" )
  resultTopGO.classic <- runTest(tgd, algorithm = "classic", statistic = "Fisher" )
  
  #Extract results
  tab[[i]] <- GenTable(tgd, Fisher.elim = resultTopGO.elim, 
                         Fisher.classic = resultTopGO.classic,
                         orderBy = "Fisher.elim" , topNodes = length(nodes(graph(tgd))))


}

#Add gene list contributing to the significantly enriched GO terms
for (m in 1:length(modules)) {
  go_results <- tab[[m]]
  foreG <- as.character(gene_module_assignment[gene_module_assignment$ModuleColor == modules[m],]$Gene)
  go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
  go_results$Genes <- NA
  for(i in 1:nrow(go_results)){
    go_id <- as.vector(go_results[i,1])
    alleges <- get(go_id, org.Hs.egGO2ALLEGS)
    genes <- unlist(mget(alleges, org.Hs.egENSEMBL))
    genes_in_cat <- intersect(as.vector(genes), as.vector(foreG))

    gene_sym_in_cat <- as.vector(unlist(mget(unlist(mget(genes_in_cat, org.Hs.egENSEMBL2EG)), org.Hs.egSYMBOL)))

    gene_sym_in_cat_str <- ""
    
    if(length(genes_in_cat) > 0){
      for(j in 1:length(gene_sym_in_cat)){
        gene_sym_in_cat_str <- paste(gene_sym_in_cat_str, gene_sym_in_cat[j], sep = ',')
      }
    }
    if (gene_sym_in_cat_str == ""){
      gene_sym_in_cat_str <- NA
    }

    go_results$Genes[i] <- gene_sym_in_cat_str

  }
  write.table(as.data.frame(go_results), file = paste("Coexpression_analysis/module_GO/GO_result_", modules[m], '_with_genenames.txt', sep = ''), sep = '\t', quote = FALSE)  
}

#Code for generating GO dot plot as in Figure 2G, Figure S2A, Figure S2B, Figure 3B, Figure 4G, Figure S7C and Figure S7D---------

#Function to filter GO results to only include top_num significantly enriched gene sets with unique gene sets
filter_go_results <- function(go_results, top_num){
  go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
  go_results <- go_results[order(go_results$Fisher.elim), ] 
  exclude_indices <- c()
  unique_gene_sets <- c()
  for(i in 1:length(go_results$GO.ID)){

    if(as.vector(go_results$Genes[i]) %in% unique_gene_sets){
      exclude_indices <- c(exclude_indices, i)
    }else{
      unique_gene_sets <- c(unique_gene_sets, as.vector(go_results$Genes[i]))
    }
  }
  if(length(exclude_indices) > 0){
    go_results_filtered <- go_results[-exclude_indices, ]
  }else{
    go_results_filtered <- go_results
  }
  if(length(go_results_filtered$GO.ID) > top_num){
    top <- go_results_filtered[1:top_num,]
  }else{
    top <- go_results_filtered
  }
  
  top$p <- -1*log10(top$Fisher.elim)
  top$Term <- str_wrap(Term(as.vector(top$GO.ID)), width = 60)
  top$Term2 <- reorder(top$Term, top$p)
  
  return(top)
}

#Function to generate plotting data for GO terms enriched in modules/timeponts/conditions
get_module_go_data <- function(modules, n_terms){
  plist <- vector('list', length(modules))
  go_list <- c()
  go_terms <- c()
  for (i in 1:length(modules)) {
    
    go_results <- read.csv(paste("GO_result_", modules[i], "_with_genenames.txt", sep = ''), 
                           header = TRUE, sep = '\t')
    go_results <- subset(go_results, go_results$Fisher.elim < 0.05)
    go_results <- subset(go_results, go_results$Significant > 1)
    
    top_go_results <- filter_go_results(go_results, n_terms)
    
    go_list <- c(go_list, as.character(top_go_results$GO.ID))
    go_terms <- c(go_terms, as.character(top_go_results$Term2))
    plist[[i]] <- go_results
    
  }
  
  p_vector <- c()
  n_vector <- c()
  for (i in 1:length(modules)) {
    p_vector <- c(p_vector, paste("p", i, sep = '_'))
    n_vector <- c(n_vector, paste("n", i, sep = '_'))
  }
  
  go_data <- data.frame(GO.ID=go_list, Term=go_terms)
  go_data[, p_vector] <- NA
  go_data[, n_vector] <- NA
  
  for (i in 1:length(modules)) {
    df <- plist[[i]]
    rownames(df) <- as.character(df$GO.ID)
    module_go_data <- df[as.character(go_data$GO.ID),]
    go_data[, (i+2)] <- module_go_data$Fisher.elim
    go_data[, (i+2+length(modules))] <- module_go_data$Significant/module_go_data$Annotated
  }
  
  return(go_data)
}

#Generate got plot for selected significantly enriched GO terms

R_NR_differentially_preserved_modules <- c("navajowhite2", "plum")
R_NR_differentially_preserved_module_go_data <- get_module_go_data(R_NR_differentially_preserved_modules, 10)

p_vector <- colnames(R_NR_differentially_preserved_module_go_data)[grep("p_", colnames(R_NR_differentially_preserved_module_go_data))]
n_vector <- colnames(R_NR_differentially_preserved_module_go_data)[grep("n_", colnames(R_NR_differentially_preserved_module_go_data))]
plot_data <- melt(R_NR_differentially_preserved_module_go_data[, c("GO.ID", "Term", p_vector)], id.vars = c("GO.ID","Term"))
plot_data$n <- unlist(R_NR_differentially_preserved_module_go_data[, n_vector])
plot_data$Module <- R_NR_differentially_preserved_modules[as.numeric(sub("p_", "", plot_data$variable))]
plot_data$p <- -1*log10(plot_data$value)
plot_data$Term2 <- reorder(plot_data$Term, 1:nrow(plot_data))
plot_data <- subset(plot_data, !is.na(plot_data$Term))

pdf("Selected_modules_GO_plot.pdf", width = 10, height = 12)
p <- ggplot(plot_data, mapping = aes(x=Module, y=Term2, size=n, color=p))
p <- p + geom_point()
p <- p + scale_x_discrete(limits=R_NR_differentially_preserved_modules, labels=c("M6", "M7", "M12", "M14", "M16"))
p <- p + xlab("Module") + ylab("GO Term")
p <- p + scale_colour_gradient(high="#990000", low="#FF9999")
p <- p + theme_bw() + theme(axis.text.y = element_text(hjust = 1, size=12, color = "black"), 
                            axis.text.x = element_text(size=14, color = "black", angle=45, hjust=1), 
                            axis.title = element_text(size = 20))
p
dev.off()



