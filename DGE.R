# Required Libraries
library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(EnhancedVolcano)
library(tximport)
library(corrplot)
library(WGCNA)
library(foreach)
library(doParallel)


# Set Parameters
pvalue_threshold <- 0.05
logfc_threshold <- log2(1.5)
metadata_path <- "Sample_Metadata.csv"
samples_root_dir <- "/path/to/samples/" # Adjust this path

# Read Metadata
metadata <- read.csv(metadata_path)
metadata$Treatment <- with(metadata, paste(Experiment, Stage, sep="_"))
metadata$Experiment_Stage <- with(metadata, paste(Experiment, Stage, sep="_"))

# Helper Function to Create Directories
create_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  return(path)
}

# Define Output Directories
output_dir <- create_dir(file.path(samples_root_dir, "Analysis_Results"))
output_dir_deseq2 <- create_dir(file.path(output_dir, "DESeq2"))
output_dir_limma <- create_dir(file.path(output_dir, "Limma"))
heatmap_dir_deseq2 <- create_dir(file.path(output_dir, "DESeq2_Heatmaps"))
heatmap_dir_limma <- create_dir(file.path(output_dir, "Limma_Heatmaps"))
overlap_dir <- create_dir(file.path(output_dir, "Overlap_Analysis"))

# Import Data Using Tximport
txi <- tximport(file.path(samples_root_dir, metadata$filename, "quant.sf"), type = "salmon", txOut = TRUE)

# Prepare DESeq DataSet
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ Treatment)
dds <- DESeq(dds)

# Normalizing Counts for Limma-Voom
v <- voom(counts(dds), model.matrix(~ Treatment, data = metadata), plot = FALSE)
fit <- lmFit(v, model.matrix(~ Treatment, data = metadata))

# Experiments, Stages, and Contrast Comparisons
experiments <- unique(metadata$Experiment)

for (experiment in experiments) {
    stages <- unique(metadata$Stage[metadata$Experiment == experiment])
    
    # DESeq2 and Limma Analysis for Each Contrast
    for (i in 1:length(stages)) {
        for (j in (i + 1):length(stages)) {
            contrast_name <- paste(experiment, stages[i], "vs", stages[j], sep = "_")
            contrast_dir_deseq2 <- create_dir(file.path(output_dir_deseq2, contrast_name))
            contrast_dir_limma <- create_dir(file.path(output_dir_limma, contrast_name))
            contrast_dir_overlap <- create_dir(file.path(overlap_dir, contrast_name))
            
            # DESeq2 Analysis
            res_deseq2 <- results(dds, contrast = c("Treatment", paste(experiment, stages[i], sep = "_"), paste(experiment, stages[j], sep = "_")))
            deseq2_sig_genes <- rownames(res_deseq2)[which(res_deseq2$padj < pvalue_threshold)]
            
            # Limma Analysis
            contrast_matrix <- makeContrasts(Contrasts = paste(experiment, stages[i], sep = "_") - paste(experiment, stages[j], sep = "_"), levels = colnames(v))
            fit2 <- contrasts.fit(fit, contrast_matrix)
            fit2 <- eBayes(fit2)
            limma_sig_genes <- rownames(fit2)[which(fit2$table$adj.P.Val < pvalue_threshold)]
            
            # Save DESeq2 and Limma results
            write.csv(as.data.frame(res_deseq2), file.path(contrast_dir_deseq2, paste0(contrast_name, "_DESeq2_results.tsv")))
            write.csv(as.data.frame(fit2$table), file.path(contrast_dir_limma, paste0(contrast_name, "_Limma_results.tsv")))
            
            # EnhancedVolcano for DESeq2
            png(file.path(contrast_dir_deseq2, paste0(contrast_name, "_Volcano_DESeq2.png")))
            EnhancedVolcano(res_deseq2,
                            lab = rownames(res_deseq2),
                            x = 'log2FoldChange',
                            y = 'pvalue',
                            title = paste("Volcano Plot DESeq2", contrast_name),
                            pCutoff = pvalue_threshold,
                            FCcutoff = logfc_threshold)
            dev.off()

            # MA Plot for DESeq2
            png(file.path(contrast_dir_deseq2, paste0(contrast_name, "_MA_DESeq2.png")))
            plotMA(res_deseq2, main = paste("MA Plot DESeq2", contrast_name), ylim = c(-2, 2))
            dev.off()

            # EnhancedVolcano for Limma
            png(file.path(contrast_dir_limma, paste0(contrast_name, "_Volcano_Limma.png")))
            EnhancedVolcano(fit2,
                            lab = rownames(fit2$table),
                            x = 'logFC',
                            y = 'adj.P.Val',
                            title = paste("Volcano Plot Limma", contrast_name),
                            pCutoff = pvalue_threshold,
                            FCcutoff = logfc_threshold)
            dev.off()
			
			# Limma MA Plot
			png(file.path(contrast_dir_limma, paste0(contrast_name, "_MA_Limma.png")))
			# Assuming fit2 contains the eBayes() output from Limma
			limmaMAPlot <- function(fit, coef = 1) {
				plot(fit$Amean, fit$coefficients[,coef], pch=20, xlab="Average Expression", ylab="Log Fold Change", main="MA Plot Limma")
				abline(h=c(-logfc_threshold, logfc_threshold), col="red")
			}
			limmaMAPlot(fit2, coef = 1) # Adjust 'coef' if necessary based on your contrasts
			dev.off()

            # Heatmaps for DESeq2
            if (length(deseq2_sig_genes) > 0) {
                vst_data <- assay(vst(dds[rownames(dds) %in% deseq2_sig_genes, ], blind = FALSE))
                png(file.path(heatmap_dir_deseq2, paste0(contrast_name, "_heatmap_DESeq2.png")))
                pheatmap(vst_data, main = paste("DESeq2 Heatmap", contrast_name))
                dev.off()
            }

            # Heatmaps for Limma
            if (length(limma_sig_genes) > 0) {
                # Limma uses logCPM - need to normalize or directly use v$E
                sig_genes_data <- v$E[limma_sig_genes, , drop = FALSE]
                norm_data <- t(scale(t(sig_genes_data)))  # Normalizing gene expression
                png(file.path(heatmap_dir_limma, paste0(contrast_name, "_heatmap_Limma.png")))
                pheatmap(norm_data, main = paste("Limma Heatmap", contrast_name))
                dev.off()
            }

            # Overlap Analysis with Venn Diagram
            venn_data <- list(DESeq2 = deseq2_sig_genes, Limma = limma_sig_genes)
            png(file.path(contrast_dir_overlap, paste0(contrast_name, "_Venn.png")))
            draw.pairwise.venn(area1 = length(deseq2_sig_genes), area2 = length(limma_sig_genes), cross.area = length(intersect(deseq2_sig_genes, limma_sig_genes)), category = c("DESeq2", "Limma"), lty = "blank")
            dev.off()
        }
    }
}



# Define a function for correlation analysis and plotting
performCorrelationAnalysis <- function(geneExpressionMatrix, metadata, outputFilePrefix) {
  # Calculate correlation between gene expression and phenotype data
  cor_20E <- cor(geneExpressionMatrix, metadata$`20E`, use = "pairwise.complete.obs")
  cor_Rvalue <- cor(geneExpressionMatrix, metadata$Rvalue, use = "pairwise.complete.obs")

  # Plot correlation matrices
  png(paste0(outputFilePrefix, "_20E_correlation.png"))
  corrplot(cor_20E, method = "circle", title = "Correlation with 20E")
  dev.off()

  png(paste0(outputFilePrefix, "_Rvalue_correlation.png"))
  corrplot(cor_Rvalue, method = "circle", title = "Correlation with Rvalue")
  dev.off()

  # Save correlation results to files
  write.csv(cor_20E, paste0(outputFilePrefix, "_20E_correlation.csv"))
  write.csv(cor_Rvalue, paste0(outputFilePrefix, "_Rvalue_correlation.csv"))
}

# Perform correlation analysis after DESeq2 and Limma analyses for each contrast
for (experiment in experiments) {
    stages <- unique(metadata$Stage[metadata$Experiment == experiment])
    
    for (i in 1:length(stages)) {
        for (j in (i + 1):length(stages)) {
            contrast_name <- paste(experiment, stages[i], "vs", stages[j], sep = "_")
            outputFilePrefix_deseq2 <- file.path(output_dir_deseq2, contrast_name, contrast_name)
            outputFilePrefix_limma <- file.path(output_dir_limma, contrast_name, contrast_name)

            # Assuming vst_data and sig_genes_data are prepared as before for DESeq2 and Limma respectively
            
            # DESeq2 Correlation Analysis
            if (length(deseq2_sig_genes) > 0) {
                performCorrelationAnalysis(vst_data[deseq2_sig_genes, ], metadata, outputFilePrefix_deseq2)
            }

            # Limma Correlation Analysis
            if (length(limma_sig_genes) > 0) {
                norm_data <- t(scale(t(sig_genes_data)))  # Assuming normalization is done
                performCorrelationAnalysis(norm_data, metadata, outputFilePrefix_limma)
            }
        }
    }
}

# Register parallel backend to speed up computations
registerDoParallel(cores = 4) # Adjust based on your system's capabilities

# After normalizing counts for Limma-Voom and defining output directories
# Select a subset of genes based on variance
selectSubsetGenesForWGCNA <- function(exprData) {
  gene_variance <- apply(exprData, 1, var) # Calculate variance for each gene
  high_variance_genes <- names(sort(gene_variance, decreasing = TRUE)[1:5000]) # Select top 5000 genes
  return(high_variance_genes)
}

# Function to perform WGCNA and save module-trait relationships
performWGCNAAndSave <- function(exprData, traits, outputDir, experimentName) {
  # Ensure numeric encoding of factors
  traits$Treatment <- as.numeric(as.factor(traits$Treatment))
  
  # WGCNA: Choose power
  powers <- c(1:10)
  sft <- pickSoftThreshold(exprData, powerVector = powers, verbose = 5)
  softPower <- sft$powerEstimate
  
  # Adjacency and TOM matrix
  adj <- adjacency(exprData, power = softPower)
  TOM <- TOMsimilarity(adj)
  
  # Gene dendrogram and module detection
  geneTree <- hclust(as.dist(1-TOM), method = "average")
  minModuleSize <- 30
  dynamicMods <- cutreeDynamic(dendro = geneTree, deepSplit = 2, minClusterSize = minModuleSize)
  moduleColors <- labels2colors(dynamicMods)
  
  # Calculate eigengenes
  MEList <- moduleEigengenes(exprData, colors = moduleColors)
  MEs <- MEList$eigengenes
  
  # Relate eigengenes to traits
  METraitsRelation <- cor(MEs, traits[,c("20E", "Rvalue")], use = "p")
  
  # Plot and save module-trait relationships
  png(file.path(outputDir, paste0("ModuleTrait_", experimentName, ".png")))
  corrplot(METraitsRelation, method = "circle")
  dev.off()
  
  # Save module-trait relationships to file
  write.csv(METraitsRelation, file.path(outputDir, paste0("ModuleTrait_", experimentName, "_correlation.csv")))
}

# Iterate over experiments for WGCNA on DESeq2 and Limma results
foreach(experiment = unique(metadata$Experiment)) %dopar% {
  experiment_metadata <- metadata[metadata$Experiment == experiment, ]
  
  # DESeq2 Analysis Subset for WGCNA
  dds_sub <- dds[, dds$Experiment == experiment]
  vst_data <- assay(vst(dds_sub))
  high_variance_genes <- selectSubsetGenesForWGCNA(vst_data)
  
  performWGCNAAndSave(vst_data[high_variance_genes, ], experiment_metadata, output_dir, paste0(experiment, "_DESeq2"))
  
  # Limma Analysis Subset for WGCNA
  logCPM <- cpm(fit, log = TRUE)
  high_variance_genes_limma <- selectSubsetGenesForWGCNA(logCPM)
  
  performWGCNAAndSave(logCPM[high_variance_genes_limma, ], experiment_metadata, output_dir, paste0(experiment, "_Limma"))
}




