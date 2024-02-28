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

# Register parallel backend to speed up computations
registerDoParallel(cores = 8) # Adjust based on system's capabilities

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
heatmap_dir_deseq2 <- create_dir(file.path(output_dir_deseq2, "DESeq2_Heatmaps"))
heatmap_dir_limma <- create_dir(file.path(output_dir_limma, "Limma_Heatmaps"))
overlap_dir <- create_dir(file.path(output_dir, "Overlap_Analysis"))
summary_dir_deseq2 <- create_dir(file.path(output_dir_deseq2, "Summary"))
summary_dir_limma <- create_dir(file.path(output_dir_limma, "Summary"))

# Import Data Using Tximport
txi <- tximport(file.path(samples_root_dir, metadata$filename, "quant.sf"), type = "salmon", txOut = TRUE)

# Prepare DESeq DataSet
dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ Treatment)
dds <- DESeq(dds)

# Normalizing Counts for Limma-Voom
v <- voom(counts(dds), model.matrix(~ Treatment, data = metadata), plot = FALSE)
fit <- lmFit(v, model.matrix(~ Treatment, data = metadata))


# FUNCTION Principal Component Analysis (PCA) for Each Experiment and All Experiments Together
performPCA <- function(dds, experiment, output_dir, is_combined = FALSE) {
  rld <- rlog(dds, blind=FALSE)
  data <- plotPCA(rld, intgroup=c("Experiment", "Stage"), returnData=TRUE)
  
  plot_title <- if (is_combined) "PCA for All Experiments" else paste("PCA for", experiment)
  p <- ggplot(data, aes(x=PC1, y=PC2, color=Experiment, shape=Stage)) +
    geom_point(size=3) + 
    ggtitle(plot_title)
  
  # Determine file name based on whether it's a combined PCA or individual
  file_name <- if (is_combined) "PCA_All_Experiments.png" else paste0("PCA_", experiment, ".png")
  pcaPlotPath <- file.path(output_dir, file_name)
  
  # Save the plot
  ggsave(pcaPlotPath, plot = p, width = 10, height = 8)
}

# Define the root directory for saving PCA plots
pca_output_dir <- create_dir(file.path(output_dir, "PCA_Plots"))

# Run PCA for each experiment and save the plots in their respective directories
lapply(unique(metadata$Experiment), function(exp) {
  dds_subset <- dds[, dds$Experiment == exp]
  performPCA(dds_subset, exp, pca_output_dir)
})

# Additionally, perform and save PCA for all experiments together
performPCA(dds, "All_Experiments", pca_output_dir, is_combined = TRUE)



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
			
            # Limma MA Plot Correct Indentation
            png(file.path(contrast_dir_limma, paste0(contrast_name, "_MA_Limma.png")))
            limmaMAPlot <- function(fit, coef = 1) {
                plot(fit$Amean, fit$coefficients[, coef], pch = 20, xlab = "Average Expression", ylab = "Log Fold Change", main = "MA Plot Limma")
                abline(h = c(-logfc_threshold, logfc_threshold), col = "red")
            }
            limmaMAPlot(fit2, coef = 1)
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



# FUNCTION for correlation analysis and plotting
# Assuming 'dds' and 'fit' are already created and 'metadata' is loaded

# Function to get top genes based on variance
getTopVarianceGenes <- function(data, topN = 500) {
  variances <- apply(data, 1, var)
  highVarianceGenes <- names(sort(variances, decreasing = TRUE)[1:topN])
  return(highVarianceGenes)
}

# FUNCTION for Correlation Analyses
performCorrelationAnalysis <- function(geneExpressionMatrix, significantGenes, metadata, phenotype, outputFilePrefix) {
  # Ensure geneExpressionMatrix only includes significantGenes with high variance
  significantGeneData <- geneExpressionMatrix[rownames(geneExpressionMatrix) %in% significantGenes,]
  
  # Calculate correlation
  corMatrix <- cor(significantGeneData, metadata[, phenotype, drop = FALSE], use = "pairwise.complete.obs")
  
  # Identify top correlated genes
  topCorrelations <- sort(corMatrix, decreasing = TRUE)
  topPositives <- head(topCorrelations, 20)
  topNegatives <- tail(topCorrelations, 20)
  
  # Plotting and Saving
  plotAndSaveCorrelations <- function(corSubset, suffix) {
    png(paste0(outputFilePrefix, "_", suffix, "_correlation.png"))
    corrplot(corSubset, method = "circle")
    dev.off()
    write.csv(as.data.frame(corSubset), paste0(outputFilePrefix, "_", suffix, "_correlation.csv"))
  }
  
  plotAndSaveCorrelations(topPositives, "Top20_Positive")
  plotAndSaveCorrelations(topNegatives, "Top20_Negative")
}


# FUNCTION to perform variance filtering and correlation analysis
performAnalysisForContrast <- function(contrastName, geneData, significantGenes, metadata, outputDir) {
  # Filter for high variance genes
  highVarianceGenes <- getTopVarianceGenes(geneData[significantGenes, ], 500)
  
  # Perform correlation analysis
  outputFilePrefix <- file.path(outputDir, contrastName)
  performCorrelationAnalysis(geneData, highVarianceGenes, metadata, '20E', outputFilePrefix)
  performCorrelationAnalysis(geneData, highVarianceGenes, metadata, 'Rvalue', outputFilePrefix)
}



# Loop through contrasts to perform analysis
for (experiment in experiments) {
  stages <- unique(metadata$Stage[metadata$Experiment == experiment])
  
  for (i in 1:length(stages)) {
    for (j in (i + 1):length(stages)) {
      contrastName <- paste(experiment, stages[i], "vs", stages[j], sep = "_")
      
      # DESeq2 Analysis
      resDESeq2 <- results(dds, contrast=c("Treatment", paste(experiment, stages[i], sep = "_"), paste(experiment, stages[j], sep = "_")))
      sigGenesDESeq2 <- rownames(resDESeq2)[which(resDESeq2$padj < pvalue_threshold)]
      vstDataDESeq2 <- assay(vst(dds, blind=FALSE))
      
      performAnalysisForContrast(contrastName, vstDataDESeq2, sigGenesDESeq2, metadata, output_dir_deseq2)
      
      # Limma Analysis
      contrastMatrix <- makeContrasts(Contrasts = paste(experiment, stages[i], sep = "_") - paste(experiment, stages[j], sep = "_"), levels = colnames(v))
      fit2 <- contrasts.fit(fit, contrastMatrix)
      fit2 <- eBayes(fit2)
      sigGenesLimma <- rownames(fit2)[which(fit2$table$adj.P.Val < pvalue_threshold)]
      logCPMDataLimma <- cpm(fit, log = TRUE)
      
      performAnalysisForContrast(contrastName, logCPMDataLimma, sigGenesLimma, metadata, output_dir_limma)
    }
  }
}



# FUNCTION to filter out low-expressed genes
filterLowExpressedGenes <- function(exprData) {
    cutoff <- log2(10) # Adjust based on your dataset, here using log2 CPM of 10 as an example
    filteredData <- exprData[rowMeans(log2(exprData + 1)) > cutoff, ]
    return(filteredData)
}

# FUNCTION to select a subset of genes based on variance
selectSubsetGenesForWGCNA <- function(exprData) {
    gene_variance <- apply(exprData, 1, var)
    high_variance_genes <- names(sort(gene_variance, decreasing = TRUE)[1:5000])
    return(high_variance_genes)
}

# FUNCTION for WGCNA analysis, adapted to include DESeq2 and Limma data processing
performWGCNAAndSave <- function(exprData, metadata, experimentName, outputDirType) {
    filteredExprData <- filterLowExpressedGenes(exprData)
    high_variance_genes <- selectSubsetGenesForWGCNA(filteredExprData)
    finalExprData <- filteredExprData[high_variance_genes, ]
    
    # Choose power
    sft <- pickSoftThreshold(finalExprData, powerVector = (1:10), verbose = 5)
    softPower <- sft$powerEstimate
    
    # Build network
    adj <- adjacency(finalExprData, power = softPower)
    TOM <- TOMsimilarity(adj)
    geneTree <- hclust(as.dist(1-TOM), method = "average")
    dynamicMods <- cutreeDynamic(dendro = geneTree, deepSplit = 2, minClusterSize = 30)
    mergedDynamicMods <- mergeCloseModules(finalExprData, dynamicMods, cutHeight = 0.25, verbose = 3)
    
    # Relate modules to traits
    traitData <- metadata[, c("20E", "Rvalue")]
    MEList <- moduleEigengenes(finalExprData, colors = mergedDynamicMods$colors)
    MEs <- MEList$eigengenes
    METraitsRelation <- cor(MEs, traitData, use = "pairwise.complete.obs")
    METraitsPvalue <- corPvalueStudent(METraitsRelation, nrow(finalExprData))
    
    # Save results
    write.csv(METraitsRelation, file.path(output_dir, outputDirType, paste0(experimentName, "_METraitsRelation.csv")))
    write.csv(METraitsPvalue, file.path(output_dir, outputDirType, paste0(experimentName, "_METraitsPvalue.csv")))
    
    png(file.path(output_dir, outputDirType, paste0(experimentName, "_ModuleTraitRelationship.png")))
    corrplot(METraitsRelation, method = "circle")
    dev.off()
}

# Iterate over experiments for WGCNA on DESeq2 and Limma results
foreach(experiment = unique(metadata$Experiment)) %dopar% {
    experiment_metadata <- metadata[metadata$Experiment == experiment, ]
    
    # DESeq2 Data Preparation
    dds_sub <- dds[, dds$Experiment == experiment]
    vst_data <- assay(vst(dds_sub, blind = FALSE))
    performWGCNAAndSave(vst_data, experiment_metadata, paste0(experiment, "_DESeq2"), output_dir_deseq2)
    
    # Limma Data Preparation
    logCPM_data <- cpm(fit, log = TRUE)
    performWGCNAAndSave(logCPM_data, experiment_metadata, paste0(experiment, "_Limma"), output_dir_limma)
}


# FUNCTION Enhanced Visualization for WGCNA
plotModuleTraitRelationships <- function(METraitsRelation, outputDir, experimentName) {
  png(file.path(outputDir, paste0(experimentName, "_ModuleTraitRelationships.png")))
  corrplot(METraitsRelation, method = "circle")
  dev.off()
}


# After the WGCNA analysis and saving of METraitsRelation.csv for each experiment
foreach(experiment = unique(metadata$Experiment)) %dopar% {
  experiment_metadata <- metadata[metadata$Experiment == experiment, ]

  # Define the directory where the METraitsRelation.csv file is saved for each experiment
  outputDirTypeDESeq2 <- file.path(output_dir, "DESeq2", paste0(experiment, "_DESeq2"))
  outputDirTypeLimma <- file.path(output_dir, "Limma", paste0(experiment, "_Limma"))

  # Read the METraitsRelation.csv file for DESeq2 and Limma results
  METraitsRelationDESeq2 <- read.csv(file.path(outputDirTypeDESeq2, "_METraitsRelation.csv"))
  METraitsRelationLimma <- read.csv(file.path(outputDirTypeLimma, "_METraitsRelation.csv"))

  # Plot and save the module-trait relationships using the defined function
  plotModuleTraitRelationships(METraitsRelationDESeq2, outputDirTypeDESeq2, paste0(experiment, "_DESeq2"))
  plotModuleTraitRelationships(METraitsRelationLimma, outputDirTypeLimma, paste0(experiment, "_Limma"))
}


# FUNCTION Summarize Significant Findings
summarizeSignificantFindings <- function(dds, contrastDirs) {
  significantGenesList <- list()
  for (dir in contrastDirs) {
    resultsPath <- file.path(dir, "DESeq2_results.tsv")
    if (file.exists(resultsPath)) {
      results <- read.csv(resultsPath)
      sigResults <- results[results$padj < pvalue_threshold, ]
      sigGenesUp <- sum(sigResults$log2FoldChange > logfc_threshold)
      sigGenesDown <- sum(sigResults$log2FoldChange < -logfc_threshold)
      significantGenesList[[dir]] <- list(Upregulated=sigGenesUp, Downregulated=sigGenesDown)
    }
  }
  return(significantGenesList)
}

# FUNCTION Add a section to summarize WGCNA Findings
summarizeWGCNAFindings <- function(outputDir, experimentNames) {
  summaryList <- list()
  for (experimentName in experimentNames) {
    relationPath <- file.path(outputDir, paste0(experimentName, "_METraitsRelation.csv"))
    if (file.exists(relationPath)) {
      relations <- read.csv(relationPath)
      significantModules <- which.max(abs(relations$correlation))
      summaryList[[experimentName]] <- significantModules
    }
  }
  return(summaryList)
}


# FUNCTION Summarize and save significant findings for DESeq2 and Limma
saveSignificantFindingsSummary <- function(outputDirType, method, contrastDirs) {
  # Use the 'summarizeSignificantFindings' function here
  significantGenesList <- summarizeSignificantFindings(dds, contrastDirs)
  summaryFilePath <- file.path(outputDirType, "Summary", paste0(method, "_SignificantFindingsSummary.csv"))
  # Convert the list to a data frame for easier CSV writing, if necessary
  # This step may require transformation of 'significantGenesList' to a suitable format
  write.csv(significantGenesList, summaryFilePath, row.names = TRUE)
}

# calls to save significant findings summaries
contrastDirsDESeq2 <- list.dirs(path = output_dir_deseq2, full.names = TRUE, recursive = FALSE)
contrastDirsLimma <- list.dirs(path = output_dir_limma, full.names = TRUE, recursive = FALSE)

saveSignificantFindingsSummary(output_dir_deseq2, "DESeq2", contrastDirsDESeq2)
saveSignificantFindingsSummary(output_dir_limma, "Limma", contrastDirsLimma)


# FUNCTION Save WGCNA findings summary
saveWGCNAFindingsSummary <- function(outputDir, experimentNames, method) {
  summaryList <- summarizeWGCNAFindings(outputDir, experimentNames)
  summaryFilePath <- file.path(outputDir, "Summary", paste0(method, "_WGCNAFindingsSummary.csv"))
  # Convert the list to a data frame for easier CSV writing, if necessary
  # This step may require transformation of 'summaryList' to a suitable format
  write.csv(summaryList, summaryFilePath, row.names = TRUE)
}

# Call to save WGCNA findings summaries
experimentNamesDESeq2 <- sub(pattern = output_dir_deseq2, replacement = "", x = contrastDirsDESeq2)
experimentNamesLimma <- sub(pattern = output_dir_limma, replacement = "", x = contrastDirsLimma)

saveWGCNAFindingsSummary(output_dir_deseq2, experimentNamesDESeq2, "DESeq2")
saveWGCNAFindingsSummary(output_dir_limma, experimentNamesLimma, "Limma")
