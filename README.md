# DGE.R

Jorge L. Perez-Moreno, PhD

## Description

The script performs the following analyses and tasks:

- **Differential Expression Analysis**: Utilizes both DESeq2 and Limma-Voom methods to identify differentially expressed genes across conditions or treatments.
- **PCA**: Generates principal component analysis plots to visualize the data variance and sample relationships.
- **Correlation Analysis**: Identifies the top positively and negatively correlated genes with specific phenotypes.
- **WGCNA**: Constructs a gene co-expression network to find modules of highly correlated genes and relate them to external traits.
- **Visualization**: Produces heatmaps, volcano plots, MA plots, PCA plots, and Venn diagrams to visualize the results.

## Expected Inputs

- **Quantification Files**: Salmon quantification files (`quant.sf`) for each sample.
- **Metadata File**: A CSV file containing metadata for samples, including conditions, treatments, and any other relevant information.

## Outputs

- **Differential Expression Results**: CSV files with DESeq2 and Limma differential expression analysis results.
- **PCA Plots**: PCA plots for individual experiments and a combined plot for all experiments.
- **Heatmaps**: Heatmaps showing expression patterns of significant genes.
- **Volcano and MA Plots**: For both DESeq2 and Limma analyses.
- **Venn Diagrams**: Illustrating overlaps between DESeq2 and Limma significant genes.
- **Correlation Plots and Files**: Detailing top correlated genes with specified phenotypes.
- **WGCNA Results**: Module-trait relationships plots and related CSV files.

## Installation of Dependencies

The following R packages are required to run the analysis pipeline:

```r
install.packages(c("DESeq2", "limma", "edgeR", "ggplot2", "pheatmap", "VennDiagram", "EnhancedVolcano", "tximport", "corrplot", "WGCNA", "foreach", "doParallel"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "limma", "edgeR"))
```

Ensure all packages are successfully installed and loaded before running the script.

## Running the Script

1. Adjust the `metadata_path` and `samples_root_dir` variables to point to your metadata CSV file and the directory containing Salmon quantification files, respectively.
2. Confirm that the metadata file format matches expected columns and format.
3. Run the script in R or RStudio environment.

## Note

This DGE pipeline is designed for flexibility and can be adapted to specific research needs. Ensure you have adequate computing resources for processing, especially for large datasets.

## Support

For questions or issues, please open an issue on this GitHub repository.

