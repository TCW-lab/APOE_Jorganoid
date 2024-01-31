# APOE_Jorganoid
Repository for APOE Jorganoid project

In /scripts/simpleaf are scripts to
  1. Generate a simpleaf index with the provided reference (hg38).
  2. Use simpleaf quant to map the sequencing reads in provided FASTQ files against the provided reference and quantifies the mapping records to generate a gene count matrix.
There are individual simpleaf quant scripts for each of the 24 samples' files.

/scripts/01-get_seurat_object.R is the main R script used in this analysis. This contains code that first reads in the outputted matrix files from the simpleaf quants code, generates a 'knee' plot to assess quality of alignment, and build a Seurat object for each of the samples. Furthermore, this contains code that combines the individual Seurat objects into one object, performs log-normalization, finding variable features, log-scaling of data to prepare the data for PCA. There is code that performs PCA, plots the results with a DimPlot() and ElbowPlot() (for scree plot), then code to run UMAP analysis.

/scripts/02-Final_QC.R is the main R script used to perform the Quality Control (QC). This contains code that first reads in the individual Seurat files created in the previous script, labels them by their APOE genotype and sample identifiers, then visualizes the before QC statistics. The script then performs the QC, trimming out certain cells / genes according to specific criteria, then visualizes the results of the QC. It then splits the combined Seurat object by sample ID, and corrects for ambient RNA in each individual's data; combining them back into one object after. The script then performs dimensionality reduction and analysis. Please refer to the R Markdown file / word document to get a detailed annotation of the entire pipeline.
