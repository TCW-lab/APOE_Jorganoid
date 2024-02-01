#Final QC

### Note: This is taken from Deepti Murthy and adapted by Andrew Gjelsteen,
### most of the comments are hers ###
### Andrew's comments are denoted by ##*  *##
#TCW Single-Cell RNA Seq Data
#Dataset Source: TCW_14311_run1

#Date: 11/10/22 - 11/13/22
#Date: 1/1/23 (Altering Cluster Resolution)

#Note, prior to this analysis:
#The seurat_object.rds" object was created using RNA and HTO information, including
#only those cells with at least 10 nUMI for any given HTO.
#HTODemux was run using default parameters.
#Across-sample Doublets and dropouts have been removed from the dataset.

#In this analysis, I remove cells with high mitochondrial percentage and 
## remove outliers, scale, normalize, and cluster.
###
#R version 4.3.1
sessionInfo()

#load libraries 
library(dplyr)
library(patchwork)
library(ggplot2)
library(gridExtra)
library(Matrix)
library(scales)
library(cowplot)
library(Seurat) #Seurat v.3.0
library(RColorBrewer)
library(BiocManager) #v 3.17
library(tibble)
library(patchwork)
library(flexmix)
library(miQC)
library(singleCellTK)
library(fgsea)
library(singleCellTK)
library(DropletUtils)
library(SeuratData)
library(tidyr)
.libPaths()
source("~/.Rprofile")

###* denotes changes made by AKG *###
###* Step 1: Read in the individual samples' Seurat object files            *###
###* Then, combine them into one Seurat object. This object will be used to *###
###* be a point of comparison for before/after quality control steps        *###
setwd("/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/")

dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs"
out<-'outputs/simpleafSeurat'

S1 <- readRDS(file = "./outputs/simpleafSeurat/sample1.rds")
S2 <- readRDS(file = "./outputs/simpleafSeurat/sample2.rds")
S3 <- readRDS(file = "./outputs/simpleafSeurat/sample3.rds")
S4 <- readRDS(file = "./outputs/simpleafSeurat/sample4.rds")
###* No Sample 5 (deemed poor quality by researchers)                       *###
S6 <- readRDS(file = "./outputs/simpleafSeurat/sample6.rds")
###* No Sample 7 (deemed poor quality by researchers)                       *###
S8 <- readRDS(file = "./outputs/simpleafSeurat/sample8.rds")
S9 <- readRDS(file = "./outputs/simpleafSeurat/sample9.rds")
S10 <- readRDS(file = "./outputs/simpleafSeurat/sample10.rds")
S11 <- readRDS(file = "./outputs/simpleafSeurat/sample11.rds")
S12 <- readRDS(file = "./outputs/simpleafSeurat/sample12.rds")
S13 <- readRDS(file = "./outputs/simpleafSeurat/sample13.rds")
S14 <- readRDS(file = "./outputs/simpleafSeurat/sample14.rds")
S15 <- readRDS(file = "./outputs/simpleafSeurat/sample15.rds")
S16 <- readRDS(file = "./outputs/simpleafSeurat/sample16.rds")
S17 <- readRDS(file = "./outputs/simpleafSeurat/sample17.rds")
S18 <- readRDS(file = "./outputs/simpleafSeurat/sample18.rds")
S19 <- readRDS(file = "./outputs/simpleafSeurat/sample19.rds")
S20 <- readRDS(file = "./outputs/simpleafSeurat/sample20.rds")
S21 <- readRDS(file = "./outputs/simpleafSeurat/sample21.rds")
S22 <- readRDS(file = "./outputs/simpleafSeurat/sample22.rds")
S23 <- readRDS(file = "./outputs/simpleafSeurat/sample23.rds")
S24 <- readRDS(file = "./outputs/simpleafSeurat/sample24.rds")
S25 <- readRDS(file = "./outputs/simpleafSeurat/sample25.rds")
S26 <- readRDS(file = "./outputs/simpleafSeurat/sample26.rds")

###* Proceed by merging all the samples into a single object, called organoid0
APOE22_samples_group <- c("S3", "S9", "S12", "S15", "S21", "S25")
APOE33_samples_group <- c("S1", "S8", "S14", "S16", "S19", "S23")
APOE33Ch_samples_group <- c("S2", "S6", "S13", "S18", "S20", "S24")
APOE44_samples_group <- c("S4", "S10", "S11", "S17", "S22", "S26")

sample_names_group22 <- c("sample3", "sample9", "sample12", "sample15", "sample21", "sample25")
sample_names_group33 <- c("sample1", "sample8", "sample14", "sample16", "sample19", "sample23")
sample_names_group33Ch <- c("sample2", "sample6", "sample13", "sample18", "sample20", "sample24")
sample_names_group44 <- c("sample4", "sample10", "sample11", "sample17", "sample22", "sample26")

# Loop through each sample and add a metadata column for sample labels
# for APOE22
for (i in seq_along(APOE22_samples_group)) {
  object_name <- APOE22_samples_group[i]
  sample_name <- sample_names_group22[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE22'")))
}

for (i in seq_along(APOE33_samples_group)) {
  object_name <- APOE33_samples_group[i]
  sample_name <- sample_names_group33[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33'")))
}

for (i in seq_along(APOE33Ch_samples_group)) {
  object_name <- APOE33Ch_samples_group[i]
  sample_name <- sample_names_group33Ch[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33Ch'")))
}

for (i in seq_along(APOE44_samples_group)) {
  object_name <- APOE44_samples_group[i]
  sample_name <- sample_names_group44[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE44'")))
}

################################################################################
#########@######################################################################


##* Merge all Seurat objects into one *##
##* Can easily process it as one object *##
##* 
organoid0 <- merge(S1, y = c(S2, S3, S4, S6, S8, S9, S10,
                             S11, S12, S13, S14, S15, S16, S17, S18, S19,
                             S20, S21, S22, S23, S24, S25, S26), 
                   add.cell.ids = c("S1", "S2", "S3", "S4", "S6", "S8", "S9", "S10",
                                    "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19",
                                    "S20", "S21", "S22", "S23", "S24", "S25", "S26"),
                   project = "APOE_Jorganoid"
)

organoid0 <- readRDS('outputs/simpleafSeurat/organoid0.rds')

### Note: A lot of the mitochondrial genes, when undergoing the mapping of 
### ENSEMBL_IDs to gene symbols, somehow got excluded from this process.
### In order to deal with this, we will define a list of all the ENSEMBL IDs AND
### their gene symbols for
### all 37 Human Mitochondrial genes
mt_gene_names <- c(
  "ENSG00000198888", # ND1
  "ENSG00000198763", # ND2
  "ENSG00000198804", # COX1
  "ENSG00000198712", # COX2
  "ENSG00000228253", # ATP8
  "ENSG00000198899", # ATP6
  "ENSG00000198938", # COX3
  "ENSG00000198840", # ND3
  "ENSG00000212907", # ND4L
  "ENSG00000198886", # ND4
  "ENSG00000198786", # ND5
  "ENSG00000198695", # ND6
  "ENSG00000198727", # CYTB
  "ENSG00000210049", # MT-TF
  "ENSG00000210082", # MT-TV
  "ENSG00000209082", # MT-TL1
  "ENSG00000198888", # MT-TR
  "ENSG00000210100", # MT-TN
  "ENSG00000210107", # MT-TG
  "ENSG00000210112", # MT-TL2
  "ENSG00000210117", # MT-TS1
  "ENSG00000210127", # MT-TV
  "ENSG00000210135", # MT-TE
  "ENSG00000210144", # MT-TS2
  "ENSG00000210151", # MT-TH
  "ENSG00000210154", # MT-TD
  "ENSG00000210156", # MT-TK
  "ENSG00000210164", # MT-TM
  "ENSG00000210174", # MT-TI
  "ENSG00000210184", # MT-TT
  "ENSG00000210191", # MT-TW
  "ENSG00000210194", # MT-TC
  "ENSG00000210211", # MT-TY
  "ENSG00000210228", # MT-TA
  "ENSG00000210243", # MT-TQ
  "ENSG00000198899", # MT-RNR1 (12S rRNA)
  "ENSG00000198763",  # MT-RNR2 (16S rRNA)
  "ND1", "ND2", "COX1", "COX2", "ATP8", "ATP6", "COX3", 
  "ND3", "ND4L", "ND4", "ND5", "ND6", "CYTB",
  "TF", "TV", "TL1", "TR", "TN", "TG", "TL2", "TS1", "TV",
  "TE", "TS2", "TH", "TD", "TK", "TM", "TI", "TT", "TW", 
  "TC", "TY", "TA", "TQ", "RNR1", "RNR2"
)

### let's do ambient
### RNA correction first:

# Initialize a data frame to store cell counts before and after QC
cell_counts <- data.frame(sample = character(), 
                          cells_before_QC = integer(), 
                          cells_after_QC = integer())

# Get unique sample identifiers
samples <- unique(organoid0@meta.data$sample)

# Define mitochondrial genes (adjust according to your organism)
mt_gene_names <- grep("^MT-", rownames(organoid0), value = TRUE)

# Loop over each sample
for (sample_id in samples) {
  # Subset Seurat object for the current sample
  T_QC <- subset(organoid0, subset = sample == sample_id)
  
  # Record the number of cells before QC
  cells_before <- ncol(T_QC)
  
  # Filter genes expressed in less than 3 cells
  raw_counts_mat <- T_QC@assays$RNA@counts
  genes_to_keep <- Matrix::colSums(raw_counts_mat != 0) >= 3
  T_QC <- subset(T_QC, features = rownames(T_QC)[genes_to_keep])
  
  # Calculate upper limit, set at 99.5 percentile
  upper_limit <- quantile(T_QC@meta.data$nCount_RNA, probs = 0.995)
  
  T_QC <- CalculateBarcodeInflections(
    T_QC,
    barcode.column = "nCount_RNA",
    group.column = "orig.ident",
    threshold.low = 1000,
    threshold.high = NULL
  )
  T_QC <- SubsetByBarcodeInflections(object = T_QC)
  
  # Record the number of cells after QC
  cells_after <- ncol(T_QC)
  
  # Append the cell count information to the data frame
  cell_counts <- rbind(cell_counts, data.frame(sample = sample_id,
                                               cells_before_QC = cells_before,
                                               cells_after_QC = cells_after))
}

# View the cell counts data frame
print(cell_counts)

#### SoupX Correction



######~~~~~ Deepti's SoupX Code ~~~~~###########################################
#SoupX 
library(tools)
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(fishpond)
library(scater)
library(org.Hs.eg.db)
library(Matrix)
library(biomaRt)
library(SoupX)

## initialize a dataframe which will contain the before/after SoupX correction
## cell counts
SoupXSummary <- data.frame(
  sampleID = character(24),
  before_SoupX = numeric(24),
  after_SoupX = numeric(24),
  percent_removed = numeric(24),
  stringsAsFactors = FALSE
)


# Assuming 'organoid0' is your large combined Seurat object
sample_ids <- unique(organoid0@meta.data$sample)
# sample_ids
seurat_objects <- SplitObject(organoid0, split.by = "sample")


source <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/'

### Function to convert ENSEMBL_ID names to gene symbols
getGeneSymbols <- function(ensembl_ids) {
  require(org.Hs.eg.db)
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  return(gene_symbols)
}



# Load raw counts
custom_format <- list("counts" = c("U","S","A"))

# Function to apply SoupX correction
apply_soupX_correction <- function(seurat_obj, raw_counts_directory) {
  
  # # set a "raw counts" count for each of the individual Seurat samples' objects
  # seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
  # 
  sce <- fishpond::loadFry(raw_counts_directory,
                           outputFormat = custom_format)
  # Assuming 'sce' is your Single Cell Experiment object
  # Replace Ensembl IDs in row names with gene symbols
  ensembl_ids_sce <- rownames(sce)
  gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)
  
  # Replace Ensembl IDs with gene symbols in row names
  # If no gene symbol is found, retain the original Ensembl ID
  updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
  rownames(sce) <- updated_rownames
  
  #create the seurat object by filtering for selected cells
  S1Counts <- counts(sce)
  
  raw <- S1Counts
  #SoupX Correction
  #modify genes in raw
  genes <- intersect(rownames(seurat_obj), rownames(raw))
  length(genes)
  raw <- raw[genes,]

  #run SoupX algo (this package is what I ended up using instead of DecontX)
  sc <- SoupChannel(raw, seurat_obj@assays$RNA@counts)
  sc <- setClusters(sc, seurat_obj$initial_seurat_clustering_0.5)
  sc <- autoEstCont(sc, doPlot = FALSE)
  soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix
  
  seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
  seurat_obj@assays$RNA@counts <- soup_out #USED
  
  #calculate percent of mRNA removed
  percent_removed_soup <- (1 - (sum(seurat_obj@assays$RNA@counts) / sum(seurat_obj@assays$original_counts@counts)))*100
  
  
  return(seurat_obj)
}

# Directories for raw counts of each sample
raw_counts_directories <- c("1_S16/af_quant/", "2_S8/af_quant/", "3_S9/af_quant/", "4_S14/af_quant/", "6_S23/af_quant/",
                            "8_S17/af_quant/", "9_S21/af_quant/", "10_S7/af_quant/", "11_S1/af_quant/", "12_S20/af_quant/", "13_S11/af_quant/",
                            "14_S13/af_quant/", "15_S5/af_quant/", "16_S24/af_quant/", "17_S19/af_quant/", "18_S22/af_quant/",
                            "19_S3/af_quant/", "20_S15/af_quant/", "21_S6/af_quant/", "22_S2/af_quant/", "23_S12/af_quant/", "24_S10/af_quant/",
                            "25_S18/af_quant/", "26_S4/af_quant/") # paths to quant files

organoid1 <- readRDS('outputs/simpleafSeurat/organoid1.rds')

# Process each Seurat object
for (i in seq_along(seurat_objects)) {
  sampleID <- names(seurat_objects)[i]
  print(sampleID)
  raw_counts_dir <- paste(source, raw_counts_directories[i], sep = '') 
  print(raw_counts_dir)
  original_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  
  seurat_objects[[sampleID]] <- apply_soupX_correction(seurat_objects[[sampleID]], raw_counts_dir)
  
  after_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  percent_removed <- (1 - (after_cell_count / original_cell_count)) * 100
  
  samples_summary[i, ] <- c(sampleID, original_cell_count, after_cell_count, percent_removed)
}


raw_counts_dir <- paste(source, "1_S16/af_quant/", sep = '') 
# # set a "raw counts" count for each of the individual Seurat samples' objects
# seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
# 
sce <- fishpond::loadFry(raw_counts_dir,
                         outputFormat = custom_format)
# Assuming 'sce' is your Single Cell Experiment object
# Replace Ensembl IDs in row names with gene symbols
ensembl_ids_sce <- rownames(sce)
gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)

# Replace Ensembl IDs with gene symbols in row names
# If no gene symbol is found, retain the original Ensembl ID
updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
rownames(sce) <- updated_rownames

#create the seurat object by filtering for selected cells
S1Counts <- counts(sce)


seurat_obj <- seurat_objects[["sample1"]]
colnames(seurat_obj)
colnames(S1Counts)

colnames(S1Counts) <- paste0("S1_", colnames(S1Counts))

raw <- S1Counts
#SoupX Correction
#modify genes in raw
genes <- intersect(rownames(seurat_obj), rownames(raw))
length(genes)
raw <- raw[genes,]

#normalize data
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures



#scale all
seurat_obj <- ScaleData(seurat_obj, assay = "RNA", verbose = FALSE)


seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:50, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)



#run SoupX algo (this package is what I ended up using instead of DecontX)
sc <- SoupChannel(raw, seurat_obj@assays$RNA@counts[genes,])
seurat_obj@meta.data$

sc <- setClusters(sc, seurat_obj$initial_seurat_clustering_0.5) #### Ask Deepti (mention that this clustering is called for before the cluster is actually performed)
sc <- autoEstCont(sc, doPlot = FALSE)
soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix

seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
seurat_obj@assays$RNA@counts <- soup_out #USED

#calculate percent of mRNA removed
percent_removed_soup <- (1 - (sum(seurat_obj@assays$RNA@counts) / sum(seurat_obj@assays$original_counts@counts)))*100


return(seurat_obj)

##### Deepti's code to set original idents on objects #####
##### to keep track of how many mt / ribo genes were removed #####

gene_names <- rownames(organoid0@assays$RNA@counts)
gene_names
###* Filter gene names based on mitochondrial ENSEMBL IDs/gene symbols list *###
mt_genes <- gene_names[gene_names %in% mt_gene_names]

# Print the mitochondrial genes
print(mt_genes)
###* Assign percent.mt and percent.ribo to the organoid0 object:            *###
###* 
# Set Idents, calculate MT and Ribosomal percentages
Idents(organoid0) <- organoid0$orig.ident

# Subset mt_gene_names to include only genes present in the Seurat object
valid_mt_gene_names <- mt_gene_names[mt_gene_names %in% rownames(organoid0[["RNA"]]@counts)]

# Recalculate the percentage of mitochondrial genes
organoid0[["percent.mt"]] <- PercentageFeatureSet(organoid0, features = valid_mt_gene_names)
organoid0[["percent.ribo"]] <- PercentageFeatureSet(organoid0, pattern = "^RP[SL]")

###* Let's now plot the violin plots for percent.mt, percent.ribo, etc.     *###


# organoid0 is the Seurat object and 'genotype' is a column in meta.data
unique_genotypes <- unique(organoid0@meta.data$genotype)

# Output directory
out <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/images/"  # Replace with your actual directory path

# Loop through each genotype
library(Seurat)
library(ggplot2)
library(reshape2)

rm(S1, S2, S3, S4, S6, S8, S9, S10, 
   S11, S12, S13, S14, S15, S16, S17, S18, S19, S20, 
   S21, S22, S23, S24, S25, S26) #clear recent memory

# Define the features to plot
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# Loop through each genotype and make violin plots for each APOE genotype
for (genotype in unique_genotypes) {
  # Subset Seurat object by genotype
  organoid_subset <- subset(organoid0, subset = genotype == genotype)
  
  # Create violin plots
  vps <- VlnPlot(object = organoid_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
  
  # Modify the plots to spread them out evenly
  vps <- vps + plot_layout(guides = 'collect') & theme(legend.position = 'none')
  
  # Add the genotype as the overall figure title
  combined_plot <- vps + plot_annotation(title = genotype, theme = theme(plot.title = element_text(hjust = 0.5)))
  
  # Construct the full path for the output file
  output_file_path <- file.path(out, paste0("before_QC_violin_plots_", genotype, ".png"))
  
  # Save the combined plot to the specified directory
  ggsave(output_file_path, plot = combined_plot, width = 12, height = 6, units = "in")
}

out <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/miQC"
filename = 'organoid0.rds'
###* Saving our orgnaoid0 object. Good idea to do so, so we can restart R   *###
###* prior to next step and load in the object                              *###
saveRDS(organoid0, file.path(out, filename))

################################################################################
#### Following Deepti's Pipeline here: #########################################
#### Step 2: Remove Low Quality Cells  #########################################
################################################################################

organoid0 <- readRDS('outputs/simpleafSeurat/organoid0.rds')

organoid0@active.assay
#set assay to RNA
organoid0@active.assay <- "RNA" #8,149 cells & 33,546 genes

### Before deciding cutoffs, use vln plots to visualize the number of detected
### genes expressed, number of UMIs detected, and pt.mitochondrial genes:

vps<-VlnPlot(object = organoid0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="red") +labs(x = "organoid0")
vp2<-vps[[2]]#+geom_hline(yintercept = 75000,colour="red") +geom_hline(yintercept = 500,colour="red") +labs(x = "organoid0")
vp3<-vps[[3]]#+geom_hline(yintercept = 14,colour="red") +labs(x = "organoid0")

vps <- VlnPlot(object = organoid0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
vp1 <- vps[[1]] + labs(x = "organoid0")
vp2 <- vps[[2]] + labs(x = "organoid0")
vp3 <- vps[[3]] + labs(x = "organoid0")




pdf(paste(dir, "/QC_metric_VP_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank()) 
plot(p)
dev.off()

### We can also visualize a scatter plot of nFeature_RNA vs. nCount_RNA
### which ought to have a linear relationship:

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(organoid0, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(organoid0, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(paste(dir, "/images/QC_metrics_FP_before_cutoffs.pdf", sep = ""), width = 10, height = 5)
plot1 + plot2 
dev.off()

# Expect this to not be linear prior to following QC steps

#### Remove low quality cells based on miQC threshold ##########################

#BiocManager::install("flexmix")
organoid1 <- readRDS('outputs/simpleafSeurat/organoid0.rds')
library(miQC)
library(Seurat)
library(flexmix)
library(SeuratWrappers)
setwd("/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/")

dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs"
out <- 'outputs/miQC'

#Run miQC
#note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
#Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.

organoid1 <- RunMiQC(organoid1, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.97, model.slot = "flexmix_model") 

pdf(paste(dir, "/miQC/stringent_miQC_probability_plot.pdf", sep = ""), width = 15, height = 10)
PlotMiQC(organoid1, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
dev.off()

pdf(paste(dir, "/miQC/miQC_keep_probability_plot.pdf", sep = ""), width = 15, height = 10)
PlotMiQC(organoid1, color.by = "miQC.keep")
dev.off()



#Run MiQC Filtering
organoid1 <- subset(organoid1, miQC.keep == "keep") #118,688 cells and 36,601 genes

nrow(organoid1@assays$RNA@counts) #genes in organoid1
ncol(organoid1@assays$RNA@counts) #cells in organoid1
nrow(organoid0@assays$RNA@counts) #genes in organoid0
ncol(organoid0@assays$RNA@counts) #cells in organoid0


############ Remove Additional Outliers ########################################

#Remove Additional Outliers
#Only include cells with at least 200 genes expressed
organoid2 <- subset(organoid1, subset = nFeature_RNA > 200) #number of genes detected in each cell
#7,669 cells and 24,255 genes

#Only include cells with at least 500 molecules expressed
organoid2 <- subset(organoid2, subset = nCount_RNA > 500) #number of genes detected in each cell
#7,669 cells and 24,255 genes

#Remove the upper tail observed on the "nCount_RNA" violin plot (99.5% percentile)
ub <- quantile(organoid2[["nCount_RNA"]]$nCount_RNA, probs = 0.995) 
organoid2 <- organoid2[, organoid2[["nCount_RNA"]] < ub] #7,592 cells and 24255 genes

nrow(organoid2@assays$RNA@counts) #genes in organoid2 (36,601)
ncol(organoid2@assays$RNA@counts) #cells in organoid2 (118,094)


######## Visualize cell expression profile after QC ############################

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(paste(dir, "/miQC/QC_metrics_FP_after_cutoffs.pdf", sep = ""), width = 10, height = 5)
plot1 + plot2 
dev.off()


###### Now, let's examine how many cells were recovered: #######################
# Get the list of unique genotypes
genotypes <- unique(organoid2@meta.data$genotype)

# Initialize a data frame to store the counts
cell_counts_table <- data.frame(
  Genotype = character(),
  Cells_Before_QC = integer(),
  Cells_After_QC = integer(),
  Percent_Recovered = double(),
  stringsAsFactors = FALSE
)

# Loop over each genotype to get gene counts before and after QC
for (geno in genotypes) {
  # Subset for the current genotype
  subset_organoid0 <- subset(organoid0, subset = genotype == geno)
  subset_organoid2 <- subset(organoid2, subset = genotype == geno)
  
  # Count the number of cells before QC
  cells_before_QC <- length(subset_organoid0@meta.data$nFeature_RNA)
  
  # Count the number of genes after QC
  cells_after_QC <- length(subset_organoid2@meta.data$nFeature_RNA)
  
  # Calculate percent recovered
  percent_recovered <- ifelse(cells_before_QC > 0, (cells_after_QC / cells_before_QC) * 100, 0)
  
  # Add the counts to the data frame
  cell_counts_table <- rbind(
    cell_counts_table,
    data.frame(
      Genotype = geno,
      Cells_Before_QC = cells_before_QC,
      Cells_After_QC = cells_after_QC,
      Percent_Recovered = percent_recovered
    )
  )
}

# Print the table
print(cell_counts_table)

######### save organoid2 object before redeclaring it as organoid3.rds #########
saveRDS(organoid2,file.path(out,'organoid3.rds'))

################################################################################
### Remove Ambient RNA #########################################################
################################################################################

organoid3 <- readRDS('outputs/miQC/organoid3.rds')

### Run SoupX Here ####
#SoupX 
library(tools)
library(Seurat)
library(data.table)
library(ggplot2)
library(SingleCellExperiment)
library(fishpond)
library(scater)
library(org.Hs.eg.db)
library(Matrix)
library(biomaRt)
library(SoupX)

#### Scale, Normalize, Cluster data:
#normalize data


seurat_obj <- NormalizeData(organoid3)
seurat_obj <- FindVariableFeatures(seurat_obj)



#scale all
seurat_obj <- ScaleData(seurat_obj, assay = "RNA", verbose = FALSE)


seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:50, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)







## initialize a dataframe which will contain the before/after SoupX correction
## cell counts
SoupXSummary <- data.frame(
  sampleID = character(24),
  before_SoupX = numeric(24),
  after_SoupX = numeric(24),
  percent_removed = numeric(24),
  stringsAsFactors = FALSE
)

organoid3 <- seurat_obj

# Assuming 'organoid0' is your large combined Seurat object
sample_ids <- unique(organoid0@meta.data$sample)
# sample_ids
seurat_objects <- SplitObject(organoid3, split.by = "sample")


source <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/'

### Function to convert ENSEMBL_ID names to gene symbols
getGeneSymbols <- function(ensembl_ids) {
  require(org.Hs.eg.db)
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = ensembl_ids,
                         column = "SYMBOL",
                         keytype = "ENSEMBL",
                         multiVals = "first")
  return(gene_symbols)
}



# Load raw counts
custom_format <- list("counts" = c("U","S","A"))

# Function to apply SoupX correction
apply_soupX_correction <- function(seurat_obj, raw_counts_directory) {
  
  # # set a "raw counts" count for each of the individual Seurat samples' objects
  # seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
  # 
  sce <- fishpond::loadFry(raw_counts_directory,
                           outputFormat = custom_format)
  # Assuming 'sce' is your Single Cell Experiment object
  # Replace Ensembl IDs in row names with gene symbols
  ensembl_ids_sce <- rownames(sce)
  gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)
  
  # Replace Ensembl IDs with gene symbols in row names
  # If no gene symbol is found, retain the original Ensembl ID
  updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
  rownames(sce) <- updated_rownames
  
  #create the seurat object by filtering for selected cells
  S1Counts <- counts(sce)
  
  raw <- S1Counts
  #SoupX Correction
  #modify genes in raw
  genes <- intersect(rownames(seurat_obj), rownames(raw))
  length(genes)
  raw <- raw[genes,]
  
  #run SoupX algo (this package is what I ended up using instead of DecontX)
  sc <- SoupChannel(raw, seurat_obj@assays$RNA@counts[genes,])
  sc <- setClusters(sc, seurat_obj$seurat_clusters)
  sc <- autoEstCont(sc, doPlot = FALSE)
  soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix
  
  seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
  seurat_obj@assays$RNA@counts <- soup_out #USED
  
  #calculate percent of mRNA removed
  percent_removed_soup <- (1 - (sum(seurat_obj@assays$RNA@counts) / sum(seurat_obj@assays$original_counts@counts)))*100
  
  
  return(seurat_obj)
}

# Directories for raw counts of each sample
raw_counts_directories <- c("1_S16/af_quant/", "2_S8/af_quant/", "3_S9/af_quant/", "4_S14/af_quant/", "6_S23/af_quant/",
                            "8_S17/af_quant/", "9_S21/af_quant/", "10_S7/af_quant/", "11_S1/af_quant/", "12_S20/af_quant/", "13_S11/af_quant/",
                            "14_S13/af_quant/", "15_S5/af_quant/", "16_S24/af_quant/", "17_S19/af_quant/", "18_S22/af_quant/",
                            "19_S3/af_quant/", "20_S15/af_quant/", "21_S6/af_quant/", "22_S2/af_quant/", "23_S12/af_quant/", "24_S10/af_quant/",
                            "25_S18/af_quant/", "26_S4/af_quant/") # paths to quant files


# Process each Seurat object
for (i in seq_along(seurat_objects)) {
  sampleID <- names(seurat_objects)[i]
  print(sampleID)
  raw_counts_dir <- paste(source, raw_counts_directories[i], sep = '') 
  print(raw_counts_dir)
  original_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  
  seurat_objects[[sampleID]] <- apply_soupX_correction(seurat_objects[[sampleID]], raw_counts_dir)
  
  after_cell_count <- ncol(seurat_objects[[sampleID]]@assays$RNA@counts)
  percent_removed <- (1 - (after_cell_count / original_cell_count)) * 100
  
  samples_summary[i, ] <- c(sampleID, original_cell_count, after_cell_count, percent_removed)
}

organoid5 <- merge(seurat_objs[[1]], y = c(seurat_objs[2:len(seurat_objs)]),
                   project = "APOE_Jorganoid"
)




#### This is for running SoupX on an individual Sample's Data: (as an example)


raw_counts_dir <- paste(source, "1_S16/af_quant/", sep = '') 
# # set a "raw counts" count for each of the individual Seurat samples' objects
# seurat_obj[["original_counts"]] <- CreateAssayObject(counts = seurat_obj@assays$RNA@counts) #raw counts
# 
sce <- fishpond::loadFry(raw_counts_dir,
                         outputFormat = custom_format)
# Assuming 'sce' is your Single Cell Experiment object
# Replace Ensembl IDs in row names with gene symbols
ensembl_ids_sce <- rownames(sce)
gene_symbols_sce <- getGeneSymbols(ensembl_ids_sce)

# Replace Ensembl IDs with gene symbols in row names
# If no gene symbol is found, retain the original Ensembl ID
updated_rownames <- ifelse(is.na(gene_symbols_sce[ensembl_ids_sce]), ensembl_ids_sce, gene_symbols_sce[ensembl_ids_sce])
rownames(sce) <- updated_rownames

#create the seurat object by filtering for selected cells
S1Counts <- counts(sce)


seurat_obj <- seurat_objects[["sample1"]]
colnames(seurat_obj)
colnames(S1Counts)

colnames(S1Counts) <- paste0("S1_", colnames(S1Counts))

raw <- S1Counts
#SoupX Correction
#modify genes in raw
genes <- intersect(rownames(seurat_obj), rownames(raw))
length(genes)
raw <- raw[genes,]



# 
# #run SoupX algo (this package is what I ended up using instead of DecontX)
# sc <- SoupChannel(raw, seurat_obj@assays$RNA@counts[genes,])
# seurat_obj@meta.data$
#   
#   sc <- setClusters(sc, seurat_obj$initial_seurat_clustering_0.5) #### Ask Deepti (mention that this clustering is called for before the cluster is actually performed)
# sc <- autoEstCont(sc, doPlot = FALSE)
# soup_out <- adjustCounts(sc, roundToInt = TRUE) #adjusted counts matrix
# 
# seurat_obj[["SoupX_counts"]] <- CreateAssayObject(counts = soup_out)
# seurat_obj@assays$RNA@counts <- soup_out #USED
# 
# #calculate percent of mRNA removed
# percent_removed_soup <- (1 - (sum(seurat_obj@assays$RNA@counts) / sum(seurat_obj@assays$original_counts@counts)))*100




########## Everything below this point is from Deepti's original code #########
######### As I progress through this, I will be deleting and editing parts of ##
### this code.



#preserve original counts prior to ambient RNA removal
organoid3[["original_counts"]] <- CreateAssayObject(counts = organoid3@assays$RNA@counts) #raw counts



### Removal of Homotypic Doublets ##############################################
### First, Scale and Normalize the Data ########################################
TCW <- CellCycleScoring(organoid3, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes,
                        set.ident = TRUE,
                        search=TRUE)

TCW[["CC.Difference"]] <- TCW$G2M.Score - TCW$S.Score

all.genes <- rownames(TCW)
TCW <- NormalizeData(TCW, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")

#find variable features
TCW <- FindVariableFeatures(TCW, selection.method = "vst", nfeatures = 2000)

#scale only mt
tcw <- ScaleData(TCW, assay = "RNA", vars.to.regress = "percent.mt") #visualization scaling
tcw  <- RunPCA(tcw, features = VariableFeatures(object = tcw))

pdf(paste(dir, "/pca_dim_1&2_after_first_scaling.pdf", sep = ""), width = 10, height = 10)
VizDimLoadings(tcw, dims = 1:2, reduction = "pca")
dev.off()

pdf(paste(dir, "/dim_heatmap_after_first_scaling.pdf", sep = ""), width = 10, height = 10)
DimHeatmap(tcw, dims = 1, cells = 500, balanced = TRUE)
dev.off()

pdf(paste(dir, "/ridge_plot_cell_cycle_before_scaling.pdf", sep = ""), width = 10, height = 10)
RidgePlot(tcw, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
dev.off()

pdf(paste(dir, "/pca_plot_cell_cycle_before_scaling.pdf", sep = ""), width = 10, height = 10)
tcw <- RunPCA(tcw, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))
DimPlot(tcw)
dev.off()

saveRDS(TCW,file.path(out,'organoid4.rds'))
TCW <- readRDS('outputs/miQC/organoid4.rds')

# HERE, this gave me trouble, so I'm separating out the functions into three separate regressions
#scale all
# TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)
all.genes <- rownames(TCW)
# Scale data regressing out 'percent.mt'
TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = "percent.mt", verbose = TRUE)

# Scale data regressing out 'nCount_RNA'
TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = "nCount_RNA", verbose = TRUE)

# Scale data regressing out 'CC.Difference'
TCW <- ScaleData(TCW, features = all.genes, assay = "RNA", vars.to.regress = "CC.Difference", verbose = TRUE)

####### Run it with these 3 functions having been run first thing after lab meeting #####

TCW <- ScaleData(TCW)

rm(organoid_subset, plot1, plot2, sce)



#run sctransform
#BiocManager::install("glmGamPoi") ### H E R E is where we stopped for now.
library(glmGamPoi)

TCW <- SCTransform(TCW, method = "glmGamPoi", vars.to.regress = c("percent.mt", "nCount_RNA", "CC.Difference"), verbose = FALSE)

TCW <- RunPCA(TCW, verbose = FALSE)
TCW <- RunUMAP(TCW, dims = 1:50, verbose = FALSE)
TCW <- FindNeighbors(TCW, dims = 1:50, verbose = FALSE)
TCW <- FindClusters(TCW, resolution = 0.4, verbose = FALSE)


# split object (by samples)


rm(reffy, reffo, ref_assay, ref_metadata)

###### Removal the doublets ####################################################
#DoubletFinder
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)


##### Fix this ...

#pK identification
sweep.res.list_TCW <- paramSweep(TCW, PCs = 1:50, sct = FALSE)
sweep.stats_TCW <- summarizeSweep(sweep.res.list_TCW, GT = FALSE)
bcmvn_TCW <- find.pK(sweep.stats_TCW)

pdf(paste(dir, "/DoubletFinder_pk_param_ident_plot.pdf", sep = ""), width = 10, height = 10) 
ggplot(bcmvn_TCW, aes(pK, BCmetric, group = 1)) +
  geom_point() +
  geom_line()
dev.off()

pK <- bcmvn_TCW %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  select(pK) 
pK <- as.numeric(as.character(pK[[1]])) #optimal pK = 0.23

#pN = 0.25 (default)
#exp = 12% (high throughput) or 24% (standard) #USED 12

#Homotypic Doublet Estimation
annotations <- TCW@meta.data$SoupX_seurat_clustering_0.6
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.12*nrow(TCW@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#run doublet finder
TCW_doubletfinder <- doubletFinder_v3(TCW, 
                                      PCs = 1:50, 
                                      pN = 0.25, 
                                      pK = pK, 
                                      nExp = nExp_poi.adj,
                                      reuse.pANN = FALSE, sct = FALSE)

TCW$DF.classifications_0.25_0.23_841 <- TCW_doubletfinder$DF.classifications_0.25_0.23_841

DF_doublets <- subset(TCW_doubletfinder, DF.classifications_0.25_0.23_841 %in% "Doublet") #841




################################################################################
##################### Cell Type Annotation #####################################
################################################################################
## Read in organoid4

### Harmony Clusters

install.packages("harmony")
library(harmony)

## Run Harmony for Batch Correction:
TCW <- RunHarmony(TCW, group.by.vars = "sample")

# After Harmony integration, typically proceed with re-clustering
TCW <- FindNeighbors(TCW, reduction = "harmony")
TCW <- FindClusters(TCW)

# to visualize the results, you can run UMAP or another dimensionality reduction technique
TCW <- RunUMAP(TCW, reduction = "harmony")

DimPlot(TCW, reduction = "umap", group.by = "genotype")




TCW$harmony_clusters <- Idents(TCW)


Final_Markers_NPC = c("TOP2A", "MKI67", "FOXM1", "CENPF", "PAX6", "PCNA", "E2F2", "HOPX", "LHX2", "OTX2", "GLI3") 
Final_Markers_Neuron = c("DCX", "INSM1", "ST18", "NHLH1", "NEUROD1", "NEUROD6", "PRDM8", "MYT1L", "MAPT", "SNAP25", "SYP", "BCL11B", "DLG4", "STMN2", "SLC17A6", "GRIN2B", "RELN", "GRIA2", "GAD2", "DLX1", "DLX2", "SST", "NRXN1", "ANK3") #"DLX6-AS1" 
Final_Markers_Astrocyte = c("APOE", "GFAP", "PEA15", "S100B", "ALDH1L1", "FGFR3", "AGT", "AQP4")
Final_Markers_OPC = c("OLIG1", "OLIG2", "PCDH15", "LHFPL3") 
Final_Markers_VLMC = c("LUM", "PDGFRA", "DCN", "POSTN", "OGN", "APOD", "RSPO3") 
Final_Markers_Pericytes = c("CSPG4", "PDGFRB", "CD146", "VTN", "CSPG4", "ACTA2", "KCNJ8", "ABCC9", "ACE2", "ART3", "ATP13A5") #


## Visualize Canonical Cell Type Marker Expression
#Feature Plots

Final_Marker_List = list(Final_Markers_NPC, Final_Markers_Neuron, Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC)
Final_Marker_List_Names = c("NPC_Markers", "Neuron_Markers", "Astrocyte_Markers", "OPC_Markers", "VLMC_Markers")

count = 1
for (x in Final_Marker_List) {
  for (y in x) {
    
    name = Final_Marker_List_Names[count]
    pdf(file = paste(dir, "/Feature_Plots/", name, "_", y, "_feature_plot.pdf", sep = ""), width = 15, height = 15)
    print(FeaturePlot(TCW, reduction = "umap", features = y, label = TRUE, label.size = 6))
    dev.off()
  }
  count = count + 1
}

names(TCW@reductions)

# Feature plots should be similar to below:
#Dot Plots 

feats <- list(Final_Markers_NPC, Final_Markers_Neuron[1:13], Final_Markers_Neuron[14:25], Final_Markers_Astrocyte, Final_Markers_OPC, Final_Markers_VLMC)
marker_names <- list("NPC", "Neuron_1", "Neuron_2", "Astrocyte", "OPC", "VLMC")
count = 0

for (y in feats) {
  
  count = count + 1
  name = marker_names[count]
  
  #dotplot
  pdf(paste(dir, "/Dot_Plots/", name, "_dot_plot.pdf", sep = ""), width = 25, height = 15)
  
  print(DotPlot(
    TCW,
    features = y,
    cols = c("light grey", "blue"),
    col.min = -2.5,
    col.max = 2.5,
    dot.min = 0,
    dot.scale = 10,
    idents = NULL,
    group.by = NULL,
    split.by = NULL,
    cluster.idents = FALSE,
    scale = TRUE,
    scale.by = "radius"
    # scale.min = NA,
    # scale.max = NA
    )
  )
  
  dev.off()
} 


## Barcode rank plot

# Retrieve the total counts per cell
cell_counts <- TCW@meta.data$nCount_RNA

# Order cells by total counts
sorted_counts <- sort(cell_counts, decreasing = TRUE)

plot(log10(sorted_counts), pch = 19, cex = 0.5, main = "Barcode Rank Plot", 
                          xlab = "Cell Rank", ylab = "Log10 Total Counts")



### Find the markers of each cluster to confirm preliminary labels
#CSV of DEGs in Every Cluster

#Use FindMarkers to get the DEGs in each cluster
#Marker genes for each cluster are identified using the Wilcoxon rank sum test implemented in the FindMarkers function, 
#with a logFC threshold of 0.25 and Bonferroni correction for multiple testing
#Marker genes for each cluster were manually compared to known marker genes in the PangloDB Database
TCW$
#select clusters
Idents(TCW) <- TCW$harmony_clusters

for (x in rev(unique(TCW$harmony_clusters))){
  
  #get positive DEGs in each cluster
  cluster.markers <- FindMarkers(TCW, ident.1 = x, min.pct = 0.25, test.use = "wilcox", only.pos = TRUE)
  
  #csv file of marker genes
  write.csv(cluster.markers %>% slice_max(n=2000, order_by = avg_log2FC), paste(dir, "/Cluster_Markers/cluster_",x,"_markers.csv", sep = ""))
}


saveRDS(TCW,file.path(out,'TCW.rds'))

################################################################################
######## Annotate as a first stage with Single r ###############################
################################################################################

#SINGLER ANNOTATIONS: Automatically map query data to reference atlas

library(scRNAseq)
library(SingleR)
library(tidyverse)
library(Seurat)

remotes::install_version("SingleR", version = "1.0.1")

remotes::install_github('dviraran/SingleR')

#Reference Atlas: Polioudakis et. al.

ref_dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/ref_data"

#load count matrix
load(paste(ref_dir, "/sc_dev_cortex_geschwind/raw_counts_mat.rdata", sep = ""))
raw.counts.mat <- as.matrix(raw_counts_mat)
reffy <- CreateSeuratObject(raw_counts_mat)

#load metadata
ref_metadata <- read.csv(file = paste(ref_dir, "/sc_dev_cortex_geschwind/cell_metadata.csv", sep = ""))
rownames(ref_metadata) <- ref_metadata[,1]
ref_metadata[,1] <- NULL

#append metadata to reference object
reffy <- AddMetaData(reffy, ref_metadata)
reffy <- reffy[,!is.na(reffy$Cluster)]

#normalize and scale
reffy <- NormalizeData(reffy)
reffy <- ScaleData(reffy)

#convert to a single cell experiment
reffo <- as.SingleCellExperiment(reffy)
metadata(reffo)<- ref_metadata



#normalize data
library(scuttle)
reffo <- logNormCounts(reffo)
reffo@metadata$Cluster <- as.matrix(reffo@metadata$Cluster)

# Extract assay data from both test and ref objects
test_assay <- GetAssayData(TCW)
ref_assay <- assays(reffo)$counts



#perform SingleR labeling
results <- SingleR(test = test_assay, ref = ref_assay, labels = reffo$Cluster)

save(results, file = "TCW_Polioudakis_SingleR_Results.R")

#append resulting cell type labels to metadata
TCW$Polioudakis_SingleR <- results$labels

#save UMAP visualizations
pdf(file = paste(dir, "/SingleR/Polioudakis_Labelled_UMAP.pdf", sep = ""), width = 10, height = 10)
DimPlot(TCW, reduction = "harmony_res_0.5_umap", group.by = "Polioudakis_SingleR", label = TRUE, label.size = 6)
dev.off()

#save annotated seurat object
save(TCW, file = paste(dir, "/TCW_SingleR_Annotated.R", sep = ""))

########### Confirm the cells annotation based on both canonical and single r ##
###################### based  results ########################################## 
Idents(TCW) <- TCW$Final_Cell_Type 
unique(TCW$Final_Cell_Type)

newIdent <- "GPC"
names(newIdent) <- "NPC 2 (Cycling NPC + GPC)"
TCW <- RenameIdents(object = TCW, newIdent)

###### Quantify Cell Type Population Proportions in Dataset ####################
################################################################################


#Barplots
#BiocManager::install("dittoSeq")
library(dittoSeq)

#Color Palettes!
APOE_colors = c("#1984c5", "#c23728")
Individual_colors = c("grey", "#bd7ebe")

#retrieve harmony FCT colors
require(scales)
identities <- levels(TCW$Final_Cell_Type)
FCT_palette <- hue_pal()(length(identities))
FCT_order <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)

#official id retrieve colors
require(scales)
identities <- unique(TCW$Individual)
Ind_color_palette <- hue_pal()(length(identities))

#official id retrieve order
Ind_cluster_order <- match(unique(TCW@meta.data[["Individual"]]), metaLevels("Individual", TCW))

###

#APOE distribution
png(file = paste(dir, "/Bar_Plots/Barplot_APOE_Genotype_Distribution.png", sep = "") , width = 5, height = 5) 
dittoBarPlot(object = TCW,
             var = "APOE_Genotype", 
             group.by = "orig.ident",
             y.breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
             xlab = "TCW_14311_run1",
             color.panel = APOE_colors) + ggtitle("Barplot of APOE Genotype Distribution")
dev.off()

#APOE distribution split by cell type
png(file = paste(dir, "/Bar_Plots/Barplot_APOE_Genotype_Distribution_by_FCT.png", sep = "") , width = 5, height = 5)
dittoBarPlot(object = TCW,
             var = "Final_Cell_Type", 
             group.by = "APOE_Genotype",
             color.panel = FCT_colors,
             var.labels.reorder = FCT_order,
             y.breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +ggtitle(expression(paste("Cell Subtype Proportion Changes in ", italic('APOE'), " 44")))
dev.off()

#Cell Type split by APOE Genotype
png(file = paste(dir, "/Bar_Plots/Barplot_Cell_Type_Distribution_by_APOE_Genotype.png", sep = "") , width = 6, height = 5)
dittoBarPlot(object = TCW,
             var = "APOE_Genotype", 
             group.by = "Final_Cell_Type",
             color.panel = APOE_colors,
             y.breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +ggtitle("Barplot of APOE Genotype Distribution Across Cell Types")
dev.off()

##

#Cell Type split by Individual
png(file = paste(dir, "/Bar_Plots/Barplot_Individual_Distribution_by_Cell_Type.png", sep = "") , width = 6, height = 5)
dittoBarPlot(object = TCW,
             var = "Individual", 
             group.by = "Final_Cell_Type",
             color.panel = Individual_colors,
             y.breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +ggtitle("Barplot of Individual Distribution Across Cell Types")
dev.off()

















#### First, do we need to correct for ambient RNA?                          ####
#### Generate a knee plot of the barcode ranks
# Normalize the data if not already done
organoid3 <- NormalizeData(organoid2, normalization.method = "LogNormalize", scale.factor = 10000)

# Retrieve the total counts per cell
cell_counts <- organoid3@meta.data$nCount_RNA

# Order cells by total counts
sorted_counts <- sort(cell_counts, decreasing = TRUE)

# Create a simple Barcode Rank plot
plot(log10(sorted_counts), pch = 19, cex = 0.5, main = "Barcode Rank Plot", xlab = "Cell Rank", ylab = "Log10 Total Counts")

##### Let's produce a tSNE projection plot in order to identify which features
#### Are appearing as enriched

library(patchwork)

# Finding variable features
organoid3 <- FindVariableFeatures(organoid3, selection.method = "vst", nfeatures = 2000)
# Plotting variable features
top10 <- head(VariableFeatures(organoid3), 10)
plot1 <- VariableFeaturePlot(organoid3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
combined_plot <- plot1 + plot2
print(combined_plot)

# Scale the data
organoid3 <- ScaleData(organoid3, features = rownames(organoid3))

# PCA
organoid3 <- RunPCA(organoid3, features = VariableFeatures(object = organoid3))
DimPlot(organoid3, reduction = "pca")
# t-SNE
# First, ensure that PCA has been performed as t-SNE relies on the results of PCA
organoid3 <- RunPCA(organoid3, features = VariableFeatures(object = organoid3))

# Now run t-SNE
organoid3 <- RunTSNE(organoid3, dims = 1:10) # You can adjust the number of dims if needed

# After running t-SNE, you can plot the top genes
# First, identify the top features contributing to PC1
top_genes_pc1 <- head(VariableFeatures(organoid3), 10)

# Now create a t-SNE plot with expression levels of the top genes from PC1
FeaturePlot(organoid3, features = top_genes_pc1, reduction = "tsne")

# Create a standard t-SNE plot
tsne_plot <- DimPlot(organoid3, reduction = "tsne")

# Overlay labels for specific genes
tsne_plot <- LabelPoints(plot = tsne_plot, points = top_genes_pc1, repel = TRUE)

# Print the plot with labels
print(tsne_plot)



organoid3 <- RunTSNE(organoid3, dims = 1:10)
DimPlot(organoid3, reduction = "tsne")
# UMAP
organoid3 <- RunUMAP(organoid3, dims = 1:10)
DimPlot(organoid3, reduction = "umap")


# Finding clusters
organoid3 <- FindNeighbors(organoid3, dims = 1:10)
organoid3 <- FindClusters(organoid3, resolution = 0.5)
# Plotting clusters
DimPlot(organoid3, reduction = "umap", label = TRUE)


###### Ambient RNA removal via DropletQC pkg ##################################
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# devtools::install_github("powellgenomicslab/DropletQC")
library(DropletQC)

# Access the raw counts of the Seurat object
raw_counts <- seurat_obj@assays$RNA@counts



























################################################################################
###### Step 2: Perform QC across all samples ###################################
###### then visualize the results             ##################################
################################################################################
################################################################################
####### Loop to perform QC on each individual Seurat object ####################
################################################################################

setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')
## Begin step 2 by reading in our organoid0.rds object:
organoid0 <- readRDS('outputs/simpleafSeurat/organoid0.rds')
# List of sample identifiers
samples <- c("1", "2", "3", "4", "6", "8", "9", "10",
             "11", "12", "13", "14", "15", "16", "17", "18", "19",
             "20", "21", "22", "23", "24", "25", "26")



library(Seurat)
library(dplyr)

###* Step 2: Perform QC on each individual Seurat object, then combine them *###
###* and plot the results to compare to organoid0.rds                       *###

# Initialize a data frame to store cell counts
cell_counts <- data.frame(sample = character(),
                          cells_before_QC = integer(),
                          cells_after_QC = integer(),
                          stringsAsFactors = FALSE)

# Initialize variables for combined statistics
stats <- list(
  APOE22 = list(cells_before = 0, cells_after = 0, genes_before = vector(), genes_after = vector()),
  APOE33 = list(cells_before = 0, cells_after = 0, genes_before = vector(), genes_after = vector()),
  APOE33Ch = list(cells_before = 0, cells_after = 0, genes_before = vector(), genes_after = vector()),
  APOE44 = list(cells_before = 0, cells_after = 0, genes_before = vector(), genes_after = vector())
)

# APOE22: S3, S9, S12, S15, S21, S25
# APOE33: S1, S8, S14, S16, S19, S23
# APOE33Ch: S2, S6, S13, S18, S20, S24
# APOE44: S4, S10, S11, S17, S22, S26 *##

APOE22_samples_group <- c("S3", "S9", "S12", "S15", "S21", "S25")
APOE33_samples_group <- c("S1", "S8", "S14", "S16", "S19", "S23")
APOE33Ch_samples_group <- c("S2", "S6", "S13", "S18", "S20", "S24")
APOE44_samples_group <- c("S4", "S10", "S11", "S17", "S22", "S26")

sample_names_group22 <- c("sample3", "sample9", "sample12", "sample15", "sample21", "sample25")
sample_names_group33 <- c("sample1", "sample8", "sample14", "sample16", "sample19", "sample23")
sample_names_group33Ch <- c("sample2", "sample6", "sample13", "sample18", "sample20", "sample24")
sample_names_group44 <- c("sample4", "sample10", "sample11", "sample17", "sample22", "sample26")

# Loop through each sample and add a metadata column for sample labels
# for APOE22
for (i in seq_along(APOE22_samples_group)) {
  object_name <- APOE22_samples_group[i]
  sample_name <- sample_names_group22[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE22'")))
}

for (i in seq_along(APOE33_samples_group)) {
  object_name <- APOE33_samples_group[i]
  sample_name <- sample_names_group33[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33'")))
}

for (i in seq_along(APOE33Ch_samples_group)) {
  object_name <- APOE33Ch_samples_group[i]
  sample_name <- sample_names_group33Ch[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE33Ch'")))
}

for (i in seq_along(APOE44_samples_group)) {
  object_name <- APOE44_samples_group[i]
  sample_name <- sample_names_group44[i]
  
  eval(parse(text = paste0(object_name, "[['sample']] <- '", sample_name, "'")))
  eval(parse(text = paste0(object_name, "[['genotype']] <- 'APOE44'")))
}

out <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/simpleafSeurat'
unique(organoid3@meta.data$sample)
# Loop over each sample
for (sample_id in samples) {
  # Load the Seurat object
  T_QC <- readRDS(file = paste0("./outputs/simpleafSeurat/sample", sample_id, ".rds"))
  
  # Record the number of cells before QC
  cells_before <- ncol(T_QC)
  
  # Filter genes expressed in less than 3 cells
  # Using sparse matrix directly to avoid memory issues
  raw_counts_mat <- T_QC@assays$RNA@counts
  genes_to_keep <- Matrix::colSums(raw_counts_mat != 0) >= 3
  T_QC <- subset(T_QC, features = rownames(T_QC)[genes_to_keep])
  
  # Set Idents, calculate MT and Ribosomal percentages
  Idents(T_QC) <- T_QC$orig.ident
  # Subset mt_gene_names to include only genes present in the Seurat object
  valid_mt_gene_names <- mt_gene_names[mt_gene_names %in% rownames(T_QC[["RNA"]]@counts)]
  
  # Recalculate the percentage of mitochondrial genes
  T_QC[["percent.mt"]] <- PercentageFeatureSet(T_QC, features = valid_mt_gene_names)
  
  # Filter out cells with greater than 5% mtDNA
  T_QC <- subset(T_QC, subset = percent.mt <= 5)
  
  # Assign percent.ribo to cells
  T_QC[["percent.ribo"]] <- PercentageFeatureSet(T_QC, pattern = "^RP[SL]")
  
  # Calculate upper limit, set at 99.5 percentile
  upper_limit <- quantile(T_QC@meta.data$nCount_RNA, probs = 0.995)
  
  # Calculate Barcode Inflections for lower limit
  T_QC <- CalculateBarcodeInflections(
    T_QC,
    barcode.column = "nCount_RNA",
    group.column = "orig.ident",
    threshold.low = 1000,
    threshold.high = NULL
  )
  T_QC <- SubsetByBarcodeInflections(object = T_QC)
  
  # Filter cells based on RNA count limits
  T_QC <- subset(T_QC, subset = nCount_RNA < upper_limit)
  
  # Record the number of cells after QC
  cells_after <- ncol(T_QC)
  
  #save modified object
  filename <- paste0("QC_sample", sample_id, ".rds")
  saveRDS(T_QC, file.path(out, filename))
  
  # Append the cell count information to the data frame
  cell_counts <- rbind(cell_counts, data.frame(sample = sample_id,
                                               cells_before_QC = cells_before,
                                               cells_after_QC = cells_after))
}

# Save the cell counts data frame
output_dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/"
write.csv(cell_counts, file = paste0(output_dir, "cell_counts.csv"), row.names = FALSE)
# Convert stats to a data frame
data_frame_stats <- data.frame(APOE_group = names(stats),
                               cells_before = sapply(stats, function(x) x$cells_before),
                               cells_after = sapply(stats, function(x) x$cells_after),
                               genes_before = sapply(stats, function(x) length(x$genes_before)),
                               genes_after = sapply(stats, function(x) length(x$genes_after)))

# Save the data frame
write.csv(data_frame_stats, file = paste0(output_dir, "APOE_group_stats.csv"), row.names = FALSE)

S1 <- readRDS('outputs/simpleafSeurat/QC_sample1.rds')
S2 <- readRDS('outputs/simpleafSeurat/QC_sample2.rds')
S3 <- readRDS('outputs/simpleafSeurat/QC_sample3.rds')
S4 <- readRDS('outputs/simpleafSeurat/QC_sample4.rds')
S6 <- readRDS('outputs/simpleafSeurat/QC_sample6.rds')
S8 <- readRDS('outputs/simpleafSeurat/QC_sample8.rds')
S9 <- readRDS('outputs/simpleafSeurat/QC_sample9.rds')
S10 <- readRDS('outputs/simpleafSeurat/QC_sample10.rds')
S11 <- readRDS('outputs/simpleafSeurat/QC_sample11.rds')
S12 <- readRDS('outputs/simpleafSeurat/QC_sample12.rds')
S13 <- readRDS('outputs/simpleafSeurat/QC_sample13.rds')
S14 <- readRDS('outputs/simpleafSeurat/QC_sample14.rds')
S15 <- readRDS('outputs/simpleafSeurat/QC_sample15.rds')
S16 <- readRDS('outputs/simpleafSeurat/QC_sample16.rds')
S17 <- readRDS('outputs/simpleafSeurat/QC_sample17.rds')
S18 <- readRDS('outputs/simpleafSeurat/QC_sample18.rds')
S19 <- readRDS('outputs/simpleafSeurat/QC_sample19.rds')
S20 <- readRDS('outputs/simpleafSeurat/QC_sample20.rds')
S21 <- readRDS('outputs/simpleafSeurat/QC_sample21.rds')
S22 <- readRDS('outputs/simpleafSeurat/QC_sample22.rds')
S23 <- readRDS('outputs/simpleafSeurat/QC_sample23.rds')
S24 <- readRDS('outputs/simpleafSeurat/QC_sample24.rds')
S25 <- readRDS('outputs/simpleafSeurat/QC_sample25.rds')
S26 <- readRDS('outputs/simpleafSeurat/QC_sample26.rds')

###* At this point, let's just combine our QC'd samples into an .rds called *###
###* organoid1.rds and save it                                              *###

organoid1 <- merge(S1, y = c(S2, S3, S4, S6, S8, S9, S10,
                             S11, S12, S13, S14, S15, S16, S17, S18, S19,
                             S20, S21, S22, S23, S24, S25, S26), 
                   add.cell.ids = c("S1", "S2", "S3", "S4", "S6", "S8", "S9", "S10",
                                    "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19",
                                    "S20", "S21", "S22", "S23", "S24", "S25", "S26"),
                   project = "APOE_Jorganoid"
)


out <- '/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/'
saveRDS(organoid1,file.path(out,'organoid1.rds'))

#organoid1 <- readRDS('outputs/simpleafSeurat/organoid1.rds')

# Assuming organoid1 is your Seurat object and 'genotype' is a column in meta.data
unique_genotypes <- unique(organoid1@meta.data$genotype)

# Output directory
out <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/images/"  # Replace with your actual directory path

# Define the features to plot
features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

# Loop through each genotype
for (genotype in unique_genotypes) {
  # Subset Seurat object by genotype
  organoid_subset <- subset(organoid1, subset = genotype == genotype)
  
  # Create violin plots
  vps <- VlnPlot(object = organoid_subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
  
  # Modify the plots to spread them out evenly
  vps <- vps + plot_layout(guides = 'collect') & theme(legend.position = 'none')
  
  # Add the genotype as the overall figure title
  combined_plot <- vps + plot_annotation(title = genotype, theme = theme(plot.title = element_text(hjust = 0.5)))
  
  # Construct the full path for the output file
  output_file_path <- file.path(out, paste0("after_QC_violin_plots_", genotype, ".png"))
  
  # Save the combined plot to the specified directory
  ggsave(output_file_path, plot = combined_plot, width = 12, height = 6, units = "in")
}

################################################################################
################################################################################
########## Step 3: Dimensionality reduction and analysis #######################
################################################################################
################################################################################

###* Begin step 3 by reading in our organoid1 quality controlled object     *###
###* This is a good time to restart your R session.   
organoid1 <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/organoid1.rds')
library(Seurat)
library(miQC)

###* don't over-saturate your cache. Easier to rename all of these           *###
# Assuming organoid1 is your Seurat object
T_QC <- organoid1

# Preprocessing steps
T_QC <- NormalizeData(T_QC, verbose = FALSE)
T_QC <- FindVariableFeatures(T_QC, selection.method = "vst", nfeatures = 2000)
T_QC <- ScaleData(T_QC, features = rownames(T_QC), verbose = FALSE)

# Dimensionality reduction
T_QC <- RunPCA(T_QC, features = VariableFeatures(object = T_QC), verbose = FALSE)
T_QC <- RunUMAP(T_QC, dims = 1:50)

# Clustering
T_QC <- FindNeighbors(T_QC, dims = 1:50, k.param = 30, verbose = FALSE)
T_QC <- FindClusters(T_QC, resolution = 0.6, verbose = FALSE)
###* Show in each cluster the genotype (extract from organoid2@meta.data$geno...*###
###* Plot what is each genotype distribution in each cluster                *###
###* She was using findCluster, did not do actual cellType annotation here. *###
###* When you do UMAP, you can color after by Seurat cluster (must run her other code first) *###
###* 
###* Bring questions for Deepti for the meeting tomorrow
# Visualization
DimPlot(T_QC, reduction = "umap", group.by = "genotype")

###* this is a good point to save our .rds object, this time as organoid2   *###
saveRDS(T_QC,file.path(out,'organoid2.rds'))

# ribo_counts <- as.data.frame(APOE22$percent.ribo)
# colnames(ribo_counts) <- "percent.ribo" # Rename the column without using dplyr
organoid2 <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/simpleafSeurat/organoid2.rds')
organoid2 <- T_QC

#DimPlot(organoid2, reduction = "umap", group.by = "genotype")
###* Let's take a look at the percent.RNA of the cells in the data          *###
RNA_counts <- as.data.frame(organoid2$percent.RNA)
colnames(RNA_counts) <- "percent.RNA"
# Now you can plot without renaming
ggplot(RNA_counts, aes(x = percent.RNA)) + 
  geom_density() +
  theme_bw() +
  ggtitle("Percent RNA genes per cell")

# PCA Visualization
DimPlot(organoid2, reduction = "pca")

# Adjust UMAP dimensions
organoid2 <- RunUMAP(organoid2, dims = 1:15)
DimPlot(organoid2, reduction = "umap")

# Clustering
organoid2 <- FindNeighbors(organoid2, dims = 1:15)
organoid2 <- FindClusters(organoid2)

# Quality Control Feature Plots
plot1 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file = "qcfeatureplots_filtered_object.pdf", width = 15, height = 10)
plot1 + plot2
dev.off()


#Remove Additional Outliers (added 1/5/23)
#Only include cells with at least 200 genes expressed
T_QC <- organoid2
T_QC <- RunUMAP(organoid2, dims = 1:50)
library(miQC)
T_QC.MiQC <- T_QC
T_QC.MiQC <- subset(T_QC, miQC.keep == "keep") 
T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, k.param = 30, verbose = FALSE)
T_QC.MiQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected

#umap
T_QC.MiQC <- RunUMAP(T_QC.MiQC, assay = SCT, dims = 1:50,  n.neighbors = 30, verbose = FALSE)
#7,337 cells and 22,719 genes

#Check violin plots after outlier removal
Idents(object = T_QC.MiQC) <- "project"
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/violin_plot.pdf", width = 15, height = 10)
VlnPlot(T_QC.MiQC, features = c("nFeature_RNA", "nCount_RNA"), slot = "counts", pt.size = 0, split.by = NULL, ncol = 3) 
dev.off()

####################################################################################################################
#Re-cluster dataset after filtering
T_QC@meta.data$percent.mt
#normalize data: SCT
all.genes <- rownames(T_QC)
T_QC <- ScaleData(T_QC, features = all.genes, assay = "RNA", vars.to.regress = c("percent.mt", "CC.Difference"), verbose = FALSE)
T_QC <- NormalizeData(T_QC, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
T_QC <- FindVariableFeatures(T_QC, selection.method = "vst", nfeatures = 2000)

#SCT assay
#scale data while removing variation due to mt percentage & cell phase
T_QC <-SCTransform(T_QC, vars.to.regress = c("percent.mt", "CC.Difference"))


#find variable features: find the 2000 most variable genes
T_QC.MiQC <- FindVariableFeatures(T_QC.MiQC, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(T_QC.MiQC), 10)
?VariableFeaturePlot
plot1 <- VariableFeaturePlot(T_QC, selection.method = "vst")
plot1
head(VariableFeatures(T_QC.MiQC))

#saveRDS(T_QC.miQC, file.path(out, "T_QC.miQC"))

plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/variablefeatureplot_filtered_object.pdf", width = 15, height = 10)
plot1+plot2
dev.off()

#scale data while removing variation due to mt percentage
all.genes <- rownames(T_QC.MiQC)
T_QC.MiQC <- ScaleData(T_QC.MiQC, features = all.genes, vars.to.regress = "percent.mt")

#pca
T_QC.MiQC <- RunPCA(T_QC.MiQC, features = VariableFeatures(object = T_QC.MiQC), assay = "SCT")

#examine & visualize pca results
print(T_QC[["pca"]], dims = 1:5, nfeatures = 5)
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/pc1&pc2_filtered_object.pdf", width = 15, height = 20)
VizDimLoadings(T_QC, dims = 1:2, reduction = "pca")
dev.off()
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/pc1&pc2dimplot_filtered_object.pdf", width = 15, height = 20)
DimPlot(T_QC.MiQC, reduction = "pca")
dev.off()

###################################
#3D pca plot
#date: 1/7/23
#source: https://rpubs.com/HWH/920093
library(plotly)

# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)


# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = T_QC.MiQC, vars = c("PC_1", "PC_2", "PC_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

seurat_default_colors <- c("#F8766D", "#E68613", "#CD9600", "#ABA300", "#7CAE00", "#0CB702", "#00BE67",
                           "#00C19A", "#00BFC4", "#00B8E7", "#00A9FF", "#8494FF", "#C77CFF", "#ED68ED", "#FF61CC", "#FF68A1")

fig <- plot_ly(data = plot.data, 
               x = ~PC_1, y = ~PC_2, z = ~PC_3, 
               color = ~seurat_clusters, 
               colors = seurat_default_colors,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig_cube

#########################################
#pca heatmap #downsampled to 500 cells
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/pc1_heatmap_filtered_object.pdf")
DimHeatmap(object = T_QC.MiQC, dims = 1, cells = 500, reduction = 'pca', balanced = TRUE) 
?DimHeatmap
#+ggplot2::scale_fill_gradientn(colors = c("steelblue3", "white", "tomato2"))
dev.off()
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/pcs1-15_heatmap_filtered_object.pdf")
DimHeatmap(T_QC.MiQC, dims = 1:15, cells = 500, balanced = TRUE)
#+ggplot2::scale_fill_gradientn(colors = c("steelblue1", "white", "tomato"))
dev.off()

#pca: elbow plot
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/pca_elbow_plot.pdf", width = 30, height = 30)
ElbowPlot(T_QC.MiQC)
dev.off()

#cluster the cells based on previously identified PCs
#KNN: K-nearest neighbor network: Graph-Theoretic Approach
T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, verbose = FALSE)
T_QC.MiQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected
head(Idents(T_QC.MiQC), 5)
#22 clusters

#umap
T_QC.MiQC <- RunUMAP(T_QC.MiQC, assay = SCT, dims = 1:50, verbose = FALSE)

pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/umap_bycluster_filtered_object.pdf", width = 10, height = 10)
DimPlot(T_QC.MiQC, reduction = "umap", group.by = "seurat_clusters", pt.size = .1) +ggtitle("UMAP by Seurat Cluster") + annotate("text", x=8.5, y=-20, label= "Cells = 7,300")
dev.off()

save(T_QC.MiQC, file = "T_QC.MiQC")
load(file = "T_QC.MiQC")

T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, verbose = FALSE)
T_QC.MIQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected
head(Idents(T_QC.MC), 5)
#res = 0.2: 16 clusters
T_QC.MIQC <- RunUMAP(T_QC.MIQC, assay = SCT, dims = 1:50, verbose = FALSE)
DimPlot(T_QC.MIQC, reduction = "umap", group.by = "seurat_clusters", pt.size = .1)


# Load necessary libraries
library(ggplot2)

# Assuming you have loaded organoid0 and organoid1 from your .rds files


# Assuming organoid0 and organoid1 are loaded from your .rds files

# Extract and check percent.mt data for organoid0
if ("percent.mt" %in% names(organoid0@meta.data)) {
  organoid0_percent_mt <- organoid0@meta.data$percent.mt
  if (length(organoid0_percent_mt) > 0) {
    organoid0_df <- data.frame(percent_mt = organoid0_percent_mt, organoid = 'organoid0')
  } else {
    warning("percent.mt in organoid0 is empty.")
  }
} else {
  warning("percent.mt not found in organoid0 meta data.")
}

# Extract and check percent.mt data for organoid1
if ("percent.mt" %in% names(organoid1@meta.data)) {
  organoid1_percent_mt <- organoid1@meta.data$percent.mt
  if (length(organoid1_percent_mt) > 0 && organoid1_percent_mt <= 5) {
    organoid1_df <- data.frame(percent_mt = organoid1_percent_mt, organoid = 'organoid1')
  } else {
    warning("percent.mt in organoid1 is empty.")
  }
} else {
  warning("percent.mt not found in organoid1 meta data.")
}

# Combine the data (if both are available)
if (exists("organoid0_df") && exists("organoid1_df")) {
  combined_df <- rbind(organoid0_df, organoid1_df)
} else {
  warning("One or both of organoid data frames are not available.")
}
# Create a density plot
ggplot(combined_df, aes(x = percent_mt, fill = organoid, color = organoid)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Mitochondrial DNA per Cell",
       x = "Percent Mitochondrial DNA",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("blue", "red")) +
  scale_color_manual(values = c("blue", "red"))

organoid1 <- subset(organoid1)
organoid1 <- subset(organoid1, subset = percent.mt <= 5)

################################################################################
################################################################################
###### Deepti's original Pipeline: #############################################


out <- "outputs/simpleaf"
saveRDS(APOE22,file.path(out,'APOE22.rds'))



unique_values_genotype <- unique(organoid1@meta.data$genotype)
print(unique_values_genotype)
# Visualize the PCA, grouping by cell cycle phase
DimPlot(organoid1,
        reduction = "pca",
        group.by= "genotype")

organoid2@meta.data$genotype
# Visualize PCA, grouping by Individual 
DimPlot(T_QC,
        reduction = "pca",
        group.by= "individual")

# Visualize PCA, grouping by APOE Genotype 
DimPlot(T_QC,
        reduction = "pca",
        group.by= "genotype")

# Visualize PCA, grouping by APOE Genotype 
DimPlot(T_QC,
        reduction = "pca",
        group.by= "A")










### end of changes by akg ###

#LOAD ORIGINAL TCW OBJECT: 'seurat_obj.rds'
T_QC <- readRDS(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Data/Project_TCW_14311_B01_SCR_Lane.2020-07-20/compbio/Results_to_deliver/seurat_object.rds") 
T_QC <-  UpdateSeuratObject(object = T_QC)
#SCT assay: 7,911 cells & 46423 genes
#22 Seurat clusters

#EDIT OBJECT: RE-LABEL SAMPLES, LABEL BY APOE GENOTYPE
#set project name
T_QC@project.name <- "TCW_14311_run1"

#Add Metadata Label by APOE Genotype
#Split Dataset by APOE Genotype (APOE 33 versus APOE 44)
E33_list = c("ApoE33-1-GTCAACTCTTTAGCG", "ApoE33-2-TTCCGCCTCTCTTTG", 
             "ApoE33.MG-2-TGTCTTTCCTGCCAG", "ApoE33.MG-1-AAGTATCGTTTCGCA") 

E44_list = c("ApoE44.MG-2-CTCCTCTGCAATTAC", "ApoE44-1-TGATGGCCTATTGGG",  
             "ApoE44.MG-1-GGTTGCCAGATGTCA",  "ApoE44-2-AGTAAGTTCAGCGTA")

T_QC.E33 <- subset(T_QC, subset = hash.ID %in% E33_list)
T_QC.E44 <- subset(T_QC, subset = hash.ID %in% E44_list)

cellNames <- rownames(T_QC.E33@meta.data)

T_QC$barcode <- rownames(T_QC@meta.data)
T_QC@meta.data <- T_QC@meta.data %>% mutate(APOE_Genotype = ifelse((T_QC$barcode %in% cellNames), "APOE 33",  "APOE 44"))

#individual metadata
individual_1_list = c("ApoE33-1-GTCAACTCTTTAGCG", "ApoE44-1-TGATGGCCTATTGGG", 
                      "ApoE44.MG-1-GGTTGCCAGATGTCA", "ApoE33.MG-1-AAGTATCGTTTCGCA") 

individual_2_list = c("ApoE44.MG-2-CTCCTCTGCAATTAC", "ApoE33-2-TTCCGCCTCTCTTTG",  
                      "ApoE33.MG-2-TGTCTTTCCTGCCAG", "ApoE44-2-AGTAAGTTCAGCGTA")

T_QC.1 <- subset(T_QC, subset = hash.ID %in% individual_1_list)
T_QC.2 <- subset(T_QC, subset = hash.ID %in% individual_2_list)

cellNames <- rownames(T_QC.1@meta.data)

T_QC@meta.data <- T_QC@meta.data %>% mutate(individual = ifelse((T_QC$barcode %in% cellNames), "individual 1",  "individual 2"))

#check metadata
#matrix <- T_QC@meta.data

#save seurat object
save(T_QC, file = "T_QC")

#no batch effects, all run together
#unique(levels(T_QC$HTO_classification))#8, 37 levels
#unique(levels(T_QC$hash.ID)) #8, 10 levels

#Now, APOE 33 cells are labeled with the metadata "APOE_Genotype" as "APOE 33" and APOE 44 cells are labeled as "APOE 44".
##########

#Initial Data Distribution Before Any QC

metadata <- as.data.frame(organoid2@meta.data)

#APOE Genotype
pdf(paste(dir, "/original/original_APOE_Distribution_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
ggplot(metadata, aes(x = genotype)) + geom_bar() + theme_bw()
dev.off()

#Individual
pdf(paste(dir, "/original/original_Individual_Distribution_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
ggplot(metadata, aes(x = individual)) + geom_bar() + theme_bw()
dev.off()

#Sample ID
pdf(paste(dir, "/original/original_Sample_Distribution_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
ggplot(metadata, aes(x = hash.ID)) + geom_bar() + theme_bw() + theme(axis.text.x = element_text(angle = 90))
dev.off()

#number of cells each gene is expressed in
raw_counts_mat <- as.data.frame(T_QC@assays$RNA@counts)
raw_counts_mat <- raw_counts_mat %>% mutate(num_cells_expressing =rowSums(.!=0)) #count non-zero cells

pdf(paste(dir, "/original/original_number_of_cells_expressing_each_gene_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
ggplot(raw_counts_mat, aes(num_cells_expressing)) + geom_density()
dev.off()

#number of genes each cell expresses

pdf(paste(dir, "/original/original_number_of_genes_expressed_by_each_cell_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
ggplot(metadata, aes(nFeature_RNA)) + geom_density()
dev.off()

#save file
save(T_QC, file = paste(dir,"/original/T_QC", sep = ''))

####################################################################################################################
#New QC 8/2023

load(file = "T_QC")

#We start with just the raw umi count data.
head(T_QC) #46423 features & 7911 samples
T_QC@active.assay <- "RNA" #23,696 genes & 7,911 cells

#Convert to a SingleCellExperiment
#T_QC <- as.SingleCellExperiment(T_QC)

#####################################################
#Remove Droplets with only ambient RNA

#Seurat Barcode Rank
T_QC <- CalculateBarcodeInflections(
  T_QC,
  barcode.column = "nCount_RNA",
  group.column = "orig.ident",
  threshold.low = NULL,
  threshold.high = NULL
)

pdf(paste(dir, "/Empty_Droplet_Removal/seurat_barcode_ranks_plot.pdf", sep = ""), width = 10, height = 10)
BarcodeInflectionsPlot(T_QC)
dev.off()
#DropletUtils

#Barcode Rank Plot: log(- total UMI) vs. log ( - rank) [ranked by total UMI in decreasing order]
br.out <- barcodeRanks(raw_counts_mat)

pdf(paste(dir, "/Empty_Droplet_Removal/dropletutils_barcode_ranks_plot.pdf", sep = ""), width = 10, height = 10)
plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
o <- order(br.out$rank)
lines(br.out$rank[o], br.out$fitted[o], col="red")

abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
       legend=c("knee", "inflection"))
dev.off()


#emptyDrops function
?? emptyDrops 

set.seed(100)
empty_drops <- emptyDrops(raw_counts_mat)

#T_QC <- 




#####################################################

###
#first, keep only the genes expressed in at least 3 cells
raw_counts_mat <- as.data.frame(T_QC@assays$RNA@counts)
raw_counts_mat <- raw_counts_mat %>% mutate(num_cells_expressing =rowSums(.!=0)) #count non-zero cells

pdf(paste(dir, "/QC_number_of_cells_expressing_each_gene_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
ggplot(raw_counts_mat, aes(num_cells_expressing)) + geom_density()
dev.off()

raw_counts_mat <- subset(raw_counts_mat, num_cells_expressing >= 3)
geneNames <- rownames(raw_counts_mat) #genes in >= 3 cells
length(geneNames) #23,412 genes

T_QC <- subset(T_QC, features = geneNames) 
#23,412 genes & 7,911 cells
save(T_QC, file = "T_QC")
###

raw_counts_mat <- as.data.frame(T_QC@assays$RNA@counts)
raw_counts_mat <- raw_counts_mat %>% mutate(num_cells_expressing =rowSums(.!=0)) #count non-zero cells

pdf(paste(dir, "/QC_number_of_cells_expressing_each_gene_after_cutoffs.pdf", sep = ""), width = 10, height = 10)
ggplot(raw_counts_mat, aes(num_cells_expressing)) + geom_density()
dev.off()
###

#total number of reads (UMIs) detected per cell = library size
#density plot
pdf(paste(dir, "density_plot_nUMI_per_cell_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
umi_counts <- as.data.frame(T_QC$nCount_RNA)
umi_counts %>%
  rename("nCount_RNA" = "T_QC$nCount_RNA") %>%
  ggplot(aes(x = nCount_RNA)) + 
  geom_density() +
  theme_bw()
dev.off()

#barplot
pdf(paste(dir, "bar_plot_nUMI_per_cell_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
umi_counts <- as.data.frame(T_QC$nCount_RNA)
umi_counts <- umi_counts %>%
  rename("nCount_RNA" = "T_QC$nCount_RNA") %>%
  rownames_to_column("ID") %>%
  add_column(seq(1, nrow(umi_counts)))  %>%
  rename("cell_index" = "seq(1, nrow(umi_counts))")
ggplot(umi_counts, aes(x = cell_index, y = nCount_RNA)) + geom_bar(stat='identity')+
  theme_bw()
dev.off()

#total number of genes detected per cell
pdf(paste(dir, "density_plot_genes_per_cell_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
gene_counts <- as.data.frame(T_QC$nFeature_RNA)
gene_counts %>%
  rename("nFeature_RNA" = "T_QC$nFeature_RNA") %>%
  ggplot(aes(x = nFeature_RNA)) + 
  geom_density() +
  theme_bw()
dev.off()

#barplot
pdf(paste(dir, "bar_plot_genes_per_cell_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
gene_counts <- as.data.frame(T_QC$nFeature_RNA)
gene_counts <- gene_counts %>%
  rename("nFeature_RNA" = "T_QC$nFeature_RNA") %>%
  rownames_to_column("ID") %>%
  add_column(seq(1, nrow(gene_counts)))  %>%
  rename("cell_index" = "seq(1, nrow(gene_counts))") 
ggplot(gene_counts, aes(x = cell_index, y = nFeature_RNA)) + geom_bar(stat='identity')+
  theme_bw()
dev.off()

#Set Idents
Idents(T_QC) <- T_QC$orig.ident

#label MT percent in each cell
T_QC[["percent.mt"]] <- PercentageFeatureSet(T_QC, pattern = "^MT-")

#MT percentage per cell

pdf(paste(dir, "density_plot_mtc_percent_per_cell_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
mt_counts <- as.data.frame(T_QC$percent.mt)
mt_counts %>%
  rename("percent.mt" = "T_QC$percent.mt") %>%
  ggplot(aes(x = percent.mt)) + 
  geom_density() +
  theme_bw()
dev.off()

#Ribosomal percent per cell
T_QC[["percent.ribo"]] <- PercentageFeatureSet(T_QC, pattern =  "^RP[SL]")

pdf(paste(dir, "density_plot_ribo_percent_per_cell_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
ribo_counts <- as.data.frame(T_QC$percent.ribo)
ribo_counts %>%
  rename("percent.ribo" = "T_QC$percent.ribo") %>%
  ggplot(aes(x = percent.ribo)) + 
  geom_density() +
  theme_bw()
dev.off()

###
#QC metric violin plots
vps<-VlnPlot(object = organoid0, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0)
vps
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="red") +labs(x = "TCW_14311_run1")
vp2<-vps[[2]]+geom_hline(yintercept = 75000,colour="red") +geom_hline(yintercept = 500,colour="red") +labs(x = "TCW_14311_run1")
vp3<-vps[[3]]+geom_hline(yintercept = 14,colour="red") +labs(x = "TCW_14311_run1")

pdf(paste(dir, "/QC_metric_VP_before_cutoffs.pdf", sep = ""), width = 10, height = 10)
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
plot(p)
dev.off()

#QC metric violin plots with ribo
vps<-VlnPlot(object = T_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, slot = "counts", pt.size = 0)
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="red") +labs(x = "TCW_14311_run1")
vp2<-vps[[2]]+geom_hline(yintercept = 75000,colour="red") +geom_hline(yintercept = 500,colour="red") +labs(x = "TCW_14311_run1")
vp3<-vps[[3]]+geom_hline(yintercept = 14,colour="red") +labs(x = "TCW_14311_run1")
vp4<-vps[[4]]+geom_hline(yintercept = 30,colour="red") +labs(x = "TCW_14311_run1")

pdf(paste(dir, "/QC_metric_VP_before_cutoffs_with_ribo.pdf", sep = ""), width = 10, height = 10)
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank()) +vp4 + theme(axis.text.x = element_blank())
plot(p)
dev.off()

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(organoid2, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "organoid") & NoLegend() 
plot1
plot2 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/QC_metrics_FP_before_cutoffs.pdf", width = 10, height = 5)
plot1 + plot2
dev.off()

#with ribo
plot1 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.ribo") +labs(legend = "TCW_14311_run1") & NoLegend() 
plot3 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/QC_metrics_with_ribo_FP_before_cutoffs.pdf", width = 10, height = 5)
plot1 + plot2 + plot3
dev.off()

#DONE

#QC on QC Metrics

#Remove Cells with High MT Percentage
#MiQC FILTERING: REMOVE low quality cells and set MT threshold
#vingette: https://bioconductor.org/packages/devel/bioc/vignettes/miQC/inst/doc/miQC.html
#BiocManager::install("flexmix")
library(miQC)
library(SeuratWrappers)
library(flexmix)

#Run MiQC Algorithm
#note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
#Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.
T_QC <- RunMiQC(T_QC, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.97, model.slot = "flexmix_model") 

#Look at parameters and posterior values
#flexmix::parameters(Misc(T_QC, "flexmix_model"))
#head(flexmix::posterior(Misc(T_QC, "flexmix_model")))

pdf(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/stringent_miQC_probability_plot.pdf", width = 15, height = 10)
PlotMiQC(T_QC, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
dev.off()

pdf(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/stringent_miQC_keep_plot.pdf", width = 15, height = 10)
PlotMiQC(T_QC, color.by = "miQC.keep")
dev.off()

#Run MiQC Filtering
T_QC <- subset(T_QC, miQC.keep == "keep") 
save(T_QC, file = "T_QC")
#7,305 cells and 23,412 genes

##* Do this for the combined Seurat object *##
#Check new mito.pct vs. nCount_RNA scatter
plot1 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.mt")
pdf(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/QC_mito_FP_after_mt_cutoff.pdf", width = 15, height = 10)
plot1 
dev.off()

#Remove Additional Outliers
#Only include cells with at least 200 genes expressed
T_QC <- subset(T_QC, subset = nFeature_RNA > 200) #number of genes detected in each cell
T_QC
#7,305 cells and 23,412 genes

#Only include cells with at least 500 molecules expressed
T_QC <- subset(T_QC, subset = nCount_RNA > 500) #number of genes detected in each cell
T_QC

#7,305 cells and 23,412 genes

#Remove the upper tail observed on the "nCount_RNA" violin plot (99.5% percentile)
ub <- quantile(T_QC[["nCount_RNA"]]$nCount_RNA, probs = 0.99)
T_QC <- T_QC[, T_QC[["nCount_RNA"]] < ub]
#7,231 cells and 23,412 genes

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/QC_metrics_FP_after_cutoffs.pdf", width = 10, height = 5)
plot1 + plot2 
dev.off()

#with ribo
plot1 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.ribo") +labs(legend = "TCW_14311_run1") & NoLegend() 
plot3 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 
pdf(file = "/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/QC_metrics_with_ribo_FP_after_cutoffs.pdf", width = 10, height = 5)
plot1 + plot2 + plot3
dev.off()

#QC metric violin plots
vps<-VlnPlot(object = T_QC_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, slot = "counts", pt.size = 0) 
vp1<-vps[[1]] +labs(x = "TCW_14311_run1")
vp2<-vps[[2]]+geom_hline(yintercept = 75000,colour="red") +labs(x = "TCW_14311_run1")
vp3<-vps[[3]]+geom_hline(yintercept = 14,colour="red") +labs(x = "TCW_14311_run1")

pdf("/projectnb/tcwlab/LabMember/dmurthy/projects/Single_Cell_RNA_Seq/TCW_SC_RNA_Seq_Analysis/Files/Redone_QC/QC_metric_VP_after_cutoffs.pdf", width = 10, height = 10)
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank())
plot(p)
dev.off()

#QC metric violin plots with ribo
vps<-VlnPlot(object = T_QC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, slot = "counts", pt.size = 0)
vp1<-vps[[1]] +geom_hline(yintercept = 200,colour="red") +labs(x = "TCW_14311_run1")
vp2<-vps[[2]]+geom_hline(yintercept = 75000,colour="red") +geom_hline(yintercept = 500,colour="red") +labs(x = "TCW_14311_run1")
vp3<-vps[[3]]+geom_hline(yintercept = 14,colour="red") +labs(x = "TCW_14311_run1")
vp4<-vps[[4]]+geom_hline(yintercept = 30,colour="red") +labs(x = "TCW_14311_run1")

pdf(paste(dir, "/QC_metric_VP_after_cutoffs_with_ribo.pdf", sep = ""), width = 10, height = 10)
p <- vp1 +  theme(axis.text.x = element_blank()) +vp2 + theme(axis.text.x = element_blank()) +vp3 + theme(axis.text.x = element_blank()) +vp4 + theme(axis.text.x = element_blank())
plot(p)
dev.off()

save(T_QC, file = "T_QC")

max(T_QC$nFeature_RNA) #9016
max(T_QC$nCount_RNA) #72,529

#Number of cells per sample (hash.ID)
pdf(paste(dir, "/number_cells_per_sample.pdf", sep = ""), width = 10, height = 10)
sample_mat <- as.data.frame(T_QC$hash.ID) %>%
  rename("hash.ID" = "T_QC$hash.ID")
ggplot(sample_mat, aes(hash.ID)) + geom_bar() + theme_bw()+ theme(axis.text.x = element_text(angle = 90))
dev.off()

#Number of cells per individual
pdf(paste(dir, "/number_cells_per_individual.pdf", sep = ""), width = 10, height = 10)
sample_mat <- as.data.frame(T_QC$individual) %>%
  rename("hash.ID" = "T_QC$individual")
ggplot(sample_mat, aes(hash.ID)) + geom_bar() + theme_bw()+ theme(axis.text.x = element_text(angle = 90))
dev.off()

#Number of cells per APOE genotype
pdf(paste(dir, "/number_cells_per_APOE_genotype.pdf", sep = ""), width = 10, height = 10)
sample_mat <- as.data.frame(T_QC$APOE_Genotype) %>%
  rename("hash.ID" = "T_QC$APOE_Genotype")
ggplot(sample_mat, aes(hash.ID)) + geom_bar() + theme_bw()+ theme(axis.text.x = element_text(angle = 90))
dev.off()

################################################################################
#Scaling & Normalization

#RNA assay: raw count data (will be scaled & normalized) #use for visualization & some DEA analysis
#SCT assay: data = log of RNA counts, scale.data = normalized counts (input to PCA), counts = log normalized corrected UMI counts. #SCT transform replaces NormalizeData(), ScaleData(), and FindVariableFeatures().

#RNA assay: Scale Data & Normalize
#scale data while removing variation due to mt percentage & cell phase
all.genes <- rownames(T_QC)
T_QC <- ScaleData(T_QC, features = all.genes, assay = "RNA", vars.to.regress = c("percent.mt", "CC.Difference"), verbose = FALSE)
T_QC <- NormalizeData(T_QC, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
T_QC <- FindVariableFeatures(T_QC, selection.method = "vst", nfeatures = 2000)

#SCT assay
#scale data while removing variation due to mt percentage & cell phase
T_QC <-SCTransform(T_QC, vars.to.regress = c("percent.mt", "CC.Difference"))


#Cell cycle effect: cell cycle score is computed after normalization
T_QC <- RunPCA(T_QC, features = VariableFeatures(T_QC))
T_QC <- RunUMAP(T_QC, dims = 1:50)
T_QC <- CellCycleScoring(T_QC, s.features = cc.genes.updated.2019$s.genes,
                         g2m.features = cc.genes.updated.2019$g2m.genes,
                         set.ident = TRUE,
                         search=TRUE)

T_QC[["CC.Difference"]]<- T_QC$G2M.Score - T_QC$S.Score

# Visualize the PCA, grouping by cell cycle phase
DimPlot(organoid1,
        reduction = "pca",
        group.by= "genotype")

organoid2@meta.data$genotype
# Visualize PCA, grouping by Individual 
DimPlot(T_QC,
        reduction = "pca",
        group.by= "individual")

# Visualize PCA, grouping by APOE Genotype 
DimPlot(T_QC,
        reduction = "pca",
        group.by= "APOE Genotype")

# Visualize PCA, grouping by APOE Genotype 
DimPlot(T_QC,
        reduction = "pca",
        group.by= "A")


#PCA
T_QC <- RunPCA(T_QC, features = VariableFeatures(object = T_QC), assay = "SCT")
DimPlot(T_QC, reduction = "pca")

#PCA Plots
VizDimLoadings(T_QC, dims = 1:2, reduction = "pca")
dev.off()

plot1 <- VariableFeaturePlot(T_QC.MiQC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/variablefeatureplot_filtered_object.pdf", width = 15, height = 10)
plot1+plot2
dev.off()

#UMAP & Clustering
T_QC <- RunUMAP(T_QC, dims = 1:50)
T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, k.param = 30, verbose = FALSE)
T_QC.MiQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected

#umap
T_QC.MiQC <- RunUMAP(T_QC.MiQC, assay = SCT, dims = 1:50,  n.neighbors = 30, verbose = FALSE)

#Label 


#Harmony for Individual Effect


#TCW test
TCW <- CellCycleScoring(TCW, s.features = cc.genes.updated.2019$s.genes,
                        g2m.features = cc.genes.updated.2019$g2m.genes,
                        set.ident = TRUE,
                        search=TRUE)

TCW[["CC.Difference"]]<- TCW$G2M.Score - TCW$S.Score #correct for G2M vs S

TCW$CC.Difference

RidgePlot(TCW, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

TCW <- RunPCA(TCW, features = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes))

DimPlot(TCW, reduction = "pca")

# Plot the PCA colored by cell cycle phase
DimPlot(TCW,
        reduction = "pca",
        group.by= "pct.mito")

DimPlot(TCW, reduction = "harmony_40_nn_0.4res_umap", group.by = "Phase")

DimPlot(TCW, reduction = "no_harmony_50_nn_0.4_res_umap", group.by = "Phase")

#nCount_RNA & nFeature_RNA
VlnPlot(TCW, features = c("nCount_RNA"), pt.size = 0.5)
VlnPlot(TCW, features = c("nFeature_RNA"), pt.size = 0.5)

plot1 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.mt") + labs(legend = "TCW_14311_run1") & NoLegend() 
plot2 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +labs(legend = "TCW_14311_run1") & NoLegend() 

plot1 + plot2

####################################################################################################################
####################################################################################################################
####################################################################################################################
#OLDER QC
####################################################################################################################
#Original Clustering UMAP Visualization (prior to removal of high mito gene expressing cells)
#find variable features

#normalize data: SCT
T_QC <- SCTransform(T_QC, verbose = TRUE)

T_QC <- FindVariableFeatures(T_QC, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(T_QC), 10)

plot1 <- VariableFeaturePlot(T_QC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/variablefeatureplot_original_object.pdf", width = 15, height = 10)
plot1+plot2
dev.off()

#scale data
all.genes <- rownames(T_QC)
T_QC <- ScaleData(T_QC, features = all.genes)

#pca
T_QC <- RunPCA(T_QC, features = VariableFeatures(object = T_QC), assay = "SCT")

#examine & visualize pca results
print(T_QC[["pca"]], dims = 1:5, nfeatures = 5)
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pc1&pc2_original_object.pdf", width = 15, height = 10)
VizDimLoadings(T_QC, dims = 1:2, reduction = "pca")
dev.off()
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pc1&pc2dimplot_original_object.pdf", width = 15, height = 10)
DimPlot(T_QC, reduction = "pca")
dev.off()

#pca heatmap #downsampled to 500 cells
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pc1_heatmap_original_object.pdf", width = 15, height = 10)
DimHeatmap(T_QC, dims = 1, cells = 500, balanced = TRUE)
dev.off()
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pcs1-15_heatmap_original_object.pdf", width = 30, height = 30)
DimHeatmap(T_QC, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pca_elbow_plot_original_object.pdf", width = 30, height = 30)
ElbowPlot(T_QC)
dev.off()

#cluster the cells based on previously identified PCs
#KNN: K-nearest neighbor network: Graph-Theoretic Approach
T_QC <- FindNeighbors(T_QC, dims = 1:50, verbose = FALSE)
T_QC <- FindClusters(T_QC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected
head(Idents(T_QC), 5)
#24 clusters

#umap
T_QC <- RunUMAP(T_QC, assay = SCT, dims = 1:50, verbose = FALSE)

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/umap_bycluster_original_object.pdf", width = 10, height = 10)
DimPlot(T_QC, reduction = "umap", group.by = "seurat_clusters", pt.size = .1) +ggtitle("UMAP by Seurat Cluster") + annotate("text", x=8.5, y=-20, label= "Cells = 7,911")
dev.off()

####################################################################################################################
#MiQC FILTERING: REMOVE low quality cells and set MT threshold
#vingette: https://bioconductor.org/packages/devel/bioc/vignettes/miQC/inst/doc/miQC.html
library(miQC)
library(SeuratWrappers)
library(flexmix)

#label MT percent
T_QC[["percent.mt"]] <- PercentageFeatureSet(T_QC, pattern = "^MT-")

# Visualize QC metrics as a violin plot
Idents(T_QC)
Idents(T_QC) <- T_QC@project.name
Idents(T_QC) <- T_QC$seurat_clusters

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qc_violin_plots_original_object.pdf", width = 15, height = 10)
#pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qc_violin_plots_bycluster_original_object.pdf", width = 15, height = 10)
VlnPlot(organoid2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), slot = "counts", pt.size = 0, split.by = NULL, ncol = 3) 
dev.off()

#FeatureScatter is used to visualize feature-feature relationships
#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T_QC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qcfeatureplots_original_object.pdf", width = 15, height = 10)
plot1 + plot2
dev.off()

#Run MiQC Algorithm
#note: posterior cutoff = the posterior probability of a cell being part of the compromised distribution, a number between 0 and 1.
#Any cells below the appointed cutoff will be marked to keep. Defaults to 0.75.
T_QC.MiQC <- RunMiQC(T_QC, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA", posterior.cutoff = 0.97, 
                     model.slot = "flexmix_model") #75, 97

#Look at parameters and posterior values
flexmix::parameters(Misc(T_QC.MiQC, "flexmix_model"))
head(flexmix::posterior(Misc(T_QC.MiQC, "flexmix_model")))

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/stringent_miQC_probability_plot.pdf", width = 15, height = 10)
PlotMiQC(T_QC.MiQC, color.by = "miQC.probability") + ggplot2::scale_color_gradient(low = "grey", high = "purple")
dev.off()

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/stringent_miQC_keep_plot.pdf", width = 15, height = 10)
PlotMiQC(T_QC.MiQC, color.by = "miQC.keep")
dev.off()

#Run MiQC Filtering
T_QC.MiQC <- subset(T_QC.MiQC, miQC.keep == "keep") 
#7,337 cells and 22,719 genes

#Check QC-plots of filtered data after MiQC
Idents(T_QC.MiQC) <- T_QC.MiQC@project.name
Idents(T_QC.MiQC) <- T_QC.MiQC$seurat_clusters
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qcplots_bycluster_filtered_object.pdf", width = 15, height = 10)
#pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qcplots_filtered_object.pdf", width = 15, height = 10)
VlnPlot(T_QC.MiQC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), slot = "counts", pt.size = 0, split.by = NULL, ncol = 3)
dev.off()

#Ensure linear relationship between nFeature_RNA and nCount_RNA
plot1 <- FeatureScatter(T_QC.MiQC, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(T_QC.MiQC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qcfeatureplots_filtered_object.pdf", width = 15, height = 10)
plot1 + plot2
dev.off()

#Remove Additional Outliers (added 1/5/23)
#Only include cells with at least 200 genes expressed
T_QC.MiQC <- subset(T_QC.MiQC, subset = nFeature_RNA > 200) 
T_QC.MiQC
#7,337 cells and 22,719 genes

#Remove the upper tail observed on the "nCount_RNA" violin plot (99.5% percentile)
ub <- quantile(T_QC.MiQC[["nCount_RNA"]]$nCount_RNA, probs = 0.995)
T_QC.MiQC <- T_QC.MiQC[, T_QC.MiQC[["nCount_RNA"]] < ub]
T_QC.MiQC
#7,300 cells and 22,547 genes

#Check violin plots after outlier removal
Idents(object = T_QC.MiQC) <- "project"
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qc_violin_plots_final_filtered_object.pdf", width = 15, height = 10)
VlnPlot(T_QC.MiQC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), slot = "counts", pt.size = 0, split.by = NULL, ncol = 3) 
dev.off()

####################################################################################################################
#Re-cluster dataset after filtering

#normalize data: SCT
T_QC.MiQC <- SCTransform(T_QC.MiQC, verbose = TRUE)

#find variable features: find the 2000 most variable genes
T_QC.MiQC <- FindVariableFeatures(T_QC.MiQC, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
top10 <- head(VariableFeatures(T_QC.MiQC), 10)
top10

plot1 <- VariableFeaturePlot(T_QC.MiQC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf(file = "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/variablefeatureplot_filtered_object.pdf", width = 15, height = 10)
plot1+plot2
dev.off()

#scale data while removing variation due to mt percentage
all.genes <- rownames(T_QC.MiQC)
T_QC.MiQC <- ScaleData(T_QC.MiQC, features = all.genes, vars.to.regress = "percent.mt")

#pca
T_QC.MiQC <- RunPCA(T_QC.MiQC, features = VariableFeatures(object = T_QC.MiQC), assay = "SCT")

#examine & visualize pca results
print(T_QC[["pca"]], dims = 1:5, nfeatures = 5)
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pc1&pc2_filtered_object.pdf", width = 15, height = 20)
VizDimLoadings(T_QC, dims = 1:2, reduction = "pca")
dev.off()
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pc1&pc2dimplot_filtered_object.pdf", width = 15, height = 20)
DimPlot(T_QC.MiQC, reduction = "pca")
dev.off()

###################################
#3D pca plot
#date: 1/7/23
#source: https://rpubs.com/HWH/920093
library(plotly)

# xaxis
axx <- list(
  nticks = 4,
  range = c(-10,10) #select range of xaxis
)

# yaxis
axy <- list(
  nticks = 4,
  range = c(-10,10) #select range of yaxis
)

#zaxis
axz <- list(
  nticks = 4,
  range = c(-10,10) #select range of zaxis
)


# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = T_QC.MiQC, vars = c("PC_1", "PC_2", "PC_3", "seurat_clusters"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

seurat_default_colors <- c("#F8766D", "#E68613", "#CD9600", "#ABA300", "#7CAE00", "#0CB702", "#00BE67",
                           "#00C19A", "#00BFC4", "#00B8E7", "#00A9FF", "#8494FF", "#C77CFF", "#ED68ED", "#FF61CC", "#FF68A1")

fig <- plot_ly(data = plot.data, 
               x = ~PC_1, y = ~PC_2, z = ~PC_3, 
               color = ~seurat_clusters, 
               colors = seurat_default_colors,
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 5, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect
fig_cube

#########################################
#pca heatmap #downsampled to 500 cells
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pc1_heatmap_filtered_object.pdf")
DimHeatmap(T_QC.MiQC, dims = 1, cells = 500, balanced = TRUE) 
#+ggplot2::scale_fill_gradientn(colors = c("steelblue3", "white", "tomato2"))
dev.off()
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pcs1-15_heatmap_filtered_object.pdf")
DimHeatmap(T_QC.MiQC, dims = 1:15, cells = 500, balanced = TRUE)
#+ggplot2::scale_fill_gradientn(colors = c("steelblue1", "white", "tomato"))
dev.off()

#pca: elbow plot
pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/pca_elbow_plot.pdf", width = 30, height = 30)
ElbowPlot(T_QC.MiQC)
dev.off()

#cluster the cells based on previously identified PCs
#KNN: K-nearest neighbor network: Graph-Theoretic Approach
T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, verbose = FALSE)
T_QC.MiQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected
head(Idents(T_QC.MiQC), 5)
#22 clusters

#umap
T_QC.MiQC <- RunUMAP(T_QC.MiQC, assay = SCT, dims = 1:50, verbose = FALSE)

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/umap_bycluster_filtered_object.pdf", width = 10, height = 10)
DimPlot(T_QC.MiQC, reduction = "umap", group.by = "seurat_clusters", pt.size = .1) +ggtitle("UMAP by Seurat Cluster") + annotate("text", x=8.5, y=-20, label= "Cells = 7,300")
dev.off()

save(T_QC.MiQC, file = "T_QC.MiQC")
load(file = "T_QC.MiQC")

T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, verbose = FALSE)
T_QC.MIQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected
head(Idents(T_QC.MC), 5)
#res = 0.2: 16 clusters
T_QC.MIQC <- RunUMAP(T_QC.MIQC, assay = SCT, dims = 1:50, verbose = FALSE)
DimPlot(T_QC.MIQC, reduction = "umap", group.by = "seurat_clusters", pt.size = .1)

####################################################################################################################
#1/1/23 - 1/5/23: Find Clusters Resolution Alteration: Decreasing the resolution to obtain
#more general groupings for more representative cell type labels
#note: a resolution closer to 1 gives more clusters and closer to 0 gives less clusters

#Label previous clusters obtained using res = 0.6

T_QC.MiQC$seurat_clusters_0.6 <- T_QC.MiQC$seurat_clusters

#cluster the cells based on previously identified PCs
#KNN: K-nearest neighbor network: Graph-Theoretic Approach
T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, verbose = FALSE)
T_QC.MiQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected
head(Idents(T_QC.MiQC), 5)
#res = 0.2: 14 clusters
#res = 0.1: 10 clusters

#umap
T_QC.MiQC <- RunUMAP(T_QC.MiQC, assay = SCT, dims = 1:50, verbose = FALSE)

#Label previous clusters obtained using res = 0.1
T_QC.MiQC$seurat_clusters_0.1  <- T_QC.MiQC$seurat_clusters

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/umap_bycluster_filtered_object_res_0.1.pdf", width = 10, height = 10)
DimPlot(T_QC.MiQC, reduction = "umap", group.by = "seurat_clusters_0.1", pt.size = .1) +ggtitle("UMAP by Seurat Cluster") + annotate("text", x=8.5, y=-20, label= "Cells = 7,300")
dev.off()

#Label previous clusters obtained using res = 0.2
T_QC.MiQC$seurat_clusters_0.2  <- T_QC.MiQC$seurat_clusters

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/umap_bycluster_filtered_object_res_0.2.pdf", width = 10, height = 10)
DimPlot(T_QC.MiQC, reduction = "umap", group.by = "seurat_clusters_0.2", pt.size = .1) +ggtitle("UMAP by Seurat Cluster") + annotate("text", x=8.5, y=-20, label= "Cells = 7,300")
dev.off()

save(T_QC.MiQC, file = "T_QC.MiQC")

###############################
#note: more neighbors preserves more global structure giving less detail to local
#in general this parameter should often be in the range 5 to 50.

T_QC.MiQC <- FindNeighbors(T_QC.MiQC, dims = 1:50, k.param = 30, verbose = FALSE)
T_QC.MiQC <- FindClusters(T_QC.MiQC, resolution = 0.6, verbose = FALSE) #change resolution to experiment, alters # of clusters detected

#umap
T_QC.MiQC <- RunUMAP(T_QC.MiQC, assay = SCT, dims = 1:50,  n.neighbors = 30, verbose = FALSE)

#Label 
T_QC.MiQC$umap_35_neighbors 

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/umap_bycluster_filtered_object_30_neighbors.pdf", width = 10, height = 10)
DimPlot(T_QC.MiQC, reduction = "umap", group.by = "seurat_clusters_0.6", pt.size = .1) +ggtitle("UMAP by Seurat Cluster") + annotate("text", x=8.5, y=-20, label= "Cells = 7,300")
dev.off()



########################################################################
#Check QC Metrics Across Clusters

metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")

pdf(file = "/Users/deeptimurthy/Desktop/TCW_Lab_Work/TCW_2020_scDATA_Analysis/Final_Project_Outputs/2023_Analysis/qc_metrics_featureplots_filtered_object_no_min_cutoff.pdf", width = 10, height = 10)

FeaturePlot(T_QC.MiQC, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            sort.cell = TRUE,
            label = TRUE)

dev.off()

