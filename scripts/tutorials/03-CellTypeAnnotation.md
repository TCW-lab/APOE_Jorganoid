## Intro

We begin by loading in R libraries as well as setting our working directory.

```{r}
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
library(scRNAseq)
library(SingleR)
library(tidyverse)
dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation"
```

## Step 1: Read in Reference Data

Read in reference data to compare your data to. In this case, we are using Polioudakis et al. (2019) as our reference. They made their data publicly available through their github.

```{r}
#Reference Atlas: Polioudakis et. al.

ref_dir <- "/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/ref_data"
organoid <- readRDS('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/outputs/03-CellTypeAnnotation/organoid_post_harmony.rds')
#load count matrix
load(paste(ref_dir, "/sc_dev_cortex_geschwind/raw_counts_mat.rdata", sep = ""))

raw.counts.mat <- as.matrix(raw_counts_mat)
reffy <- CreateSeuratObject(raw.counts.mat)
```

## Step 2: Load Metadata and Process Reference Data

In this step we are appending metadata to the reference object and then convert it to a Single Cell Experiment object. This is necessary in performing SingleR analysis in the next step.

```{r}
#load metadata
ref_metadata <- read.csv(file = paste(ref_dir, "/sc_dev_cortex_geschwind/cell_metadata.csv", sep = ""))
rownames(ref_metadata) <- ref_metadata[,1]
ref_metadata[,1] <- NULL
ref_metadata[,1]
#append metadata to reference object
reffy <- AddMetaData(reffy, ref_metadata)
reffy <- reffy[,!is.na(reffy$Cluster)]

#convert to a single cell experiment
reffo <- as.SingleCellExperiment(reffy)
metadata(reffo)<- ref_metadata

#normalize data
library(scuttle)
reffo <- logNormCounts(reffo)

```

## Step 3a: Run SingleR on the 

Here we can refer to our organoid object and use SingleR to compare that to the reference data we just processed. Note the keep_rows <- !grepl(...) step here. This is to ensure that the ENSEMBL IDs that we kept in the previous tutorial are now omitted, as it is necessary to have all of our rownames in one format. The reference data uses gene symbols, so, because most of our gene names are in gene symbol format, we can just exclude the ENSEMBL ID names that were unsuccessfully converted over.

```{r}
sce <- as.SingleCellExperiment(organoid)
keep_rows <- !grepl("^ENSG", rownames(sce))
# keep only gene symbols.
sce <- sce[keep_rows, ]
sce_counts <- counts(sce)
reffo_counts <- counts(reffo)

# intersect names of sce / reffo
common_row_names <- intersect(rownames(sce), rownames(reffo))

# subset sce to only be rownames in common w reffo
sce_filtered <- sce[common_row_names, ]

# Check to ensure 'sce_filtered' contains only the common rows
rownames(sce_filtered) 

sce_counts_dense <- as.matrix(sce_counts)
reffo_counts_dense <- as.matrix(reffo_counts) 
# Run SingleR
singleR_results <- SingleR(sc_data = sce_counts_dense,
                           ref_data = reffo_counts_dense, 
                           types = factor(reffo@colData$Cluster),
                           numCores = 14)
filename = 'organoid_SingleR.rds'
saveRDS(singleR_results, file.path(dir, filename))
```


## Step 3b: Run SingleR on a Split Seurat Object

Alternatively, if our organoid object is too large, we can split it by sample ID (defined in the meta.data in the previous tutorial script, 02-scRNAseq-QC) and run it as a QSub job. Please refer to the scripts folder of the Github repository to see how to accomplish this. The general framework of it is here (note that the samples vector is sliced in order to run it in separate QSub batch jobs - in my case I have 24 so I am running it across 3-4 batch jobs, with 6-8 samples each time):

```{r}
samples <- unique(organoid@meta.data$sample)

for (sample_id in samples[18:24]) {
  sce <- as.SingleCellExperiment(subset(organoid, subset = sample == sample_id))
  keep_rows <- !grepl("^ENSG", rownames(sce))
  # keep only gene symbols.
  sce <- sce[keep_rows, ]
  sce_counts <- counts(sce)
  reffo_counts <- counts(reffo)
  
  # intersect names of sce / reffo
  common_row_names <- intersect(rownames(sce), rownames(reffo))
  
  # subset sce to only be rownames in common w reffo
  sce_filtered <- sce[common_row_names, ]
  
  # Check to ensure 'sce_filtered' contains only the common rows
  rownames(sce_filtered) 
  
  sce_counts_dense <- as.matrix(sce_counts)
  reffo_counts_dense <- as.matrix(reffo_counts)
  # Run SingleR
  singleR_results <- SingleR(sc_data = sce_counts_dense,
                             ref_data = reffo_counts_dense, types = factor(reffo@colData$Cluster),
                             numCores = 14)
  filename = paste0('organoid_SingleR_',sample_id,'.rds')
  saveRDS(singleR_results, file.path(dir, filename))
}

```

