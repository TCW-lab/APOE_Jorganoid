setwd('/projectnb/tcwlab/LabMember/akg/APOE_Jorganoid_project/')

out<-'outputs/01-get_seurat_object'
dir.create(out)

library(tools)
library(Seurat)
library(data.table)
library(ggplot2)

mat<-ReadMtx('outputs/15_S5/af_quant/alevin/quants_mat.mtx',
        cells ='outputs/15_S5/af_quant/alevin/quants_mat_cols.txt' ,
        features = 'outputs/15_S5/af_quant/alevin/quants_mat_rows.txt',feature.column=1
          )
head(rownames(mat))
nrow(mat) #191 363

#get the number of umis per cell
cells_count<-data.table(cell=rownames(mat),
                        count=rowSums(mat))
cells_count[,cell_rank:=rank(-count)]
#plot the 'knee' plot
ggplot(cells_count,aes(x=cell_rank,y=count))+geom_line()+
  scale_x_log10()+scale_y_log10()+theme_bw()+
  geom_hline(yintercept = 1000) + ggtitle('S5')
#here test different threshold to find the best one that is between the 2 knees, or at the middle the first knee
#we expecrt cutoff between 100 and 2000 UMIS (can be more if the sequencing depth or cell RNA amount is bigger)
#we expect to have between 2000 and 25k cells
cell_cutoff=1000
sum(rowSums(mat)>cell_cutoff) #number of cells 

cells<-rownames(mat)[rowSums(mat)>cell_cutoff] 
#cellS5<-cells_count[count>cell_cutoff]$cell #same
matf<-mat[cells,]

S5<-CreateSeuratObject(t(as.matrix(matf)),project = 'Jorganoid')
# An object of class Seurat 
# 174657 features across 5460 samples within 1 assay 
# Active assay: RNA (174657 features, 0 variable features)
# 1 layer present: counts
S5
saveRDS(S5,file.path(out,'S5.rds'))

S5$sample<-'S5'

#did for all


#Merge the Seurat objects
S1 <- readRDS(file = "./outputs/01-get_seurat_object/S1.rds")
S2 <- readRDS(file = "./outputs/01-get_seurat_object/S2.rds")
S3 <- readRDS(file = "./outputs/01-get_seurat_object/S3.rds")
S4 <- readRDS(file = "./outputs/01-get_seurat_object/S4.rds")
S5 <- readRDS(file = "./outputs/01-get_seurat_object/S5.rds")
S6 <- readRDS(file = "./outputs/01-get_seurat_object/S6.rds")
S7 <- readRDS(file = "./outputs/01-get_seurat_object/S7.rds")
S8 <- readRDS(file = "./outputs/01-get_seurat_object/S8.rds")
S9 <- readRDS(file = "./outputs/01-get_seurat_object/S9.rds")
S10 <- readRDS(file = "./outputs/01-get_seurat_object/S10.rds")
S11 <- readRDS(file = "./outputs/01-get_seurat_object/S11.rds")
S12 <- readRDS(file = "./outputs/01-get_seurat_object/S12.rds")
S13 <- readRDS(file = "./outputs/01-get_seurat_object/S13.rds")
S14 <- readRDS(file = "./outputs/01-get_seurat_object/S14.rds")
S15 <- readRDS(file = "./outputs/01-get_seurat_object/S15.rds")
S16 <- readRDS(file = "./outputs/01-get_seurat_object/S16.rds")
S17 <- readRDS(file = "./outputs/01-get_seurat_object/S17.rds")
S18 <- readRDS(file = "./outputs/01-get_seurat_object/S18.rds")
S19 <- readRDS(file = "./outputs/01-get_seurat_object/S19.rds")
S20 <- readRDS(file = "./outputs/01-get_seurat_object/S20.rds")
S21 <- readRDS(file = "./outputs/01-get_seurat_object/S21.rds")
S22 <- readRDS(file = "./outputs/01-get_seurat_object/S22.rds")
S23 <- readRDS(file = "./outputs/01-get_seurat_object/S23.rds")
S24 <- readRDS(file = "./outputs/01-get_seurat_object/S24.rds")

# install.packages('SeuratData')
library(remotes)
# remotes::install_github("hoxo-m/githubinstall")
library(githubinstall)

# githubinstall("SatijaLab/SeuratData")
library(SeuratData)


# unloadNamespace("fastmap")

# detach("package:fastmap", unload = TRUE)
# install.packages("fastmap")
# devtools::install_github('satijalab/seurat-data')



APOE_Jorganoid <- merge(S1, y = c(S2, S3, S4, S5, S6, S7, S8, S9, S10,
                            S11, S12, S13, S14, S15, S16, S17, S18, S19,
                            S20, S21, S22, S23, S24), 
                  add.cell.ids = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                                  "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19",
                                  "S20", "S21", "S22", "S23", "S24"),
                  project = "APOE_Jorganoid"
                  )
APOE_Jorganoid <- NormalizeData(APOE_Jorganoid)
APOE_Jorganoid <- FindVariableFeatures(APOE_Jorganoid)
APOE_Jorganoid <- ScaleData(APOE_Jorganoid)

APOE_Jorganoid <- RunPCA(APOE_Jorganoid, features = VariableFeatures(object = APOE_Jorganoid))
DimPlot(APOE_Jorganoid, reduction = "pca")
ElbowPlot(object = APOE_Jorganoid)


APOE_Jorganoid <- RunUMAP(APOE_Jorganoid, dims = 1:10)
DimPlot(APOE_Jorganoid, reduction = "umap")

saveRDS(APOE_Jorganoid,file.path(out,'APOE_Jorganoid.rds'))

#---- end of data pre-processing ----#
# APOE_Jorganoid <- readRDS(file = "./outputs/01-get_seurat_object/APOE_Jorganoid.rds")

ggplot(APOE_Jorganoid@reductions$umap@cell.embeddings, aes(x=UMAP_1, y=UMAP_2, color=orig.ident)) + geom_point()

#--- Data Analysis for CellRanger pipeline ---#



#--- end of CellRanger ---#
#--- Beginning of downstream data analysis ---#


ggplot()
