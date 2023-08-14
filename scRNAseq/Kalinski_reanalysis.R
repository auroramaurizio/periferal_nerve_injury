suppressMessages(library(ggbeeswarm))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
library(scales)
library(viridis)
library(SeuratWrappers)
library(slingshot)
require(BiocStyle)
library(data.table)
library(SingleCellExperiment)
#print(version) R version 3.6.1
#Bonanomi
#packageVersion("Seurat")
sessionInfo()


###############################################
# Examine Uninjured post injury samples from kalinski et al 2020.
# Analysis of the immune response to sciatic nerve injury identifies efferocytosis as a key mechanism of nerve debridement
###############################################



rep1_path='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Kalinski/GSE153762_RAW/rep1'

list.files(rep1_path)
rep1 <- Read10X(data.dir = rep1_path)
rep1 <- CreateSeuratObject(counts = rep1, project = "rep1", min.cells = 5, min.features = 200)
rep1[["percent.mt"]] <- PercentageFeatureSet(rep1, pattern = "^mt-")
rep1[["percent.Rpl"]] <- PercentageFeatureSet(rep1, pattern = "Rpl")
rep1$stim <- "3D_kalinski_rep1"
options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(rep1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 7500
max_percent_MT = 20



rep2_path='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Kalinski/GSE153762_RAW/rep2'

list.files(rep2_path)
rep2 <- Read10X(data.dir = rep2_path)
rep2 <- CreateSeuratObject(counts = rep2, project = "rep2", min.cells = 5, min.features = 200)
rep2[["percent.mt"]] <- PercentageFeatureSet(rep2, pattern = "^mt-")
rep2[["percent.Rpl"]] <- PercentageFeatureSet(rep2, pattern = "Rpl")
rep2$stim <- "3D_kalinski_rep2"
options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(rep2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 7500
max_percent_MT = 20



rep3_path='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Kalinski/GSE153762_RAW/rep3'

list.files(rep3_path)
rep3 <- Read10X(data.dir = rep3_path)
rep3 <- CreateSeuratObject(counts = rep3, project = "rep3", min.cells = 5, min.features = 200)
rep3[["percent.mt"]] <- PercentageFeatureSet(rep3, pattern = "^mt-")
rep3[["percent.Rpl"]] <- PercentageFeatureSet(rep3, pattern = "Rpl")
rep3$stim <- "3D_kalinski_rep3"
options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(rep3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 7500
max_percent_MT = 20

# normalize and find varible features in the objects
object_clean_new.list <- lapply(X = c(rep1, rep2, rep3 ), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# find integration anchors
integrated <- FindIntegrationAnchors(object.list = object_clean_new.list, dims = 1:35)
features.to.integrate = integrated@anchor.features
integrated <- IntegrateData(anchorset = integrated, dims = 1:35, features.to.integrate = features.to.integrate)
DefaultAssay(integrated) <- "integrated"
table(integrated$stim)

res = 1

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = res)
DimPlot(integrated, reduction = "umap", group.by = "stim")
DimPlot(integrated, reduction = "umap", split.by = "stim")

DefaultAssay(integrated) = "RNA"

saveRDS(integrated, "integrated_all_CT_Kalinski.Rds")
saveRDS(integrated, "integrated_all_CT_Kalinski_bis.Rds")

FeaturePlot(integrated, "Rgs5")

integrated_sub = subset(integrated, idents = c("20","10","32","25"))
integrated_sub <- NormalizeData(integrated_sub)
integrated_sub<- FindVariableFeatures(integrated_sub, selection.method = "vst", nfeatures = 2000)
res = 1
# Run the standard workflow for visualization and clustering
integrated_sub <- ScaleData(integrated_sub, verbose = FALSE)
integrated_sub <- RunPCA(integrated_sub, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated_sub <- RunUMAP(integrated_sub, reduction = "pca", dims = 1:30)
integrated_sub <- FindNeighbors(integrated_sub, reduction = "pca", dims = 1:30)
integrated_sub <- FindClusters(integrated_sub, resolution = res)

#identify EC after removing pericytes

integrated_sub_new = subset(integrated_sub, idents = c("0","1","2","3","4","5","6","7","8"))

saveRDS(integrated_sub_new, "integrated_EC_Kalinski_no_pericytes.Rds")



