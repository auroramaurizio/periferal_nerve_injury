suppressMessages(library(ggbeeswarm))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
library(pheatmap)
library(reshape2)
library(scales)
library(viridis)
library(SeuratWrappers)
library(slingshot)
require(BiocStyle)
library(SingleCellExperiment)
#print(version) R version 3.6.1
#Bonanomi
packageVersion("Seurat")

############################################################
# Examine 9D post injury beads sample from Carr et al. 2019
############################################################

#load full dataset from Carr 9D post injury

raw_counts_beads <-read.table(gzfile("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/fibroblasts_Carr/GSE120678_RAW/GSM3408137_Inj_Sciatic_Beads.txt.gz"), header=TRUE)

rownames(raw_counts_beads) = raw_counts_beads$GENE

raw_counts_beadss =  raw_counts_beads[-c(1)]
head(raw_counts_beadss[1:2,1:2])

Beads <- CreateSeuratObject(counts = raw_counts_beadss, project = "Beads", min.cells = 5, min.features = 200)

Beads[["percent.mt"]] <- PercentageFeatureSet(Beads, pattern = "^mt-")
Beads[["percent.Rpl"]] <- PercentageFeatureSet(Beads, pattern = "Rpl")

Beads$stim <- "Carr_9D"

options(repr.plot.width=8, repr.plot.height=6)

VlnPlot(Beads, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 6000
max_percent_MT = 20

Beads.subset <- subset(Beads, 
                       subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


FeatureScatter(Beads.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)

FeatureScatter(Beads.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Beads.subset <- NormalizeData(Beads.subset, verbose = FALSE)
Beads.subset <- FindVariableFeatures(Beads.subset, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
Beads.subset <- ScaleData(Beads.subset, 
                          vars.to.regress = c("percent.mt", "nFeature_RNA"))

Beads.subset <- RunPCA(Beads.subset, npcs = 15, verbose = FALSE)
print(Beads.subset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Beads.subset, dims = 1:3, reduction = "pca")
DimPlot(Beads.subset, reduction = "pca")
DimHeatmap(Beads.subset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Beads.subset, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(Beads.subset)

nPC = 15 
res = 5

Beads.subset <- FindNeighbors(Beads.subset, dims = 1:nPC)
Beads.subset <- FindClusters(Beads.subset, resolution = res)
Beads.subset <- RunUMAP(Beads.subset, dims = 1:nPC)
Beads.subset <- RunTSNE(Beads.subset, dims = 1:nPC)

DimPlot(Beads.subset, reduction = "umap")
DimPlot(Beads.subset, reduction = "tsne")

FeaturePlot(Beads.subset, c("Pecam1"), order = T, reduction = "umap", label = T) 
FeaturePlot(Beads.subset, c("Cdh5"), order = T, reduction = "umap", label = T) 


#among all different clusters select only the EC ones 
EC_Carr_D9 = subset(Beads.subset, idents = c("0", "9", "18", "25", "33", "41"))


EC_Carr_D9$stim <- "Carr_9D"

saveRDS(EC_Carr_D9, "EC_Carr_9D.Rds")


