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


###############################################
# Examine Uninjured post injurey samples from Toma et al. 2020 GSE147285_RAW
# Peripheral Nerve Single-Cell Analysis Identifies Mesenchymal Ligands that Promote Axonal Growth
# Sciatic nerve, GSM4423509
###############################################

TomaUninj <-read.table(gzfile("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Toma/GSE147285_RAW/GSM4423509_Uninj_Sciatic.txt.gz"), header=TRUE)

rownames(TomaUninj) = TomaUninj$GENE
TomaUnin =  TomaUninj[-c(1)]

TU <- CreateSeuratObject(counts = TomaUnin, project = "TOMA_Uninj", min.cells = 5, min.features = 200)

TU[["percent.mt"]] <- PercentageFeatureSet(TU, pattern = "^mt-")
TU[["percent.Rpl"]] <- PercentageFeatureSet(TU, pattern = "Rpl")

TU$stim <- "Toma_Uninj"

VlnPlot(TU, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 5000
max_percent_MT = 20

TU.subset <- subset(TU, 
                    subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


FeatureScatter(TU.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(TU.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

TU.subset <- NormalizeData(TU.subset, verbose = FALSE)
TU.subset <- FindVariableFeatures(TU.subset, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
TU.subset <- ScaleData(TU.subset, 
                       vars.to.regress = c("percent.mt", "nFeature_RNA"))

TU.subset <- RunPCA(TU.subset, npcs = 20, verbose = FALSE)
print(TU.subset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(TU.subset, dims = 1:3, reduction = "pca")
DimPlot(TU.subset, reduction = "pca")
DimHeatmap(TU.subset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(TU.subset, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(TU.subset)

nPC = 20
res = 4

TU.subset <- FindNeighbors(TU.subset, dims = 1:nPC)
TU.subset <- FindClusters(TU.subset, resolution = res)

TU.subset <- RunUMAP(TU.subset, dims = 1:nPC)
TU.subset <- RunTSNE(TU.subset, dims = 1:nPC)

DimPlot(TU.subset, reduction = "umap")
DimPlot(TU.subset, reduction = "tsne")


#EC
FeaturePlot(TU.subset, "Pecam1", order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, c("Vwf"), order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, c("Cdh5"), order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, c("Mest"), order = T, reduction = "umap", label = T) 

#prolif
FeaturePlot(TU.subset, c("Top2a")  , order = T, reduction = "umap") 

#pericites
FeaturePlot(TU.subset, "Rgs5", order = T, reduction = "umap", label = T) 

#fibro
FeaturePlot(TU.subset, "Pdgfra", order = T, reduction = "umap", label = T) 

#schwann non myel
FeaturePlot(TU.subset, "Ngfr", order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, "L1cam", order = T, reduction = "umap", label = T) 

#schwann  myelinating
FeaturePlot(TU.subset, "Mbp", order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, "Plp1", order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, "L1cam", order = T, reduction = "umap", label = T) 


#t nk cells
FeaturePlot(TU.subset, "Cxcr6", order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, "Nkg7", order = T, reduction = "umap", label = T) 


#macrophages
FeaturePlot(TU.subset, "Aif1", order = T, reduction = "umap", label = T) 


#among all different clusters select only the EC ones 
FeaturePlot(TU.subset, "nFeature_RNA", order = T, reduction = "umap", label = T) 
FeaturePlot(TU.subset, "percent.mt", order = T, reduction = "umap", label = T) 

#among all different clusters select only the EC ones 
EC_TU.subset = subset(TU.subset, idents = c("23", "6", "24", "20", "16", "11", "12", "14", "25"))

#among all different clusters select only the Schwann ones 
TU_shwann.subset = subset(TU.subset, idents = c("0", "1", "9"))

# all Toma uninjured cells
saveRDS(TU.subset, "Toma_Uninjured.Rds")
# only EC cells among all Toma uninjured cells
saveRDS(EC_TU.subset, "EC_Toma_Uninjured.Rds")
# only Schwann cells among all Toma uninjured cells
saveRDS(TU_shwann.subset, "Toma_Uninjured_Schwann.Rds")



###############################################
# Examine 3D post injury samples from Toma et al.
# Sciatic nerve, 3 days post injury, (myelin removal beads) GSM4423506
###############################################


Toma3d <-read.table(gzfile("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Toma/GSE147285_RAW/GSM4423506_Inj_Sciatic_3d.txt.gz"), header=TRUE)

rownames(Toma3d) = Toma3d$GENE

Toma3D =  Toma3d[-c(1)]
head(Toma3D[1:2,1:2])

T3D <- CreateSeuratObject(counts = Toma3D, project = "TOMA3D", min.cells = 5, min.features = 200)


T3D[["percent.mt"]] <- PercentageFeatureSet(T3D, pattern = "^mt-")
T3D[["percent.Rpl"]] <- PercentageFeatureSet(T3D, pattern = "Rpl")

T3D$stim <- "Toma_3D"



#options(repr.plot.width=8, repr.plot.height=6)

VlnPlot(T3D, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 6000
max_percent_MT = 20

T3D.subset <- subset(T3D, 
                     subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


FeatureScatter(T3D.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)

FeatureScatter(T3D.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

T3D.subset <- NormalizeData(T3D.subset, verbose = FALSE)
T3D.subset <- FindVariableFeatures(T3D.subset, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
T3D.subset <- ScaleData(T3D.subset, 
                        vars.to.regress = c("percent.mt", "nFeature_RNA"))

T3D.subset <- RunPCA(T3D.subset, npcs = 15, verbose = FALSE)
print(T3D.subset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(T3D.subset, dims = 1:3, reduction = "pca")
DimPlot(T3D.subset, reduction = "pca")
DimHeatmap(T3D.subset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(T3D.subset, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(T3D.subset)

nPC = 15 
res = 4

T3D.subset <- FindNeighbors(T3D.subset, dims = 1:nPC)
T3D.subset <- FindClusters(T3D.subset, resolution = res)

T3D.subset <- RunUMAP(T3D.subset, dims = 1:nPC)
T3D.subset <- RunTSNE(T3D.subset, dims = 1:nPC)

DimPlot(T3D.subset, reduction = "umap")

DimPlot(T3D.subset, reduction = "tsne")

DimPlot(T3D.subset, label=T, cells.highlight=check, cols.highlight = c("cyan"), cols= "grey", reduction = "tsne")

T3D.subset$stim <- "Toma_3D"


FeaturePlot(T3D.subset, c("Mest", "Mki67","Birc5","Cdk1","Top2a")  , order = T, reduction = "tsne") 


FeaturePlot(T3D.subset, c("Mest", "Mki67","Birc5","Cdk1","Top2a")  , order = T, reduction = "umap") 


#endothelial
FeaturePlot(T3D.subset, "Cdh5", order = T, reduction = "umap", label = T) 
FeaturePlot(T3D.subset, "Cdh5", order = T, reduction = "tsne", label = T) 
FeaturePlot(T3D.subset, "Pecam1", order = T, reduction = "umap", label = T) 
FeaturePlot(T3D.subset, "Pecam1", order = T, reduction = "tsne", label = T) 

#pericites
FeaturePlot(T3D.subset, "Rgs5", order = T, reduction = "tsne", label = T) 
FeaturePlot(T3D.subset, "Rgs5", order = T, reduction = "umap", label = T) 

#fibro
FeaturePlot(T3D.subset, "Pdgfra", order = T, reduction = "tsne", label = T) 
FeaturePlot(T3D.subset, "Pdgfra", order = T, reduction = "umap", label = T) 

#schwann
FeaturePlot(T3D.subset, "Ngfr", order = T, reduction = "tsne", label = T) 
FeaturePlot(T3D.subset, "Ngfr", order = T, reduction = "umap", label = T) 

#tcells
FeaturePlot(T3D.subset, "Cxcr6", order = T, reduction = "tsne", label = T) 
FeaturePlot(T3D.subset, "Cxcr6", order = T, reduction = "umap", label = T) 

#macrophages
FeaturePlot(T3D.subset, "Aif1", order = T, reduction = "tsne", label = T) 
FeaturePlot(T3D.subset, "Aif1", order = T, reduction = "umap", label = T) 


#among all different clusters select only the EC ones 
EC_T3D.subset = subset(T3D.subset, idents = c("7", "34", "35", "9", "20"))

FeaturePlot(EC_T3D.subset, c("Mest", "Mki67","Birc5","Cdk1","Top2a")  , order = T, reduction = "umap") 

FeaturePlot(T3D.subset, "nFeature_RNA", order = T, reduction = "tsne", label = T) 

saveRDS(EC_T3D.subset, "EC_Toma_3D.Rds")
