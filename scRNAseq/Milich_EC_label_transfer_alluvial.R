suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))
library(scales)
library(ggplot2)
library(ggalluvial)
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
library(dittoSeq)


# Re-analyze the dataset - Milich


uninj_sample1 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955359_qc_filtered_feature_bc_matrix_uninj_sample1.txt", header = T)
uninj_s1 <- CreateSeuratObject(counts = uninj_sample1, project = "uninj_sample1", min.cells = 5, min.features = 200)

uninj_s1[["percent.mt"]] <- PercentageFeatureSet(uninj_s1, pattern = "^mt-")
uninj_s1[["percent.Rpl"]] <- PercentageFeatureSet(uninj_s1, pattern = "Rpl")

uninj_s1$stim <- "uninj"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(uninj_s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##########
uninj_sample2 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955360_qc_filtered_feature_bc_matrix_uninj_sample2.txt", header = T)
uninj_s2 <- CreateSeuratObject(counts = uninj_sample2, project = "uninj_sample2", min.cells = 5, min.features = 200)


uninj_s2[["percent.mt"]] <- PercentageFeatureSet(uninj_s2, pattern = "^mt-")
uninj_s2[["percent.Rpl"]] <- PercentageFeatureSet(uninj_s2, pattern = "Rpl")

uninj_s2$stim <- "uninj"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(uninj_s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##########

uninj_sample3 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955361_qc_filtered_feature_bc_matrix_uninj_sample3.txt", header = T)
uninj_s3 <- CreateSeuratObject(counts = uninj_sample3, project = "uninj_sample3", min.cells = 5, min.features = 200)

uninj_s3[["percent.mt"]] <- PercentageFeatureSet(uninj_s3, pattern = "^mt-")
uninj_s3[["percent.Rpl"]] <- PercentageFeatureSet(uninj_s3, pattern = "Rpl")
uninj_s3$stim <- "uninj"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(uninj_s3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##########

Onedpi_sample1 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955362_qc_filtered_feature_bc_matrix_1dpi_sample1.txt", header = T)
Onedpi_s1 <- CreateSeuratObject(counts = Onedpi_sample1, project = "1dpi_s1", min.cells = 5, min.features = 200)


Onedpi_s1[["percent.mt"]] <- PercentageFeatureSet(Onedpi_s1, pattern = "^mt-")
Onedpi_s1[["percent.Rpl"]] <- PercentageFeatureSet(Onedpi_s1, pattern = "Rpl")

Onedpi_s1$stim <- "1dpi"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Onedpi_s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##########

Onedpi_sample2 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955363_qc_filtered_feature_bc_matrix_1dpi_sample2.txt", header = T)
Onedpi_s2 <- CreateSeuratObject(counts = Onedpi_sample2, project = "1dpi_s2", min.cells = 5, min.features = 200)


Onedpi_s2[["percent.mt"]] <- PercentageFeatureSet(Onedpi_s2, pattern = "^mt-")
Onedpi_s2[["percent.Rpl"]] <- PercentageFeatureSet(Onedpi_s2, pattern = "Rpl")

Onedpi_s2$stim <- "1dpi"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Onedpi_s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



##########

Onedpi_sample3 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955364_qc_filtered_feature_bc_matrix_1dpi_sample3.txt", header = T)
Onedpi_s3 <- CreateSeuratObject(counts = Onedpi_sample3, project = "1dpi_s3", min.cells = 5, min.features = 200)


Onedpi_s3[["percent.mt"]] <- PercentageFeatureSet(Onedpi_s3, pattern = "^mt-")
Onedpi_s3[["percent.Rpl"]] <- PercentageFeatureSet(Onedpi_s3, pattern = "Rpl")

Onedpi_s3$stim <- "1dpi"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Onedpi_s3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##########

Threedpi_sample1 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955365_qc_filtered_feature_bc_matrix_3dpi_sample1.txt", header = T)
Threedpi_s1 <- CreateSeuratObject(counts = Threedpi_sample1, project = "3dpi_s1", min.cells = 5, min.features = 200)


Threedpi_s1[["percent.mt"]] <- PercentageFeatureSet(Threedpi_s1, pattern = "^mt-")
Threedpi_s1[["percent.Rpl"]] <- PercentageFeatureSet(Threedpi_s1, pattern = "Rpl")

Threedpi_s1$stim <- "3dpi"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Threedpi_s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

##########

Threedpi_sample2 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955366_qc_filtered_feature_bc_matrix_3dpi_sample2.txt", header = T)
Threedpi_s2 <- CreateSeuratObject(counts = Threedpi_sample2, project = "3dpi_s2", min.cells = 5, min.features = 200)


Threedpi_s2[["percent.mt"]] <- PercentageFeatureSet(Threedpi_s2, pattern = "^mt-")
Threedpi_s2[["percent.Rpl"]] <- PercentageFeatureSet(Threedpi_s2, pattern = "Rpl")

Threedpi_s2$stim <- "3dpi"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Threedpi_s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#########

Sevendpi_sample1 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955367_qc_filtered_feature_bc_matrix_7dpi_sample1.txt", header = T)
Sevendpi_s1 <- CreateSeuratObject(counts = Sevendpi_sample1, project = "7dpi_s1", min.cells = 5, min.features = 200)

Sevendpi_s1[["percent.mt"]] <- PercentageFeatureSet(Sevendpi_s1, pattern = "^mt-")
Sevendpi_s1[["percent.Rpl"]] <- PercentageFeatureSet(Sevendpi_s1, pattern = "Rpl")

Sevendpi_s1$stim <- "7dpi"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Sevendpi_s1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##########


Sevendpi_sample2 <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610/GSM4955368_qc_filtered_feature_bc_matrix_7dpi_sample2.txt", header = T)
Sevendpi_s2 <- CreateSeuratObject(counts = Sevendpi_sample2, project = "7dpi_s2", min.cells = 5, min.features = 200)

Sevendpi_s2[["percent.mt"]] <- PercentageFeatureSet(Sevendpi_s2, pattern = "^mt-")
Sevendpi_s2[["percent.Rpl"]] <- PercentageFeatureSet(Sevendpi_s2, pattern = "Rpl")

Sevendpi_s2$stim <- "7dpi"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Sevendpi_s2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##########


# normalize and find varible features in the objects
object_clean_new.list <- lapply(X = c(uninj_s1, uninj_s2, uninj_s3, Onedpi_s1, Onedpi_s2, Onedpi_s3, Threedpi_s1, Threedpi_s2, Sevendpi_s1, Sevendpi_s2), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# find integration anchors
integrated_milich <- FindIntegrationAnchors(object.list = object_clean_new.list, dims = 1:20, anchor.features = 3000)

features.to.integrate = integrated_milich@anchor.features

integrated_milich <- IntegrateData(anchorset = integrated_milich, dims = 1:20, features.to.integrate = features.to.integrate)

table(integrated_milich$stim)
res = 0.5
# Run the standard workflow for visualization and clustering
integrated_milich <- ScaleData(integrated_milich, verbose = FALSE)
integrated_milich <- RunPCA(integrated_milich, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated_milich <- RunUMAP(integrated_milich, reduction = "pca", dims = 1:20)
#integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:20)
integrated_milich <- FindNeighbors(integrated_milich, reduction = "pca", dims = 1:20)
integrated_milich <- FindClusters(integrated_milich, resolution = res)
DimPlot(integrated_milich, reduction = "umap", split.by = "stim", ncol = 3)
#DimPlot(integrated_bon1, reduction = "tsne")

var_genes <- VariableFeatures(integrated_milich)
as(as.matrix(integrated_milich@assays$RNA@scale.data), 'sparseMatrix')
DimPlot(integrated_milich, reduction = "umap")
DefaultAssay(integrated_milich) = "RNA"


setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/")

saveRDS(integrated_milich, "Milich_only_integrated_stim.Rds")


Milich = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/Milich_only_integrated_stim.Rds")
Metadata = read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/GSE162610_barcode_metadata.tsv", sep = "\t")


Milich_EC <- subset(Milich, cells = cId_Endo_Meta)

DimPlot(Milich_EC, split.by = "stim", ncol = 2)

DefaultAssay(Milich_EC) = "RNA"
table(Milich_EC$stim)

D0 <- subset(x = Milich_EC, subset = stim == "1_uninj")
D1 <- subset(x = Milich_EC, subset = stim == "1dpi")
D3 <- subset(x = Milich_EC, subset = stim == "3dpi")
D7 <- subset(x = Milich_EC, subset = stim == "7dpi")

object_clean_new1.list <- lapply(X = c(D0,D1,D3,D7), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# find integration anchors
integrated <- FindIntegrationAnchors(object.list = object_clean_new1.list, dims = 1:30)

features.to.integrate = integrated@anchor.features

integrated <- IntegrateData(anchorset = integrated, dims = 1:30, features.to.integrate = features.to.integrate)

DefaultAssay(integrated) <- "integrated"

res = 0.1

DefaultAssay(integrated) = "integrated"
# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30)
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:30)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30)
integrated <- FindClusters(integrated, resolution = res)



DimPlot(integrated, split.by = "stim", ncol = 3)

DimPlot(integrated)

#Milich_only_integrated_stim = Milich_only_integrated_stimEC
#DimPlot(Milich_only_integrated_stim)

#Milich_only_integrated_stimEC = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Milich/Milich_only_integrated_stimEC.Rds")
#DimPlot(Milich_only_integrated_stimEC, split.by = "orig.ident", ncol = 2)
#DefaultAssay(Milich_only_integrated_stimEC)
EC1 = readRDS("/Users/maurizio.aurora/integrated_without_cluster8.rds")

#integrated = Milich_only_integrated_stimEC
DefaultAssay(integrated) = "integrated"

TI.anchors <- FindTransferAnchors(reference = EC1, query = integrated,
                                  dims = 1:30)

predictions <- TransferData(anchorset = TI.anchors, refdata = Idents(EC1),
                            dims = 1:30)

table(predictions$predicted.id)



Milich_only_integrated_stimEC_copy = integrated

Milich_only_integrated_stimEC_copy <- AddMetaData(Milich_only_integrated_stimEC_copy, metadata = predictions)

Idents(Milich_only_integrated_stimEC_copy) = Milich_only_integrated_stimEC_copy$predicted.id


pt <- table(Idents(Milich_only_integrated_stimEC_copy), Milich_only_integrated_stimEC_copy$stim)
pt <- as.data.frame(pt)
pt$Cluster <- as.character(pt$Var1)
pt$Group=paste(pt$Cluster,pt$Var2)
pt$Cluster <- as.character(pt$Var1)
colnames(pt) = c("Clusters", "Var2", "Freq", "Cluster")

pt0 =pt[grepl("1_uninj", pt$Var2),]

pt1 =pt[grepl("1dpi", pt$Var2),]

pt2 =pt[grepl("3dpi", pt$Var2),]

pt3 =pt[grepl("7dpi", pt$Var2),]

pt0$freq  = pt0$Freq/(sum(pt0$Freq)*100)
pt2$freq  = pt2$Freq/(sum(pt2$Freq)*100)
pt1$freq  = pt1$Freq/(sum(pt1$Freq)*100)
pt3$freq  = pt3$Freq/(sum(pt3$Freq)*100)


pt4 = rbind(pt0, pt1, pt2, pt3)

positions <- c("1_uninj", "1dpi", "3dpi", "7dpi")

unique(pt4$Clusters)


ggplot(pt4, aes(x = Var2, y = freq, fill = Cluster)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Sample") +
  ylab("Proportion") 


pdf("barplot_Milich_cond_no_eight_integrated_alluvial_only_their_cells_test.pdf", 15, 10)
gg <- ggplot(pt4,
             aes(x = Var2, stratum = Clusters, alluvium = Clusters,
                 y = freq,
                 fill = Clusters))  + scale_fill_manual(values=c('ARTERIAL' = '#0066FF',
                                                                 'BARR_END_CAP' = '#336666', 
                                                                 'CAPILLARY_PLVAP-' = '#399933',
                                                                 'CAPILLARY_PLVAP+' = '#99CC33',
                                                                 'TIP_1' = '#6600CC',
                                                                 'TIP_2' = '#FF99CC',
                                                                 'TIP_3' = '#FF00FF',
                                                                 'VENOUS_PLVAP-' = '#990000',
                                                                 'VENOUS_PLVAP+' = '#FF6666')) 

gg + 
  geom_flow(alpha = 0.2) + 
  geom_stratum() + 
  geom_lode()+
  scale_x_discrete(limits = positions)#+



dev.off()


