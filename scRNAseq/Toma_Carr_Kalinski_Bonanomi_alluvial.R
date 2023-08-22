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
packageVersion("Seurat")
library(dittoSeq)
library("ggalluvial")



# upload our dataset, only EC

integrated_EC <-readRDS("integrated_EC.Rds")
integrated_EC$CellTypes = Idents(integrated_EC)
integrated_EC <- RunPCA(integrated_EC, npcs = 35, verbose = FALSE)
integrated_EC <- RunUMAP(integrated_EC, reduction = "pca", dims = 1:35, verbose = FALSE, return.model=TRUE)

sample1 <- subset(x = integrated_EC, subset = stim == "1")
sample2 <- subset(x = integrated_EC, subset = stim == "2")

# upload EC subsets of public peripheral nerve injury datasets.

# The studies:
# Carr et al. 2019 scRNAseq dataset. GSE120678
# Toma et al. 2020 dataset. GSE147285
# Kalinski et al. 2020 dataset. GSE153762

TU = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_12INTERMEDIATE7_scRNA_injury/7_bioinfo/public_data/Toma/reanalysis/Seurat_objects/EC_Toma_Uninjured.Rds")
TI = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_12INTERMEDIATE7_scRNA_injury/7_bioinfo/public_data/Toma/reanalysis/Seurat_objects/EC_Toma_3D.Rds")
KI = readRDS("/Users/maurizio.aurora/integrated_EC_Kalinski_no_pericytes.Rds")
CI = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_12INTERMEDIATE7_scRNA_injury/7_bioinfo/public_data/fibroblasts_Carr/reanalysis/Seurat_objects/EC_Carr_9D.Rds")


DefaultAssay(sample1) = "RNA"
DefaultAssay(sample2) = "RNA"
DefaultAssay(TU) = "RNA"
DefaultAssay(TI) = "RNA"
DefaultAssay(K) = "RNA"
DefaultAssay(CI) = "RNA"

sample1$stim <- "Bonanomi_Uninj"
sample2$stim <- "Bonanomi_7D"
TU$stim <- "TOMA_Uninj"
TI$stim <- "TOMA_3D"
KI$stim <- "KALINSKI_3D"
CI$stim <- "CARR_9D"

#integrate the datasets according to the Seurat standard procedure 

# normalize and find varible features in the objects
object_clean_new.list <- lapply(X = c(sample1, sample2, CI, KI, TU, TI), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# find integration anchors
integrated <- FindIntegrationAnchors(object.list = object_clean_new.list, dims = 1:15)

features.to.integrate = integrated@anchor.features

integrated <- IntegrateData(anchorset = integrated, dims = 1:15, features.to.integrate = features.to.integrate)

DefaultAssay(integrated) <- "integrated"

table(integrated$stim)

res = 2

# Run the standard workflow for visualization and clustering
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 15, verbose = FALSE)
# t-SNE and Clustering
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:15)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:15)
integrated <- FindClusters(integrated, resolution = res)
DimPlot(integrated, reduction = "umap", split.by = "stim", ncol = 3)

DimPlot(integrated, group.by = "stim", split.by = "stim", ncol = 3)

ARTERIAL <-WhichCells(object=integrated_EC, idents="ARTERIAL")
IMMATURE <-WhichCells(object=integrated_EC, idents="IMMATURE")
PROLIFERATING <-WhichCells(object=integrated_EC, idents="PROLIFERATING")
TIP <-WhichCells(object=integrated_EC, idents="TIP")
eight <-WhichCells(object=integrated_EC, idents="INTERMEDIATE")
BARR_END_CAP <-WhichCells(object=integrated_EC, idents="BARR_END_CAP")
VENOUS_PLVAP_M <-WhichCells(object=integrated_EC, idents="VENOUS_PLVAP-")
VENOUS_PLVAP_P <-WhichCells(object=integrated_EC, idents="VENOUS_PLVAP+")
CAPILLARY_PLVAP_M <-WhichCells(object=integrated_EC, idents="CAPILLARY_PLVAP-")
CAPILLARY_PLVAP_P <-WhichCells(object=integrated_EC, idents="CAPILLARY_PLVAP+")

DimPlot(integrated, label=T, cells.highlight=IMMATURE, cols.highlight = c("cyan"), cols= "grey", reduction = "umap")

#rename clusters according to cell assignments

new.cluster.ids <- c('INTERMEDIATE',
                     'IMMATURE',
                     'VENOUS_PLVAP-',
                     'VENOUS_PLVAP+',
                     'VENOUS_PLVAP+',
                     'VENOUS_PLVAP+',
                     'PROLIFERATING',
                     'ARTERIAL',
                     'CAPILLARY_PLVAP-',
                     'CAPILLARY_PLVAP+',
                     'PROLIFERATING',
                     'BARR_END_CAP',
                     'BARR_END_CAP',
                     'TIP',
                     'BARR_END_CAP',
                     'ARTERIAL',
                     'INTERMEDIATE',
                     'INTERMEDIATE',
                     'INTERMEDIATE',
                     'CAPILLARY_PLVAP-',
                     'VENOUS_PLVAP+')


names(new.cluster.ids) <- levels(integrated)
integrated_rn <- RenameIdents(integrated, new.cluster.ids)


old_names = factor(integrated_rn$stim)
levels(old_names) = c("7D","Intact","9D","3D","3D","Intact")
new.names = as.character(old_names)
integrated_rn$condition<- new.names


DimPlot(integrated_rn, label=T, cells.highlight=ARTERIAL, cols.highlight = c("cyan"), cols= "grey", reduction = "umap")


DimPlot(integrated_rn)
#select the desired assy
DefaultAssay(integrated_rn) = "integrated"

# remove intermediate cell types as they are not well carachterized
integrated_rn = subset(integrated_rn, idents = c("INTERMEDIATE"), invert = TRUE)
DimPlot(integrated_rn)
old_names = factor(integrated_rn$condition)
levels(old_names) = c("3D","7D", "9D","0_Intact")
new.names = as.character(old_names)
integrated_rn$cond<- new.names

DimPlot(integrated_rn)
#saveRDS(integrated_rn, "integrated_without_clusterINTERMEDIATE.rds")

#integrated_rn = readRDS("/Users/maurizio.aurora/integrated_without_clusterINTERMEDIATE.rds")

DimPlot(integrated_rn)
pt <- table(Idents(integrated_rn), integrated_rn$condition)
pt <- as.data.frame(pt)

pt$Cluster <- as.character(pt$Var1)
colnames(pt) = c("Clusters", "Var2", "Freq", "Cluster")
pt =pt[!grepl("INTERMEDIATE", pt$Cluster),]

pt0 =pt[grepl("Intact", pt$Var2),]
pt1 =pt[grepl("3D", pt$Var2),]
pt2 =pt[grepl("7D", pt$Var2),]
pt3 =pt[grepl("9D", pt$Var2),]

pt0$freq  = pt0$Freq/(sum(pt0$Freq)*100)
pt2$freq  = pt2$Freq/(sum(pt2$Freq)*100)
pt1$freq  = pt1$Freq/(sum(pt1$Freq)*100)
pt3$freq  = pt3$Freq/(sum(pt3$Freq)*100)


pt4 = rbind(pt0, pt1, pt2, pt3)

positions <- c("Intact", "3D", "7D", "9D")



# Fig1 G

pdf("Alluvial_EC_0_3_7_9.pdf", 15, 10)
gg <- ggplot(pt4,
             aes(x = Var2, stratum = Cluster, alluvium = Cluster,
                 y = freq,
                 fill = Cluster))  + scale_fill_manual(values=c('ARTERIAL' = '#0066FF',
                                                                'BARR_END_CAP' = '#336666',
                                                                'CAPILLARY_PLVAP-' = '#399933',
                                                                'CAPILLARY_PLVAP+' = '#99CC33',
                                                                'IMMATURE' = '#6600CC',
                                                                'PROLIFERATING' = '#FF99CC',
                                                                'TIP' = '#FF00FF',
                                                                'INTERMEDIATE' = 'grey',
                                                                'VENOUS_PLVAP-' = '#990000',
                                                                'VENOUS_PLVAP+' = '#FF6666'))

gg +
  geom_flow(alpha = 0.2) +
  geom_stratum() +
  geom_lode()+
  scale_x_discrete(limits = positions)
dev.off()