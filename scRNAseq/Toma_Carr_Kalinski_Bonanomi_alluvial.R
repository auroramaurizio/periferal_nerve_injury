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
library("ggalluvial")

###


integrated_3TIP <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3TIP/integrated_EC_3TIP.Rds")

DefaultAssay(integrated_3TIP) = "RNA"

DefaultAssay(integrated_3TIP) = "integrated"
res= 1
subs <- FindClusters(integrated_3TIP, resolution = res)

DimPlot(subs, split.by = "stim")

new.cluster.ids.lit <- c('VENOUS_PLVAP+',
                         'TIP_1',
                         'BARR_END_CAP',
                         'VENOUS_PLVAP+',
                         'ARTERIAL',
                         'CAPILLARY_PLVAP-',
                         'TIP_2',
                         'VENOUS_PLVAP-',
                         '8',
                         'TIP_3',
                         'CAPILLARY_PLVAP+'
)


names(new.cluster.ids.lit) <- levels(subs)
integrated_3TIP_new_eight <- RenameIdents(subs, new.cluster.ids.lit)

DimPlot(integrated_3TIP_new_eight)

DefaultAssay(integrated_3TIP_new_eight) = "integrated"

integrated = subset(integrated_3TIP_new_eight, idents = c("8"), invert = TRUE)

DimPlot(integrated)


table(Idents(integrated_3TIP_new_eight))

integrated_3TIP_new_eight$CellTypes = Idents(integrated_3TIP_new_eight)
integrated_3TIP_new_eight <- RunPCA(integrated_3TIP_new_eight, npcs = 35, verbose = FALSE)
integrated_3TIP_new_eight <- RunUMAP(integrated_3TIP_new_eight, reduction = "pca", dims = 1:35, verbose = FALSE, return.model=TRUE)

sample1 <- subset(x = integrated_3TIP_new_eight, subset = stim == "1")
sample2 <- subset(x = integrated_3TIP_new_eight, subset = stim == "2")

table(Idents(sample1))
table(Idents(sample2))

DimPlot(TU)


FeaturePlot(TU, "Cdh5", order = T, label = T)
TU = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Toma/reanalysis/Seurat_objects/EC_Toma_Uninjured.Rds")
TI = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Toma/reanalysis/Seurat_objects/EC_Toma_3D.Rds")
KI = readRDS("/Users/maurizio.aurora/integrated_EC_Kalinski_no_pericytes.Rds")
CI = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/fibroblasts_Carr/reanalysis/Seurat_objects/EC_Carr_9D.Rds")



FeaturePlot(CI, "Cdh5")
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

DimPlot(integrated)

#PERICYTES <-WhichCells(object=sample2, idents="PERI")
ARTERIAL <-WhichCells(object=integrated_3TIP_new_eight, idents="ARTERIAL")
TIP_1 <-WhichCells(object=integrated_3TIP_new_eight, idents="TIP_1")
TIP_2 <-WhichCells(object=integrated_3TIP_new_eight, idents="TIP_2")
TIP_3 <-WhichCells(object=integrated_3TIP_new_eight, idents="TIP_3")
eight <-WhichCells(object=integrated_3TIP_new_eight, idents="8")
BARR_END_CAP <-WhichCells(object=integrated_3TIP_new_eight, idents="BARR_END_CAP")
VENOUS_PLVAP_M <-WhichCells(object=integrated_3TIP_new_eight, idents="VENOUS_PLVAP-")
VENOUS_PLVAP_P <-WhichCells(object=integrated_3TIP_new_eight, idents="VENOUS_PLVAP+")
CAPILLARY_PLVAP_M <-WhichCells(object=integrated_3TIP_new_eight, idents="CAPILLARY_PLVAP-")
CAPILLARY_PLVAP_P <-WhichCells(object=integrated_3TIP_new_eight, idents="CAPILLARY_PLVAP+")

DimPlot(integrated, label=T, cells.highlight=TIP_1, cols.highlight = c("cyan"), cols= "grey", reduction = "umap")

new.cluster.ids <- c('8',
                     'TIP_1',
                     'VENOUS_PLVAP-',
                     'VENOUS_PLVAP+',
                     'VENOUS_PLVAP+',
                     'VENOUS_PLVAP+',
                     'TIP_2',
                     'ARTERIAL',
                     'CAPILLARY_PLVAP-',
                     'CAPILLARY_PLVAP+',
                     'TIP_2',
                     'BARR_END_CAP',
                     'BARR_END_CAP',
                     'TIP_3',
                     'BARR_END_CAP',
                     'ARTERIAL',
                     '8',
                     '8',
                     '8',
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

integrated_rn = subset(integrated_rn, idents = c("8"), invert = TRUE)

table(integrated_rn$stim)

DimPlot(integrated_rn)
old_names = factor(integrated_rn$condition)
levels(old_names) = c("3D","7D", "9D","0_Intact")
new.names = as.character(old_names)
integrated_rn$cond<- new.names

DimPlot(integrated_rn)
#saveRDS(integrated_rn, "integrated_without_cluster8_.rds")

#integrated_rn = readRDS("/Users/maurizio.aurora/integrated_without_cluster8.rds")

DimPlot(integrated_rn)
pt <- table(Idents(integrated_rn), integrated_rn$condition)
pt <- as.data.frame(pt)

pt$Cluster <- as.character(pt$Var1)
colnames(pt) = c("Clusters", "Var2", "Freq", "Cluster")
pt =pt[!grepl("8", pt$Cluster),]

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
                                                                'TIP_1' = '#6600CC',
                                                                'TIP_2' = '#FF99CC',
                                                                'TIP_3' = '#FF00FF',
                                                                '8' = 'grey',
                                                                'VENOUS_PLVAP-' = '#990000',
                                                                'VENOUS_PLVAP+' = '#FF6666'))

gg +
  geom_flow(alpha = 0.2) +
  geom_stratum() +
  geom_lode()+
  scale_x_discrete(limits = positions)
dev.off()