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




#function to plot gene signature sum/mean/median over UMAP.
#alternative to seurat module score

Plot_sign <- function(Seraut.object, signature, operator = sum, titlename) {
  x <- Seraut.object
  DefaultAssay(x) <- "RNA"
  
  x[["Sign_exp"]] <- apply(FetchData(object = x, 
                                     vars = signature),
                           1,
                           operator)
  FP <- FeaturePlot(x, reduction = "umap", 
                    features = 'Sign_exp', 
                    label = T, 
                    pt.size = 2,
                    order = T,
                    repel = T,
                    cols = c("lightgrey", "blue")) +
    theme(plot.title = element_text(color="blue", size=20, face="bold.italic"),
          plot.subtitle = element_text(color="dodgerblue2", size=8, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'dodgerblue4', size=7, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 10),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=10),
          axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 10),
          legend.text = element_text(face = "bold", color = "dodgerblue2", size = 6),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    #labs(title= titlename, subtitle = paste('',toString(signature), sep=''), 
    labs(title=titlename)
  #, 
  #     x = "tSNE 1", y = "tSNE 2") 
  return(FP)
}




############## Intact sample preproc


sample1_Td='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/COUNTS_reads_mapped_over_mm10_plus_Tomato/1/filtered_feature_bc_matrix'

CellRanger_sample1_Td.count <- Read10X(data.dir = sample1_Td)

S1_Td <- CreateSeuratObject(counts = CellRanger_sample1_Td.count, project = "1", min.cells = 5, min.features = 200)
S1_Td[["percent.mt"]] <- PercentageFeatureSet(S1_Td, pattern = "^mt-")
S1_Td[["percent.Rpl"]] <- PercentageFeatureSet(S1_Td, pattern = "Rpl")
S1_Td$stim <- "1"
VlnPlot(S1_Td, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 7000
max_percent_MT = 20

S1_Td.subset <- subset(S1_Td, 
                       subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)

FeatureScatter(S1_Td.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(S1_Td.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


######################### 7D post Injury sample preproc

sample2_Td='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/COUNTS_reads_mapped_over_mm10_plus_Tomato/2/filtered_feature_bc_matrix'

CellRanger_sample2_Td.count <- Read10X(data.dir = sample2_Td)
S2_Td <- CreateSeuratObject(counts = CellRanger_sample2_Td.count, project = "2", min.cells = 5, min.features = 200)

S2_Td[["percent.mt"]] <- PercentageFeatureSet(S2_Td, pattern = "^mt-")
S2_Td[["percent.Rpl"]] <- PercentageFeatureSet(S2_Td, pattern = "^Rpl-")

S2_Td$stim <- "2"

VlnPlot(S2_Td, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

min_nFeature_RNA = 200 
max_nFeature_RNA = 8000
max_percent_MT = 20

S2_Td.subset <- subset(S2_Td, 
                       subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)

FeatureScatter(S2_Td.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(S2_Td.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


########  INTEGRATED  ##########################################

#Integrate intact and injury samples and run the standard workflow for visualization and clustering

#seurat version 3.2.2

integrated_bon1T <- FindIntegrationAnchors(object.list = list(S1_Td.subset, S2_Td.subset), dims = 1:20)
features.to.integrate = integrated_bon1T@anchor.features
integrated_bon1T <- IntegrateData(anchorset = integrated_bon1T, dims = 1:20, features.to.integrate = features.to.integrate)
DefaultAssay(integrated_bon1T) <- "integrated"
integrated_bon1T <- CellCycleScoring(integrated_bon1T, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
integrated_bon1T <- ScaleData(integrated_bon1T, vars.to.regress = c("percent.mt", "nFeature_RNA", "S.Score", "G2M.Score"),  features = features.to.integrate, verbose = FALSE)
integrated_bon1T <- RunPCA(integrated_bon1T, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated_bon1T <- RunUMAP(integrated_bon1T, reduction = "pca", dims = 1:20)
integrated_bon1T <- FindNeighbors(integrated_bon1T, reduction = "pca", dims = 1:20)
integrated_bon1T <- FindClusters(integrated_bon1T, resolution = 0.3)
DimPlot(integrated_bon1T, reduction = "umap")

#FigS1 K

# check reporter gene 
FeaturePlot(integrated_bon1T, "TdTomato", order = T)

# check EC
FeaturePlot(integrated_bon1T, "Pecam1", order = T)
FeaturePlot(integrated_bon1T, "Cdh5", order = T)

# check pericytes
FeaturePlot(integrated_bon1T, "Rgs5", order = T)

# check LEC
FeaturePlot(integrated_bon1T, "Lyve1", order = T)

# check MES
FeaturePlot(integrated_bon1T, "Pdgfra", order = T)
FeaturePlot(integrated_bon1T, "Pdgfrb", order = T)

# check NK
FeaturePlot(integrated_bon1T, "Nkg7", order = T)

# check cell proliferation signature genes

dividing = c("Cdk1", "Birc1", "Mki67", "Top2a")


Plot_sign(integrated_bon1T,
          signature= dividing, 
          operator = sum, titlename = "dividing")

# assign names to macroclusters

new.cluster.ids.lit1 <- c('EC', 
                          'EC',
                          'EC',
                          'FRC',
                          'FRC',
                          'EC',
                          'EC',
                          'LEC',
                          'TCELLS')


names(new.cluster.ids.lit1) <- levels(integrated_bon1T)
renamed <- RenameIdents(integrated_bon1T, new.cluster.ids.lit1)

renamed <- RenameIdents(integrated_bon1T, new.cluster.ids.lit1)


# remove T cells 
renamed <- subset(object = renamed, idents = "TCELLS", invert = TRUE)


# FigS1 UMAPs

# FigS1 I

pdf("all_cells_integrated_newcols_newred.pdf", 12, 10)

DimPlot(renamed, cols = c('EC' = '#f56947',
                          'FIBRO' = '#9FE2BF',
                          'PERI'= '#0099CC',
                          'LEC'= '#ffcc00'), pt.size = 5) +
  theme(legend.position = "none")
dev.off()

#Fig S1 J

pdf("all_cells_integrated_group_by_stim.pdf", 12, 10)

DimPlot(renamed, group.by = "stim", pt.size = 5) +
  theme(legend.position = "none") +
  theme(plot.title = element_blank())

dev.off()


# FigS1 FeaturePlots

pdf("FP_all_cells_integrated_newcols_Lyve1.pdf")
FeaturePlot(renamed, "Lyve1", order = T, pt.size = 5)
dev.off()

pdf("FP_all_cells_integrated_newcols_Pecam1.pdf")
FeaturePlot(renamed, "Pecam1", order = T, pt.size = 5)
dev.off()

pdf("FP_all_cells_integrated_newcols_Pdgfra.pdf")
FeaturePlot(renamed, "Pdgfra", order = T, pt.size = 5)
dev.off()

pdf("FP_all_cells_integrated_newcols_Tomato.pdf")
FeaturePlot(renamed, "TdTomato", order = T, pt.size = 5)
dev.off()

pdf("FP_all_cells_integrated_newcols_Cdh5.pdf")
FeaturePlot(renamed, "Cdh5", order = T, pt.size = 5)
dev.off()

pdf("FP_all_cells_integrated_newcols_Rgs5.pdf")
FeaturePlot(renamed, "Rgs5", order = T, pt.size = 5)
dev.off()


