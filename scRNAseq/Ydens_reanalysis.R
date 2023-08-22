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
library(SingleCellExperiment)
#print(version) R version 3.6.1
#packageVersion("Seurat")
sessionInfo()


###########################################
# Examine the three datasets from Ydens et al.
# steady state
# day1 post injury
# day5 post injury
###########################################


#ststD1D5_rawcounts <-read.table("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/macrophages/GSE144707_countTable_aggrNerveStStD1D5.txt", header=T)

##########
# steady state
##########


stst_rawcounts <-read.table(gzfile("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/macrophages_Ydens/GSE144707_RAW/GSM4294086_countTable_nerveStSt.txt.gz"), header = T)

rownames(stst_rawcounts) = stst_rawcounts$gene
stst_rawcounts$gene <- NULL

steady <- CreateSeuratObject(counts = stst_rawcounts, project = "stst", min.cells = 5, min.features = 200)

steady[["percent.mt"]] <- PercentageFeatureSet(steady, pattern = "^mt-")
steady[["percent.Rpl"]] <- PercentageFeatureSet(steady, pattern = "Rpl")

steady$stim <- "steady"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(steady, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


min_nFeature_RNA = 300 
max_nFeature_RNA = 2000
max_percent_MT = 3

steady.subset <- subset(steady, 
             subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


VlnPlot(steady.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(steady.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(steady.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(steady.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



##########
# D1
##########

D1_rawcounts <-read.table(gzfile("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/macrophages_Ydens/GSE144707_RAW/GSM4294087_countTable_nerveD1.txt.gz"), header = T)

rownames(D1_rawcounts) = D1_rawcounts$gene
D1_rawcounts$gene <- NULL

D1 <- CreateSeuratObject(counts = D1_rawcounts, project = "D1", min.cells = 5, min.features = 200)


D1[["percent.mt"]] <- PercentageFeatureSet(D1, pattern = "^mt-")
D1[["percent.Rpl"]] <- PercentageFeatureSet(D1, pattern = "Rpl")

D1$stim <- "D1"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(D1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


min_nFeature_RNA = 300 
max_nFeature_RNA = 2000
max_percent_MT = 3

D1.subset <- subset(D1, 
             subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


FeatureScatter(D1.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(D1.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


VlnPlot(D1.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


##########
# D5
##########

D5_rawcounts <-read.table(gzfile("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/macrophages_Ydens/GSE144707_RAW/GSM4294088_countTable_nerveD5.txt.gz"), header = T)

rownames(D5_rawcounts) = D5_rawcounts$gene
D5_rawcounts$gene <- NULL

D5 <- CreateSeuratObject(counts = D5_rawcounts, project = "D5", min.cells = 5, min.features = 200)

D5[["percent.mt"]] <- PercentageFeatureSet(D5, pattern = "^mt-")
D5[["percent.Rpl"]] <- PercentageFeatureSet(D5, pattern = "Rpl")

D5$stim <- "D5"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(D5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


min_nFeature_RNA = 300 
max_nFeature_RNA = 2000
max_percent_MT = 3

D5.subset <- subset(D5, 
             subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


FeatureScatter(D5.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(D5.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(D5.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)





###################
# merge the objects
###################

ststD1D5 <- merge(steady, y = c(D1, D5), add.cell.ids = c("stst", "D1", "D5"), project = "ssD1D5")

table(ststD1D5$orig.ident)
#D1   D5 stst 
#2793 2059  728

ststD1D5[["percent.mt"]] <- PercentageFeatureSet(ststD1D5, pattern = "^mt-")
ststD1D5[["percent.Rpl"]] <- PercentageFeatureSet(ststD1D5, pattern = "Rpl")

#ststD1D5$stim <- "ssD1D5"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(ststD1D5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


min_nFeature_RNA = 300 
max_nFeature_RNA = 2000
max_percent_MT = 2.5

ststD1D5.subset <- subset(ststD1D5, 
             subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


VlnPlot(ststD1D5.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

FeatureScatter(ststD1D5.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(ststD1D5.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

VlnPlot(ststD1D5.subset, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


ststD1D5.subset <- NormalizeData(ststD1D5.subset, verbose = FALSE)
ststD1D5.subset <- FindVariableFeatures(ststD1D5.subset, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
ststD1D5.subset <- ScaleData(ststD1D5.subset, 
                vars.to.regress = c("percent.mt", "nFeature_RNA"))

ststD1D5.subset <- RunPCA(ststD1D5.subset, npcs = 15, verbose = FALSE)
print(ststD1D5.subset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ststD1D5.subset, dims = 1:3, reduction = "pca")
DimPlot(ststD1D5.subset, reduction = "pca")
DimHeatmap(ststD1D5.subset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(ststD1D5.subset, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(ststD1D5.subset)

nPC = 15 
res = 0.5

ststD1D5.subset <- FindNeighbors(ststD1D5.subset, dims = 1:nPC)
ststD1D5.subset <- FindClusters(ststD1D5.subset, resolution = res)

ststD1D5.subset <- RunUMAP(ststD1D5.subset, dims = 1:nPC)
ststD1D5.subset <- RunTSNE(ststD1D5.subset, dims = 1:nPC)

DimPlot(ststD1D5.subset, reduction = "umap")


#library(scales)
#show_col(hue_pal()(3))

DimPlot(ststD1D5.subset, reduction = "tsne", group.by = "orig.ident", cols = c('stst' = '#619CFF', 'D1' = '#00BA38',  'D5' = '#F8766D')) + ggtitle('Steady state + Day 1 + Day 5')

DimPlot(ststD1D5.subset, reduction = "tsne", group.by = "orig.ident", cols = c('stst' = '#619CFF', 'D1' = 'lightgrey',  'D5' = 'lightgrey')) + ggtitle('Steady state')

DimPlot(ststD1D5.subset, reduction = "tsne", group.by = "orig.ident", cols = c('stst' = 'lightgrey', 'D1' = '#00BA38',  'D5' = 'lightgrey')) + ggtitle('Day 1')

DimPlot(ststD1D5.subset, reduction = "tsne", group.by = "orig.ident", cols = c('stst' = 'lightgrey', 'D1' = 'lightgrey',  'D5' = '#F8766D')) + ggtitle('Day 5')



FeaturePlot(ststD1D5.subset, c("Mgl2", "Retnla", "Arg1", "Vegfa", "Ccl12", "Clec10a","S100a4","H2-Aa"), order = TRUE , reduction = "tsne", ncol = 4, label = FALSE ) & NoLegend()


DimPlot(ststD1D5.subset, reduction = "tsne")

thresh.use = 0.25
min.pct = 0.25
min.diff.pct = -Inf
test.use ="wilcox"

cluster.markers = FindAllMarkers(ststD1D5.subset, thresh.use = thresh.use, test.use=test.use, min.pct=min.pct, min.diff.pct=min.diff.pct, only.pos=TRUE)


top10 <- cluster.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
DoHeatmap(ststD1D5.subset, features = top10$gene,angle = 45, size = 2) + NoLegend()

Epineurial_macrophages= c("Ccl8", "Ccl24", "Cd209", "Clec10a", "Fcna", "Folr2", "Fxyd2", "Retnla", "Tslp")

Endoneurial_macrophages=c("Ccr2","Il1rl1","Cxcl1","Selm","Pla2g2d","Qpct","Tnfsf9","Xist","H2-Aa","H2-Ab1","Cd74")


Plot_sign <- function(Seraut.object, signature, operator = sum, titlename) {
  x <- Seraut.object
  DefaultAssay(x) <- "RNA"
  
  x[["Sign_exp"]] <- apply(FetchData(object = x, 
                                     vars = signature),
                           1,
                           operator)
  FP <- FeaturePlot(x, reduction = "tsne", 
                    features = 'Sign_exp', 
                    label = T, 
                    pt.size = 2,
                    order = T,
                    split.by="stim",
                    cols = c("lightgrey", "blue")) +
    theme(plot.title = element_text(color="blue", size=15, face="bold.italic"),
          plot.subtitle = element_text(color="dodgerblue2", size=8, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'dodgerblue4', size=10, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 10),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=10),
          axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 10),
          legend.text = element_text(face = "bold", color = "dodgerblue2", size = 10),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(title=titlename)
  return(FP)
}

Plot_sign(ststD1D5.subset,
                     signature= Epineurial_macrophages, 
                     operator = median, titlename = "Epineurial macrophages")


Plot_sign(ststD1D5.subset,
                     signature= Endoneurial_macrophages, 
                     operator = median, titlename = "Endoneurial macrophages")


DimPlot(ststD1D5.subset, split.by= "orig.ident", ncol = 2, reduction = "tsne")

DimPlot(ststD1D5.subset)

new.cluster.ids <- c('D1', 
                     'D5',
                     'D5',
                     'D1',
                     'stst',
                     'stst',
                     'D1',
                     'D5')

names(new.cluster.ids) <- levels(ststD1D5.subset)
ststD1D5.subset_rn <- RenameIdents(ststD1D5.subset, new.cluster.ids)

DimPlot(ststD1D5.subset_rn, reduction = "tsne", cols = c('stst' = '#619CFF', 'D1' = '#00BA38',  'D5' = '#F8766D')) + ggtitle('Steady state + Day 1 + Day 5')

DimPlot(ststD1D5.subset_rn, group.by = "stim", cols = c('steady' = '#619CFF', 'D1' = '#00BA38',  'D5' = '#F8766D'), reduction = "umap") + ggtitle('Steady state + Day 1 + Day 5')

saveRDS(ststD1D5.subset_rn, "Ydens_stst_D1_D5_merged.Rds")

