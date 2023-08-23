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

############################################################
# Examine 9D post injury beads sample from Carr et al. 2019
############################################################



# Carr's dataset comes with counts data, as well as a shiny app and Rdata.

# shiny app: https://github.com/millerkaplanlab/MouseSciaticNerve


# Still we reanalized it from scratch obtaining matching results.

###########################################
# Examine uninjured sample from Carr et al.
###########################################

raw_counts_uninj <-read.table(gzfile("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/fibroblasts_Carr/GSE120678_RAW/GSM3408139_Uninj_Sciatic_Mesenchymal.txt.gz"))

#raw_counts_uninj[grep("^Pdgfr", row.names(raw_counts_uninj)),]

Uninj <- CreateSeuratObject(counts = raw_counts_uninj, project = "Uninj", min.cells = 5, min.features = 200)

Uninj[["percent.mt"]] <- PercentageFeatureSet(Uninj, pattern = "^mt-")
Uninj[["percent.Rpl"]] <- PercentageFeatureSet(Uninj, pattern = "Rpl")

Uninj$stim <- "Uninj"

options(repr.plot.width=8, repr.plot.height=6)
VlnPlot(Uninj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


min_nFeature_RNA = 200 
max_nFeature_RNA = 7500
max_percent_MT = 5

Uninj.subset <- subset(Uninj, 
                       subset = nFeature_RNA > min_nFeature_RNA & nFeature_RNA < max_nFeature_RNA & percent.mt < max_percent_MT)


FeatureScatter(Uninj.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(Uninj.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

Uninj.subset <- NormalizeData(Uninj.subset, verbose = FALSE)
Uninj.subset <- FindVariableFeatures(Uninj.subset, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
Uninj.subset <- ScaleData(Uninj.subset, 
                          vars.to.regress = c("percent.mt", "nFeature_RNA"))

Uninj.subset <- RunPCA(Uninj.subset, npcs = 15, verbose = FALSE)
print(Uninj.subset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Uninj.subset, dims = 1:3, reduction = "pca")
DimPlot(Uninj.subset, reduction = "pca")
DimHeatmap(Uninj.subset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(Uninj.subset, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(Uninj.subset)

nPC = 15 
res = 0.5

Uninj.subset <- FindNeighbors(Uninj.subset, dims = 1:nPC)
Uninj.subset <- FindClusters(Uninj.subset, resolution = res)

Uninj.subset <- RunUMAP(Uninj.subset, dims = 1:nPC)
Uninj.subset <- RunTSNE(Uninj.subset, dims = 1:nPC)

DimPlot(Uninj.subset, reduction = "umap")

DimPlot(Uninj.subset, reduction = "tsne")

Uninj.subset #679 #879
#(B–E) t-SNEs overlaid for expression of (B) Pdgfra, Pdgfrb, and Col1a1; (C) the perineurial gene Itgb4; (D) the epineurial genes Pcolce2, Ly6c1, Comp, and Cxcl13;
#and (E) the endoneurial genes Sox9, Etv1, Col26a1, and Osr2. Expression levels are indicated according to the color keys.


features <- c("Pdgfra", "Pdgfrb", "Col1a1","Itgb4")
epineurial_genes <- c("Pcolce2","Ly6c1","Comp", "Cxcl13")
endoneurial_genes <- c("Sox9","Etv1","Col26a1", "Osr2", "Wif1")
perineurial_genes <- c("Itgb4")
dividing_genes <- c("Top2a")
differentiating_genes <- c("Dlk1","Tnc","Plagl1")


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
                    cols = c("lightgrey", "blue")) +
    theme(plot.title = element_text(color="blue", size=10, face="bold.italic"),
          plot.subtitle = element_text(color="dodgerblue2", size=8, face="italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = 'dodgerblue4', size=10, hjust =1), 
          axis.title.x = element_text(face = "bold", color = "dodgerblue2", size = 10),
          axis.text.y = element_text(angle = 0, face = "bold", color = 'dodgerblue4', size=10),
          axis.title.y = element_text(face = "bold", color = "dodgerblue2", size = 10),
          legend.text = element_text(face = "bold", color = "dodgerblue2", size = 10),
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    #labs(title= titlename, subtitle = paste('',toString(signature), sep=''), 
    labs(title=titlename)
  #, 
  #     x = "tSNE 1", y = "tSNE 2") 
  return(FP)
}



p1 <- Plot_sign(Uninj.subset,
                signature= epineurial_genes, 
                operator = median, titlename = "Epineurial genes")

p2 <- Plot_sign(Uninj.subset,
                signature= endoneurial_genes, 
                operator = median, titlename = "Endoneurial genes")

p3 <- Plot_sign(Uninj.subset,
                signature= perineurial_genes, 
                operator = median, titlename = "Perineurial genes")

p4 <- Plot_sign(Uninj.subset,
                signature= differentiating_genes, 
                operator = median, titlename = "Differentiating genes")

p5 <- Plot_sign(Uninj.subset,
                signature= dividing_genes, 
                operator = median, titlename = "Dividing genes")

differentiating
plot_grid(p1, p2, p3, nrow = 2)




#FindMarkers will find markers between two different identity groups - you have to specify both identity groups. This is useful for comparing the differences between two specific groups.

#FindAllMarkers will find markers differentially expressed in each identity group by comparing it to all of the others - you don't have to manually define anything. Note that markers may bleed over between closely-related groups - they are not forced to be specific to only one group. This is what most people use (and likely what you want).

#FindConservedMarkers will find markers that are conserved between two groups - this can be useful if you want to find markers that are conserved between a treated and untreated condition for a specific cell type or group of cells. It means they are differentially expressed compared to other groups, but have similar expression between the two groups you're actually comparing.

thresh.use = 0.25
min.pct = 0.25
min.diff.pct = -Inf
test.use ="wilcox"

cluster.markers = FindAllMarkers(Uninj.subset, thresh.use = thresh.use, test.use=test.use, min.pct=min.pct, min.diff.pct=min.diff.pct, only.pos=TRUE)

levels(Uninj.subset)

new.cluster.ids <- c('MES_ENDONEUR',
                     'MES_EPINEUR',
                     'MES_ENDONEUR',
                     'MES_EPINEUR'
                     'MES_PERINEUR')

Uninj.subset$stim <- "Carr_intact"

names(new.cluster.ids) <- levels(Uninj.subset)
Uninj.subset_rn <- RenameIdents(Uninj.subset, new.cluster.ids)

saveRDS(Uninj.subset_rn, "MES_Carr_intact.Rds")



###########################################
# Examine Injured sample from Carr et al.
###########################################


#load full dataset from Carr D9post injury
#GSM3408137_Inj_Sciatic_Beads


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

# Proliferative cell clusters are present in all clusters: EC, mesenchymal cells 
# and Shwann cells after injury as observed in other studies.

FeaturePlot(Beads.subset, c("Mest", "Mki67","Birc5","Cdk1","Top2a"), order = T, reduction = "umap", label = T)

# schwann cells express "Ngfr" and "Sox10 ...
FeaturePlot(Beads.subset, c("Ngfr", "Sox10"), order = T, reduction = "umap", label = T)

# mesenchymal cells express "Pdgfra" ...
FeaturePlot(Beads.subset, c("Pdgfra"), order = T, reduction = "umap", label = T)

# endothelial cells express "Pecam1", "Cdh5", "Cd31", "Esam" ...
FeaturePlot(Beads.subset, c("Pecam1", "Cdh5", "Cd31", "Esam"), order = T, reduction = "umap", label = T)

# macrophages express "Cd68","Aif1", "Emr1" ...
FeaturePlot(Beads.subset, c("Cd68","Aif1", "Emr1"), order = T, reduction = "umap", label = T)



#among all different clusters select only the EC ones 
EC_Carr_D9 = subset(Beads.subset, idents = c("0", "9", "18", "25", "33", "41"))

EC_Carr_D9$stim <- "Carr_9D"

saveRDS(EC_Carr_D9, "EC_Carr_9D.Rds")


#among all different clusters select only the Schwann ones 
Schwann_Carr_D9 = subset(Beads.subset, idents = c("20", "21", "10", "23", "39", "34"))

Schwann_Carr_D9$stim <- "Carr_9D"

saveRDS(Schwann_Carr_D9, "Schwann_Carr_D9.Rds")
