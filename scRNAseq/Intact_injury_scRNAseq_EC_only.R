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

sample1_Td='/Users/maurizio.aurora/Documents/pep/1_TdTomatoBonanomi/outs/filtered_feature_bc_matrix'

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

S1_Td.subset <- CellCycleScoring(S1_Td.subset, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

FeatureScatter(S1_Td.subset, feature1 = "nCount_RNA", feature2 = "percent.mt",cols = NULL)
FeatureScatter(S1_Td.subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

all.genesS1_Td.subset <- rownames(S1_Td.subset)

S1_Td.subset <- NormalizeData(S1_Td.subset, verbose = FALSE)
S1_Td.subset <- FindVariableFeatures(S1_Td.subset, selection.method = "vst", nfeatures = 2000)
S2_Td.subset <- CellCycleScoring(S2_Td.subset, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

#ScaleDate and regress for MTpercent and nFeature_RNA
S1_Td.subset <- ScaleData(S1_Td.subset,
                          vars.to.regress = c("percent.mt", "nFeature_RNA", "S.Score", "G2M.Score"),
                          features = all.genesS1_Td.subset)


S1_Td.subset <- RunPCA(S1_Td.subset, npcs = 30, verbose = FALSE)
S1_Td.subset <- FindNeighbors(S1_Td.subset, dims = 1:nPC)
S1_Td.subset <- FindClusters(S1_Td.subset, resolution = res)
S1_Td.subset <- RunUMAP(S1_Td.subset, dims = 1:nPC)

FeaturePlot(S1_Td.subset, c("Rgs5"), label = T)
FeaturePlot(S1_Td.subset, c("Pecam1"), label = T)

#keep only EC and Pericyte clusters
S1_Td_cc_EC_peri <- subset(S1_Td.subset, idents = c(1,3,4,5))


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

all.genesS2_Td.subset <- rownames(S2_Td.subset)

S2_Td.subset <- NormalizeData(S2_Td.subset, verbose = FALSE)
S2_Td.subset <- FindVariableFeatures(S2_Td.subset, selection.method = "vst", nfeatures = 2000)
S2_Td.subset <- CellCycleScoring(S2_Td.subset, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)

#ScaleDate and regress for MTpercent and nFeature_RNA
S2_Td.subset <- ScaleData(S2_Td.subset,
                          vars.to.regress = c("percent.mt", "nFeature_RNA", "G2M.Score", "S.Score"),
                          features = all.genesS2_Td.subset)



S2_Td.subset@commands$ScaleData.RNA

S2_Td.subset <- RunPCA(S2_Td.subset, npcs = 30, verbose = FALSE)
# Examine and visualize PCA results a few different ways
print(S2_Td.subset[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(S2_Td.subset, dims = 1:3, reduction = "pca")
DimPlot(S2_Td.subset, reduction = "pca")
DimHeatmap(S2_Td.subset, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(S2_Td.subset, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(S2_Td.subset)

nPC = 20
res = 0.1

S2_Td.subset <- FindNeighbors(S2_Td.subset, dims = 1:nPC)
S2_Td.subset <- FindClusters(S2_Td.subset, resolution = res)
S2_Td.subset <- RunUMAP(S2_Td.subset, dims = 1:nPC)
S2_Td.subset <- RunTSNE(S2_Td.subset, dims = 1:nPC)

FeaturePlot(S2_Td.subset, "Pecam1", order = T, label = T)
FeaturePlot(S2_Td.subset, "Rgs5", order = T, label = T)

#keep only EC and Pericyte clusters
S1_Td_cc_EC_peri = subset(S2_Td.subset, idents = c("0", "2"))

saveRDS(S1_Td_cc_EC_peri, file = "S2_Td_cc.rds")


########  INTEGRATED  ##########################################

#Integrate intact and injury samples and run the standard workflow for visualization and clustering


integrated_3TIP <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3TIP/integrated_EC_3TIP.Rds")

integrated_3TIP@commands$FindIntegrationAnchors

# normalize and find varible features in the objects
object_clean_new.list.bis <- lapply(X = c(S1_Td_cc_EC_peri, S2_Td_cc_EC_peri ), FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# find integration anchors
integrated1 <- FindIntegrationAnchors(object.list = object_clean_new.list.bis, dims = 1:30)

features.to.integrate1 = integrated1@anchor.features

integrated1 <- IntegrateData(anchorset = integrated1, dims = 1:30, features.to.integrate = features.to.integrate1)

DefaultAssay(integrated1) <- "integrated"


ElbowPlot(integrated1)

#res = 3
res = 1
# Run the standard workflow for visualization and clustering
integrated1 <- ScaleData(integrated, verbose = FALSE)
integrated1 <- RunPCA(integrated1, npcs = 35, verbose = FALSE)
# t-SNE and Clustering
ElbowPlot(integrated1)
integrated1 <- RunUMAP(integrated1, reduction = "pca", dims = 1:35)
integrated1 <- RunTSNE(integrated1, reduction = "pca", dims = 1:35)
integrated1 <- FindNeighbors(integrated1, reduction = "pca", dims = 1:35)
integrated1 <- FindClusters(integrated1, resolution = res)
DimPlot(integrated1, reduction = "umap", group.by = "stim")
DimPlot(integrated1, reduction = "umap", split.by = "stim")

DimPlot(integrated1, group.by = "stim")
DimPlot(integrated1)

table(integrated1$stim)
DefaultAssay(integrated1) <- "RNA"

FeaturePlot(integrated1, "Plvap", order = TRUE)



DefaultAssay(integrated_3TIP) = "integrated"

#integrated_3TIP <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3TIP/integrated_EC_3TIP.Rds")

# select EC and exclude pericytes from following analyses

subs = subset(integrated_3TIP, idents = c("0","1", "2", "3", "4", "5", "6","7","9","10"))

DimPlot(subs, split.by = "stim")

new.cluster.ids.lit <- c('VENOUS_PLVAP+',
                         'INTERMEDIATE',
                         'BARR_END_CAP',
                         'VENOUS_PLVAP+',
                         'ARTERIAL',
                         'CAPILLARY_PLVAP-',
                         'PROLIFERATING',
                         'VENOUS_PLVAP-',
                         'TIP',
                         'CAPILLARY_PLVAP+'
)


# refine cluster calling

DefaultAssay(subs) = "integrated"
res= 1
sub <- FindClusters(subs, resolution = res)
DimPlot(sub, split.by = "stim")

new.cluster.ids.lit <- c('VENOUS_PLVAP+',
                         'INTERMEDIATE',
                         'BARR_END_CAP',
                         'VENOUS_PLVAP+',
                         'ARTERIAL',
                         'CAPILLARY_PLVAP-',
                         'PROLIFERATING',
                         'VENOUS_PLVAP-',
                         'UNDETERMINED',
                         'TIP',
                         'CAPILLARY_PLVAP+'
)


names(new.cluster.ids.lit) <- levels(sub)
integrated_3TIP_new_undet <- RenameIdents(sub, new.cluster.ids.lit)

integrated_3TIP_new_undet$CellTypes = Idents(integrated_3TIP_new_undet)



# Fig1 UMAPs

#Fig1 E
pdf("EC_integrated_group_by_stim_cl8_no_legend.pdf", 12, 10)
DimPlot(integrated_3TIP_new_undet, group.by = "stim", pt.size = 5) + NoLegend()
dev.off()

#Fig1 D
pdf("EC_integrated_newcols_cl8_no_legend.pdf", 12, 10)
DimPlot(integrated_3TIP_new_undet, cols = c('VENOUS_PLVAP+' = '#FF6666',
                                            'VENOUS_PLVAP-' = '#990000',
                                            'BARR_END_CAP' = '#336666',
                                            'IMMATURE' = '#6600CC',
                                            'PROLIFERATING' = '#FF99CC',
                                            'TIP' = '#FF00FF',
                                            'INTERMEDIATE' = 'grey',
                                            'CAPILLARY_PLVAP-' = '#399933',
                                            'ARTERIAL' = '#0066FF',
                                            'CAPILLARY_PLVAP+' = '#99CC33'), pt.size = 5) + NoLegend()

dev.off()



# Fig1 FeaturePlots

#Fig1 F
p_list = FeaturePlot(object = integrated_3TIP_new_undet, features = c("Plvap", "Cldn5"),
                     cols= c("grey", "red", "blue"), blend = TRUE, order = TRUE, pt.size = 5, combine = FALSE)

p_list2<- append(p_list, list(legend = guide_area()), 4)

layout1<-"
ABCD
"


pdf("Plvap_Cldn5_explained.pdf", 40, 10)
wrap_plots(p_list , design = layout1)
dev.off()



# Fig1 heatmaps


angiogenesis = read.xlsx(
  "/Users/maurizio.aurora/Downloads/heatmap TIP cell db-2.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = T,
  rowNames = F,
)

annotation_row_angiogenesis = angiogenesis$gene
pathways_angiogenesis = rep("angiogenesis", 22)


TIP = read.xlsx(
  "/Users/maurizio.aurora/Downloads/heatmap TIP cell db-2.xlsx",
  sheet = 2,
  startRow = 1,
  colNames = T,
  rowNames = F,
)

annotation_row_TIP = TIP$gene
pathways_TIP = rep("TIP", 28)


immature = read.xlsx(
  "/Users/maurizio.aurora/Downloads/heatmap TIP cell db-2.xlsx",
  sheet = 3,
  startRow = 1,
  colNames = T,
  rowNames = F,
)

annotation_row_immature = immature$gene
pathways_immature = rep("immature", 26)


proliferating = read.xlsx(
  "/Users/maurizio.aurora/Downloads/heatmap TIP cell db-2.xlsx",
  sheet = 4,
  startRow = 1,
  colNames = T,
  rowNames = F,
)

annotation_row_proliferating = proliferating$gene
pathways_proliferating = rep("proliferating", 22)


tumor_angiogenesis = read.xlsx(
  "/Users/maurizio.aurora/Downloads/heatmap TIP cell db-2.xlsx",
  sheet = 5,
  startRow = 1,
  colNames = T,
  rowNames = F,
)

annotation_row_tumor_angiogenesis = tumor_angiogenesis$gene
pathways_tumor_angiogenesis  = rep("tumor_angiogenesis", 24)


oxidative_phosphororilation = read.xlsx(
  "/Users/maurizio.aurora/Downloads/heatmap TIP cell db-2.xlsx",
  sheet = 6,
  startRow = 1,
  colNames = T,
  rowNames = F,
)

annotation_row_oxidative_phosphororilation = oxidative_phosphororilation$gene
pathways_oxidative_phosphororilation = rep("oxidative_phosphororilation", 22)

annotation_row = c(annotation_row_proliferating,
                   annotation_row_immature,
                   annotation_row_TIP,
                   annotation_row_angiogenesis,
                   annotation_row_tumor_angiogenesis,
                   annotation_row_oxidative_phosphororilation)

pathways  = c(pathways_proliferating,
              pathways_immature,
              pathways_TIP,
              pathways_angiogenesis,
              pathways_tumor_angiogenesis,
              pathways_oxidative_phosphororilation)


DefaultAssay(integrated_3TIP_new_undet) = "RNA"

# remove INTERMEDIATE celltype

sub_obj <- subset(object = integrated_3TIP_new_undet, idents = "INTERMEDIATE", invert = TRUE)

samples <- sub_obj
avgexp = AverageExpression(samples, return.seurat = T)
counts <- GetAssayData(avgexp, assay="RNA", slot="data")
genes <- annotation_row
counts <- as.matrix(counts[rownames(counts) %in% genes, ])
cluster <- colnames(counts)
colnames(counts) = colnames(counts)
sample <- colnames(counts)


sample <- sample
cluster <- cluster
stim <- rep("integrated",9)

df <- data.frame(sample, cluster)
head(df)
df = df[order(df$sample),]


annotation_column = as.data.frame(df[,c('cluster')])
colnames(annotation_column)= "cluster"
rownames(annotation_column) <- df$cluster
head(annotation_column)


ann_colors = list(
  pathways = c("proliferating" = '#531554',
               "immature" = '#c83d81',
               "TIP" = '#D991EE',
               "angiogenesis" = '#FF8300',
               "tumor_angiogenesis" = 'red',
               "oxidative_phosphororilation" = 'yellow'),
  cluster = c(  'ARTERIAL' = '#0066FF',
                'BARR_END_CAP' = '#336666',
                'CAPILLARY_PLVAP-' = '#399933',
                'CAPILLARY_PLVAP+' = '#99CC33',
                'IMMATURE' = '#6600CC',
                'PROLIFERATING' = '#FF99CC',
                'TIP' = '#FF00FF',
                'VENOUS_PLVAP-' = '#990000',
                'VENOUS_PLVAP+' = '#FF6666'))


annotation_row = annotation_row
counts_ordered = counts[,row.names(annotation_column)]
head(counts_ordered)

sub_samp_ordered <- counts_ordered[annotation_row,,drop=FALSE]

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
minH = -1; maxH=1
myb = seq(minH, maxH, by = 0.01)
myc <- crp(length(myb))



length(pathways)
df_row <- data.frame(annotation_row, pathways)
row_annotation = df_row[2]
row.names(df_row) = df_row$annotation_row


draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

toMatch = c("Pgf","Adam12","Plxnd1","Col4a1","Apln","Cxcr4","Lcp2","Esm1","Dll4","Fscn1","Rps15a","Rps16","Rpl6",
            "Plvap","Hspe1","Ccna2", "Top2a","Aurkb","Cdk1","Chga","F2r","Cd93","Smad1","Cd276","Idh2","Acat1",
            "Acat1","Cox4i1","Ndufa8","Uqcrh")


matches <- grep(paste(toMatch,collapse="|"),
                df_row$annotation_row)

matched <- grep(paste(toMatch,collapse="|"),
                df_row$annotation_row, value = T)

ha = rowAnnotation(foo = anno_mark(at = matches, labels = matched))

pdf("heatmap_TIP_conditions_together_small_annotations.pdf", 20, 10)
ComplexHeatmap::pheatmap(sub_samp_ordered,
                         show_rownames = F,
                         show_colnames = F,
                         cellheight = 3,
                         cellwidth = 50,
                         cluster_rows = F,
                         cluster_cols = F,
                         #width = 5,
                         #height = 5,
                         annotation_col = annotation_column,
                         annotation_colors = ann_colors,
                         scale = 'row',
                         color = myc,
                         annotation_row = df_row[2],
                         fontsize_row = 1,
                         right_annotation = ha)
dev.off()


# Figure S1 heatmap


annotation_row = c(annotation_row_arterial,
                   annotation_row_barrier,
                   annotation_row_capillaryplvap_m,
                   annotation_row_capillaryplvap_p,
                   annotation_row_venous_plvap_m,
                   annotation_row_venous_plvap_p)

CellTypes  = c(pathways_arterial,
               pathways_barrier,
               pathways_capillary_plvap_m,
               pathways_capillary_plvap_p,
               pathways_venous_plvap_m,
               pathways_venous_plvap_p)



subs = subset(sub_obj, idents = c('VENOUS_PLVAP+',
                                  'VENOUS_PLVAP-',
                                  'BARR_END_CAP',
                                  'CAPILLARY_PLVAP-',
                                  'ARTERIAL',
                                  'CAPILLARY_PLVAP+'))


sample <- subs
avgexp = AverageExpression(sample, return.seurat = T)
counts <- GetAssayData(avgexp, assay="RNA", slot="data")
genes <- annotation_row
counts <- as.matrix(counts[rownames(counts) %in% genes, ])
cluster <- colnames(counts)
colnames(counts) = colnames(counts)
sample <- colnames(counts)

stim <- rep("integrated",6)

df <- data.frame(sample, cluster)
df = df[order(df$sample),]
annotation_column = as.data.frame(df[,c('cluster')])
colnames(annotation_column)= "cluster"
rownames(annotation_column) <- df$cluster

annotation_row = annotation_row
counts_ordered = counts[,row.names(annotation_column)]
counts_ordered = counts_ordered[annotation_row,]

head(counts_ordered)

sub_samp_ordered <- counts_ordered

sub_samp_ordered1 <- counts_ordered[rownames(df_row), ]

all_celltypes_but_tip = readLines("/Users/maurizio.aurora/Downloads/all_celltypes_but_tip.txt")

toMatch = all_celltypes_but_tip

matches <- grep(paste(toMatch,collapse="|"),
                df_row$annotation_row)

matched <- grep(paste(toMatch,collapse="|"),
                df_row$annotation_row, value = T)

ha = rowAnnotation(foo = anno_mark(at = matches, labels = matched))

# Fig S1 L

pdf('heatmap_conditions_together_small_all_ct_but_tip_small_annotation.pdf')
pheatmap(sub_samp_ordered1,
         show_rownames = F,
         show_colnames = F,
         cellheight = 2.5,
         cluster_rows = F,
         cluster_cols = F,
         annotation_col = annotation_column,
         annotation_colors = ann_colors,
         scale = 'row',
         color = myc,
         right_annotation = ha)
dev.off()

#Fig S3 H
pdf("Chga_FP_split.pdf", 10, 5)
FeaturePlot(sub_obj, "Chga", split.by = "stim",  pt.size =3, order = T)
dev.off()

#Fig S3 I
pdf("Cd276_FP_split.pdf", 10, 5)
FeaturePlot(sub_obj, "Cd276", split.by = "stim",  pt.size =5)
dev.off()

#KDR vlnplots figure S2 D

#sub_obj = readRDS("integrated_3TIP_new.RDS")

new.cluster.ids.lit <- c('EC',
                         'EC',
                         'EC',
                         'EC',
                         'EC',
                         'EC',
                         'EC',
                         'EC',
                         'EC'
)

names(new.cluster.ids.lit) <- levels(sub_obj)
temp <- RenameIdents(sub_obj, new.cluster.ids.lit)


DefaultAssay(temp) = "RNA"
pdf("VlnPlot_Kdr_intact_D7_nopoints.pdf", 12, 10)
VlnPlot(temp, "Kdr", split.by = "stim", cols = c("1" = '#4C9900',"2" = 'purple'), pt.size = 0)+
  theme(text = element_text(face = "bold"), title = element_text(size=25,face="bold"))+
  theme(legend.text = element_text(size=25))+
  theme(legend.key.height= unit(1, 'cm'),
        legend.key.width= unit(1, 'cm'))
dev.off()




