#install.packages("devtools")
#library(devtools)
#install_github("neurogenomics/EWCE")

###################################################################################
#for examples see the vignette https://nathanskene.github.io/EWCE/articles/EWCE.html#application-to-transcriptomic-data-1
###################################################################################

library(Seurat)
library(devtools)
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)


#the main function we require is
#help(generate.celltype.data)


#it takes as input
#exp
#annotLevels

#generate exp from my scdataset
#Numerical matrix with row for each gene and column for each cell.
#Row names are MGI/HGNC gene symbols.
#Column names are cell IDs which can be cross referenced against the annot data frame.

integrated_3TIP <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3TIP/integrated_EC_Bonanomi_newnames.RDS")
integrated_3TIP <- subset( integrated_3TIP, idents = c("INTERMEDIATE"), invert = TRUE)
DimPlot(integrated_3TIP)
DefaultAssay(integrated_3TIP) = "RNA"

#https://www.cellphonedb.org/faq-and-troubleshooting

# take raw data. In seurat 3 raw.data are called counts
count_raw =as.matrix(GetAssayData(object = integrated_3TIP, slot = "counts"))

# think about performing some kind of normalization

# gene names are stored as rownames
#rownames(count_raw)
# cell names are stored as colnames
#colnames(count_raw)

#A different number of reads is found across each cell.
#We suggest using scTransform to normalise for differences due to cell size,
#then linearly scale. Note that this might be slow.

# devtools::install_github(repo = 'ChristophH/sctransform')

#Note that this is optional (and was not used in the original EWCE publication)
#so by all means ignore this scTransform step.

#library(sctransform)
#scT = sctransform::vst(count_raw, return_cell_attr = TRUE)

#count_raw$exp_scT = correct_counts(scT, count_raw) # umi_corrected

#count_raw$exp_scT_normed = Matrix::t(Matrix::t(count_raw$exp_scT)*(1/Matrix::colSums(count_raw$exp_scT)))

#generate annotLevels
#List with arrays of strings containing the cell type names associated with each column in exp

# create meta data object
meta_data <- cbind(rownames(integrated_3TIP@meta.data),as.data.frame(Idents(integrated_3TIP)))
# create meta data "level2class" definition. It is more precise with EC subset definition
colnames(meta_data)<- c("cell_id","level1class") # was written level2class argh
# create "level1class". all EC subtypes will be called EC
meta_data$level1class <- gsub("IMMATURE", "IMMATURE", meta_data$level1class) # was written level2class argh
meta_data$level1class <- gsub("PROLIFERATING", "PROLIFERATING", meta_data$level1class)
meta_data$level1class <- gsub("TIP", "TIP", meta_data$level1class)
meta_data$level1class <- gsub("VENOUS_PLVAP\\-", "VENOUS_PLVAP\\-", meta_data$level1class)
meta_data$level1class <- gsub("VENOUS_PLVAP\\+", "VENOUS_PLVAP\\+", meta_data$level1class)
meta_data$level1class <- gsub("BARR_END_CAP", "BARR_END_CAP", meta_data$level1class)
meta_data$level1class <- gsub("ARTERIAL", "ARTERIAL", meta_data$level1class)
meta_data$level1class <- gsub("CAPILLARY_PLVAP\\-", "CAPILLARY_PLVAP\\-", as.character(meta_data$level1class))
meta_data$level1class <- gsub("CAPILLARY_PLVAP\\+", "CAPILLARY_PLVAP\\+", meta_data$level1class)
#check if everything went well
unique(meta_data$level1class)
level1class = meta_data$level1class

level2class = meta_data$level1class


#1) Drop genes which do not show significant evidence of varying between level 
exp_DROPPED = drop.uninformative.genes(exp=count_raw,level2annot = level2class)
#2) Calculate cell type averages and specificity for each gene 
annotLevels = list(level1class=level1class,level2class=level2class)
#3) Drop all genes which do not have 1:1 mouse:human orthologs
fNames = filter.genes.without.1to1.homolog(fNames)
fNames = generate.celltype.data(exp=exp_DROPPED,annotLevels=annotLevels,groupName="EWCE_ndn")

#integrated_3TIP <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3TIP/integrated_EC_3TIP.Rds")
#integrated = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/Subclustering_mm10_M16_TdTomato_cc_2/seurat_objects/all_cell_types_scBonanomi_full.Rds")

#load(fNames[2])

# use the List generated using generate.celltype.data and a differential expression
# table and determines the probability of cell-type enrichment in the up & down regulated genes
#help(ewce_expression_data)

#load dge tables
crush_D14_vs_crush_D7 = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

crush_D14_vs_crush_D7$MGI.symbol = rownames(crush_D14_vs_crush_D7)



crush_D14_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 2,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

crush_D14_vs_intact$MGI.symbol = rownames(crush_D14_vs_intact)



crush_D7_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 3,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)

crush_D7_vs_intact$MGI.symbol = rownames(crush_D7_vs_intact)

load(file = "~//CellTypeData_EWCE_ndn.rda") #now everything is in ctd

#load(file = "CellTypeData_DB.rda") #now everything is in ctd

setwd("/Users/maurizio.aurora/Documents/GitHub/MyOwnGithub/nerve_injury/bulkRNAseq")
#ctd[[1]]$plotting
#ctd[[2]]$plotting

############################################################################################
############################################################################################
############################################################################################


FDR_crush_D7_vs_intact = crush_D7_vs_intact[crush_D7_vs_intact$padj < 0.05, ]
FDR_tt_results_crush_D7_vs_intact = ewce_expression_data(sct_data=ctd,tt=FDR_crush_D7_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                         reps = 10000)


write.xlsx(FDR_tt_results_crush_D7_vs_intact$joint_results,
           file = "FDR_D7_vs_intact_barplot.xlsx",
           row.names = T,
           asTable = T)

FDR_tt_results_crush_D7_vs_intact$joint_results$SD_from_mean = FDR_tt_results_crush_D7_vs_intact$joint_results$sd_from_mean
FDR_tt_results_crush_D7_vs_intact$joint_results$SD_from_mean [which(FDR_tt_results_crush_D7_vs_intact$joint_results$SD_from_mean < 0)] = 0
FDR_tt_results_crush_D7_vs_intact$joint_results$pvalue = NA
FDR_tt_results_crush_D7_vs_intact$joint_results$pvalue[FDR_tt_results_crush_D7_vs_intact$joint_results$p<0.05]<-'*'


S3_test_D7 <- ggplot(FDR_tt_results_crush_D7_vs_intact$joint_results)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

pdf("FDR_crush_D7_vs_intact_dotplot_minSDzero_pvalue0.05_t.pdf",8,10)
S3_test_D7
dev.off()

############################################################################################
############################################################################################
############################################################################################

ctd[[1]]


FDR_crush_D14_vs_intact = crush_D14_vs_intact[crush_D14_vs_intact$padj < 0.05, ]
FDR_tt_results_crush_D14_vs_intact = ewce_expression_data(sct_data=ctd,tt=FDR_crush_D14_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                          reps = 10000)


write.xlsx(FDR_tt_results_crush_D14_vs_intact$joint_results,
           file = "FDR_D14_vs_intact_barplot.xlsx",
           row.names = T,
           asTable = T)

FDR_tt_results_crush_D14_vs_intact$joint_results$SD_from_mean = FDR_tt_results_crush_D14_vs_intact$joint_results$sd_from_mean
FDR_tt_results_crush_D14_vs_intact$joint_results$SD_from_mean [which(FDR_tt_results_crush_D14_vs_intact$joint_results$SD_from_mean < 0)] = 0
FDR_tt_results_crush_D14_vs_intact$joint_results$pvalue = NA
FDR_tt_results_crush_D14_vs_intact$joint_results$pvalue[FDR_tt_results_crush_D14_vs_intact$joint_results$p<0.05]<-'*'


S3_test_D14 <- ggplot(FDR_tt_results_crush_D14_vs_intact$joint_results)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

pdf("FDR_crush_D14_vs_intact_dotplot_minSDzero_pvalue0.05_t.pdf",8,10)
S3_test_D14
dev.off()


############################################################################################
############################################################################################
############################################################################################



FDR_crush_D14_vs_crush_D7 = crush_D14_vs_crush_D7[crush_D14_vs_crush_D7$padj < 0.05, ]
FDR_tt_results_crush_D14_vs_crush_D7 = ewce_expression_data(sct_data=ctd,tt=FDR_crush_D14_vs_crush_D7,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                            reps = 10000)

write.xlsx(FDR_tt_results_crush_D14_vs_crush_D7$joint_results,
           file = "FDR_D14_vs_D7_barplot.xlsx",
           row.names = T,
           asTable = T)

FDR_tt_results_crush_D14_vs_crush_D7$joint_results$SD_from_mean = FDR_tt_results_crush_D14_vs_crush_D7$joint_results$sd_from_mean
FDR_tt_results_crush_D14_vs_crush_D7$joint_results$SD_from_mean [which(FDR_tt_results_crush_D14_vs_crush_D7$joint_results$SD_from_mean < 0)] = 0
FDR_tt_results_crush_D14_vs_crush_D7$joint_results$pvalue = NA
FDR_tt_results_crush_D14_vs_crush_D7$joint_results$pvalue[FDR_tt_results_crush_D14_vs_crush_D7$joint_results$p<0.05]<-'*'


S3_test_D14_D7 <- ggplot(FDR_tt_results_crush_D14_vs_crush_D7$joint_results)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

pdf("FDR_crush_D14_vs_D7_dotplot_minSDzero_pvalue0.05_t.pdf",8,10)
S3_test_D14_D7
dev.off()

############################################################################################
############################################################################################
############################################################################################



d14d7 = FDR_tt_results_crush_D14_vs_crush_D7$joint_results
d14intact = FDR_tt_results_crush_D14_vs_intact$joint_results
d7intact = FDR_tt_results_crush_D7_vs_intact$joint_results

d14d7$cond = "crushD14_vs_crushD7"
d14intact$cond = "crush_D14_vs_intact"
d7intact$cond = "crush_D7_vs_intact"

combined = rbind(d7intact, d14intact, d14d7 )

combined$cond_f <- factor(combined$cond, levels=c("crush_D7_vs_intact", "crush_D14_vs_intact", "crushD14_vs_crushD7"))


combined_plot <- ggplot(combined)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))+
  theme(axis.text.x = element_text(colour = c("blue", "red")))



pdf("EWCE_tomato_3.pdf",13,8)
combined_plot + facet_grid(. ~ cond_f)
dev.off()


############################################################################################
############################################################################################
############################################################################################


combined = rbind(d7intact, d14d7 )

combined$cond_f <- factor(combined$cond, levels=c("crush_D7_vs_intact","crushD14_vs_crushD7"))


combined_plot <- ggplot(combined)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))+
  theme(axis.text.x = element_text(colour = c("blue", "red")))


pdf("EWCE_tomato_2.pdf",10,8)
combined_plot + facet_grid(. ~ cond_f)
dev.off()

############################################################################################
############################################################################################
############################################################################################

packageVersion("Seurat")




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

TU = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Toma/reanalysis/Seurat_objects/EC_Toma_Uninjured.Rds")
TI = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Toma/reanalysis/Seurat_objects/EC_Toma_3D.Rds")
KI = readRDS("/Users/maurizio.aurora/integrated_EC_Kalinski_no_pericytes.Rds")
CI = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/fibroblasts_Carr/reanalysis/Seurat_objects/EC_Carr_9D.Rds")

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

DimPlot(integrated, label=T, cells.highlight=TIP_3, cols.highlight = c("cyan"), cols= "grey", reduction = "umap")

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

257 - 246

old_names = factor(integrated_rn$stim)
levels(old_names) = c("7D","Intact","9D","3D","3D","Intact")
new.names = as.character(old_names)
integrated_rn$condition<- new.names

DimPlot(integrated_rn)
DimPlot(integrated_rn, label=T, cells.highlight=TIP_1, cols.highlight = c("cyan"), cols= "grey", reduction = "umap")



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

pdf("prippp.pdf", 15, 10)
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

