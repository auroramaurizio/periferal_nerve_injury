#install.packages("devtools")
#library(devtools)
# install_github("neurogenomics/EWCE")

###################################################################################
# for examples see the vignette 
# https://nathanskene.github.io/EWCE/articles/EWCE.html#application-to-transcriptomic-data-1
###################################################################################

library(Seurat)
library(devtools)
library(EWCE)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)


# the main function we require is
# help(c)


# it takes as input
# exp
# annotLevels

# generate exp from my scRNAseq dataset, only EC, remove INTERMEDIATE cluster
# Numerical matrix with row for each gene and column for each cell.
# Row names are gene symbols.
# Column names are cell IDs which can be cross referenced against the annot data frame.

integrated_3TIP <-readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/EC_subset/3TIP/integrated_EC_Bonanomi_newnames.RDS")
integrated_3TIP <- subset( integrated_3TIP, idents = c("INTERMEDIATE"), invert = TRUE)
DefaultAssay(integrated_3TIP) = "RNA"

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
fNames = generate.celltype.data(exp=exp_DROPPED,annotLevels=annotLevels,groupName="EWCE_ndn")
fNames = filter.genes.without.1to1.homolog(fNames)


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


KO_vs_WT = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/only_KDRCherry_Combat_no_outliers/DGE_results_Cherry/DGE_results_Cherry.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)


KO_vs_WT$MGI.symbol = rownames(KO_vs_WT)

#load(file = "CellTypeData_DB.rda") #now everything is in ctd

#setwd("/Users/maurizio.aurora/Documents/GitHub/MyOwnGithub/nerve_injury/bulkRNAseq")
#ctd[[1]]$plotting
#ctd[[2]]$plotting

############################################################################################
############################################################################################
############################################################################################


crush_D7_vs_intact = crush_D7_vs_intact[crush_D7_vs_intact$padj < 0.05, ]
tt_results_crush_D7_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D7_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                         reps = 10000)


write.xlsx(tt_results_crush_D7_vs_intact$joint_results,
           file = "D7_vs_intact_barplot.xlsx",
           row.names = T,
           asTable = T)

tt_results_crush_D7_vs_intact$joint_results$SD_from_mean = tt_results_crush_D7_vs_intact$joint_results$sd_from_mean
tt_results_crush_D7_vs_intact$joint_results$SD_from_mean [which(tt_results_crush_D7_vs_intact$joint_results$SD_from_mean < 0)] = 0
tt_results_crush_D7_vs_intact$joint_results$pvalue = NA
tt_results_crush_D7_vs_intact$joint_results$pvalue[tt_results_crush_D7_vs_intact$joint_results$p<0.05]<-'*'


S3_test_D7 <- ggplot(tt_results_crush_D7_vs_intact$joint_results)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

pdf("crush_D7_vs_intact_dotplot_minSDzero_pvalue0.05_t.pdf",8,10)
S3_test_D7
dev.off()

############################################################################################
############################################################################################
############################################################################################



crush_D14_vs_intact = crush_D14_vs_intact[crush_D14_vs_intact$padj < 0.05, ]
tt_results_crush_D14_vs_intact = ewce_expression_data(sct_data=ctd,tt=crush_D14_vs_intact,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                          reps = 10000)


write.xlsx(tt_results_crush_D14_vs_intact$joint_results,
           file = "D14_vs_intact_barplot.xlsx",
           row.names = T,
           asTable = T)

tt_results_crush_D14_vs_intact$joint_results$SD_from_mean = tt_results_crush_D14_vs_intact$joint_results$sd_from_mean
tt_results_crush_D14_vs_intact$joint_results$SD_from_mean [which(tt_results_crush_D14_vs_intact$joint_results$SD_from_mean < 0)] = 0
tt_results_crush_D14_vs_intact$joint_results$pvalue = NA
tt_results_crush_D14_vs_intact$joint_results$pvalue[tt_results_crush_D14_vs_intact$joint_results$p<0.05]<-'*'


S3_test_D14 <- ggplot(tt_results_crush_D14_vs_intact$joint_results)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

pdf("crush_D14_vs_intact_dotplot_minSDzero_pvalue0.05_t.pdf",8,10)
S3_test_D14
dev.off()


############################################################################################
############################################################################################
############################################################################################



crush_D14_vs_crush_D7 = crush_D14_vs_crush_D7[crush_D14_vs_crush_D7$padj < 0.05, ]
tt_results_crush_D14_vs_crush_D7 = ewce_expression_data(sct_data=ctd,tt=crush_D14_vs_crush_D7,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                                            reps = 10000)

write.xlsx(tt_results_crush_D14_vs_crush_D7$joint_results,
           file = "D14_vs_D7_barplot.xlsx",
           row.names = T,
           asTable = T)

tt_results_crush_D14_vs_crush_D7$joint_results$SD_from_mean = tt_results_crush_D14_vs_crush_D7$joint_results$sd_from_mean
tt_results_crush_D14_vs_crush_D7$joint_results$SD_from_mean [which(tt_results_crush_D14_vs_crush_D7$joint_results$SD_from_mean < 0)] = 0
tt_results_crush_D14_vs_crush_D7$joint_results$pvalue = NA
tt_results_crush_D14_vs_crush_D7$joint_results$pvalue[tt_results_crush_D14_vs_crush_D7$joint_results$p<0.05]<-'*'


S3_test_D14_D7 <- ggplot(tt_results_crush_D14_vs_crush_D7$joint_results)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

pdf("crush_D14_vs_D7_dotplot_minSDzero_pvalue0.05_t.pdf",8,10)
S3_test_D14_D7
dev.off()

############################################################################################
############################################################################################
############################################################################################

KO_vs_WT = KO_vs_WT[KO_vs_WT$padj < 0.05, ]
tt_results_KO_vs_WT = ewce_expression_data(sct_data=ctd,tt=KO_vs_WT,annotLevel=1,ttSpecies="mouse",sctSpecies="mouse",sortBy="log2FoldChange",
                                               reps = 10000)

write.xlsx(tt_results_KO_vs_WT$joint_results,
           file = "KO_vs_WT_barplot.xlsx",
           row.names = T,
           asTable = T)


tt_results_KO_vs_WT$joint_results$SD_from_mean = tt_results_KO_vs_WT$joint_results$sd_from_mean
tt_results_KO_vs_WT$joint_results$SD_from_mean [which(tt_results_KO_vs_WT$joint_results$SD_from_mean < 0)] = 0
tt_results_KO_vs_WT$joint_results$pvalue = NA
tt_results_KO_vs_WT$joint_results$pvalue[tt_results_KO_vs_WT$joint_results$p<0.05]<-'*'

S3_test <- ggplot(tt_results_KO_vs_WT$joint_results)+
  geom_point(aes(x=Direction, y=CellType, color=Direction,size = SD_from_mean))+
  scale_color_manual(values = c("Up" = 'red','Down' = 'blue'))+
  theme_classic() + scale_x_discrete(position = "top") +
  geom_text(aes(x=Direction, y=CellType, label = pvalue),size = 10,  color = "black") +
  scale_size(range = c(5, 20))+
  theme(text = element_text(size=20, face = "bold"))

# Fig3 E

pdf("crush_KO_vs_WT_dotplot_minSDzero_pvalue0.05.pdf",8,10)
S3_test
dev.off()
###########

############################################################################################
############################################################################################
############################################################################################



d14d7 = tt_results_crush_D14_vs_crush_D7$joint_results
d14intact = tt_results_crush_D14_vs_intact$joint_results
d7intact = tt_results_crush_D7_vs_intact$joint_results

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


# Fig 2C
pdf("EWCE_tomato_2.pdf",10,8)
combined_plot + facet_grid(. ~ cond_f)
dev.off()

############################################################################################
############################################################################################
############################################################################################

