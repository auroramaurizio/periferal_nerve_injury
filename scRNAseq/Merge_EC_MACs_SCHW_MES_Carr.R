library(Seurat)
library(ggplot2)
packageVersion("Seurat")

# The signaling crosstalk between ECs and other cellTypes involved in nerve injury was investigated by merging
# ECs from our intact/7 dpi scRNA-seq datasets, together with Schwann cells, macrophages, and Mesenchymal cells 
# from other studies {Toma, 2020; Ydens, 2020; Carr, 2019 }.
# In parallel the same analysis was performed with MES cells from our study instead of Carr's.
# In this script we use Carr's MES.


# Merge intact datasets

merged_int <- merge(EC_Intact, y = c(MES_Intact,MACROPHAGES_Intact, SCHWANN_Intact),
                    add.cell.ids = c("EC", "MES", "MACROPHAGES", "SCHWANN"), project = "Intact")



merged_int[["percent.mt"]] <- PercentageFeatureSet(merged_int, pattern = "^mt-")
merged_int[["percent.Rpl"]] <- PercentageFeatureSet(merged_int, pattern = "Rpl")
merged_int <- NormalizeData(merged_int, verbose = FALSE)
merged_int <- FindVariableFeatures(merged_int, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
merged_int <- ScaleData(merged_int, vars.to.regress = c("percent.mt", "nFeature_RNA"))
merged_int <- RunPCA(merged_int, npcs = 20, verbose = FALSE)

nPC = 20
res = 0.5

merged_int <- FindNeighbors(merged_int, dims = 1:nPC)
merged_int <- RunUMAP(merged_int, dims = 1:nPC)
merged_int <- RunTSNE(merged_int, dims = 1:nPC)

DimPlot(merged_int, label = T, repel = T)
DefaultAssay(merged_int) = "RNA"


getwd()

#saveRDS(merged_int, "Intact_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")
#merged_int = readRDS("Intact_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")


# Merge injury datasets


merged_inj <- merge(EC_Injury, y = c(MES_Injury,MACROPHAGES_Injury, SCHWANN_Injury),
                add.cell.ids = c("EC", "MES", "MACROPHAGES", "SCHWANN"), project = "Injury")

merged_inj[["percent.mt"]] <- PercentageFeatureSet(merged_inj, pattern = "^mt-")
merged_inj[["percent.Rpl"]] <- PercentageFeatureSet(merged_inj, pattern = "Rpl")
merged_inj <- NormalizeData(merged_inj, verbose = FALSE)
merged_inj <- FindVariableFeatures(merged_inj, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
merged_inj <- ScaleData(merged_inj, vars.to.regress = c("percent.mt", "nFeature_RNA"))
merged_inj <- RunPCA(merged_inj, npcs = 20, verbose = FALSE)

nPC = 20
res = 0.5

merged_inj <- FindNeighbors(merged_inj, dims = 1:nPC)
merged_inj <- RunUMAP(merged_inj, dims = 1:nPC)
merged_inj <- RunTSNE(merged_inj, dims = 1:nPC)

DimPlot(merged_inj, label = T, repel = T)
DefaultAssay(merged_inj) = "RNA"

#saveRDS(merged_inj, "Injury_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")
#merged_inj = readRDS("Injury_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")

# Intact + Injury datasets

DefaultAssay(merged_inj) = "RNA"
DefaultAssay(merged_int) = "RNA"

EC_Intact <- subset(x = merged_int, subset = stim == "1")
EC_Injury <- subset(x = merged_inj, subset = stim == "2")  
MACROPHAGES_Intact <- subset(x = merged_int, subset = stim == "macs_Intact") 
MACROPHAGES_Injury <- subset(x = merged_inj, subset = stim == "macs_Injury")
SCHWANN_Intact <- subset(x = merged_int, subset = stim == "schwann_Intact") 
SCHWANN_Injury <- subset(x = merged_inj, subset = stim == "schwann_Injury") 
MES_Intact <- subset(x = merged_int, subset = stim == "mes_Intact")   
MES_Injury <- subset(x = merged_inj, subset = stim == "mes_Injury")  

# Merge intact and injury datasets

merged <- merge(EC_Intact, y = c(MES_Intact,MACROPHAGES_Intact, SCHWANN_Intact, EC_Injury,MES_Injury,MACROPHAGES_Injury, SCHWANN_Injury), 
                add.cell.ids = c("EC_Intact", "MES_Intact", "MACROPHAGES_Intact", "SCHWANN_Intact", "EC_Injury", "MES_Injury", "MACROPHAGES_Injury", "SCHWANN_Injury"), project = "IntactInjury")


merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^mt-")
merged[["percent.Rpl"]] <- PercentageFeatureSet(merged, pattern = "Rpl")
merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
#ScaleDate and regress for MTpercent and nFeature_RNA
merged <- ScaleData(merged, vars.to.regress = c("percent.mt", "nFeature_RNA"))
merged <- RunPCA(merged, npcs = 20, verbose = FALSE)

nPC = 20 
res = 0.5

merged <- FindNeighbors(merged, dims = 1:nPC)
merged <- RunUMAP(merged, dims = 1:nPC)
merged <- RunTSNE(merged, dims = 1:nPC)

#saveRDS(merged, merged_all_cell_types.Rds)

#merged = readRDS("merged_all_cell_types.Rds")

merged@meta.data$stim[merged@meta.data$stim == "1"] <- "EC_Intact"
merged@meta.data$stim[merged@meta.data$stim == "2"] <- "EC_Injury"

# Semaphorin DotPlot

# upload semaphorin genes to plot
sema = readLines("/Users/maurizio.aurora/Downloads/sema.txt")


merged$stim <- factor(merged$stim, 
                                  levels=c("EC_Intact", "EC_Injury",
                                            "mes_Intact", "mes_Injury", 
                                            "macs_Intact", "macs_Injury",
                                            "schwann_Intact", "schwann_Injury"))



DefaultAssay(merged) = "RNA"
table(Idents(merged))
new.cluster.ids.lit <- c('EC', 
                         'EC',
                         'EC',
                         'EC',
                         'MACS',
                         'MES_DIFF',
                         'MES_ENDO',
                         'MES_EPI',
                         'MES_PERI',
                         'SCHWANN',
                         'EC',
                         'EC',
                         'EC',
                         'EC',
                         'EC'
)


names(new.cluster.ids.lit) <- levels(merged)
object_new <- RenameIdents(merged, new.cluster.ids.lit)



test = c("EC", "MES_DIFF",
         "MES_ENDO", "MES_EPI", 
         "MES_PERI", "MACS",
         "SCHWANN")

levels(object_new) <- test

object_new$dataset <- factor(object_new$stim, 
                                     levels=c("EC_Intact", "EC_Injury",
                                              "mes_Intact", "mes_Injury", 
                                              "macs_Intact", "macs_Injury",
                                              "schwann_Intact", "schwann_Injury"))


# Fig S3 Semaphorin dotplot
DefaultAssay(object_new) = "RNA"
pdf("Semadotplot_splitbystim_newcols_green_violet_test.pdf", 8, 7)
DotPlot(object_new, 
        split.by = "dataset",
        features = sema, 
        dot.scale = 8,
        cols = c("limegreen","darkviolet","limegreen","darkviolet","limegreen","darkviolet","limegreen","darkviolet"), 
        assay = "RNA") + 
  coord_flip() +
  RotatedAxis() 
dev.off()

#KDR dotplot for referees

pdf("DotPlot_Kdr_splitbystim_newcols_green_violet.pdf", 8, 3)
DotPlot(object_new, 
        split.by = "dataset",
        features = "Kdr", 
        dot.scale = 8,
        cols = c("limegreen","darkviolet","limegreen","darkviolet","limegreen","darkviolet","limegreen","darkviolet"), 
        assay = "RNA") + 
  coord_flip() +
  RotatedAxis() 
dev.off()


#KDR vlnplots for referees

pdf("VlnPlot_Kdr_splitbystim_newcols_green_violet_median_line.pdf", 8, 4)
VlnPlot(object_new, "Kdr", group.by = "dataset", pt.size = 0, cols = c("limegreen","darkviolet","limegreen","darkviolet","limegreen","darkviolet","limegreen","darkviolet")) +
  stat_summary(fun.y=median, geom="point", color="black", shape = 95,  size = 12)+
  theme(legend.position="none")
dev.off()


pdf("VlnPlot_Kdr_median_line.pdf", 8, 4)
VlnPlot(object_new, "Kdr",  pt.size = 0) +
  stat_summary(fun.y=median, geom="point", color="black", shape = 95,  size = 15)+
  theme(legend.position="none")
dev.off()


pdf("DotPlot_Kdr.pdf", 8, 4)
DotPlot(object_new, 
        features = "Kdr", 
        dot.scale = 8,
        assay = "RNA") + 
  coord_flip() +
  RotatedAxis() 
dev.off()


