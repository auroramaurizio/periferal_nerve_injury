library("Seurat")

getwd()
#upload integrated object containing only EC subtypes
integrated_EC <-readRDS("integrated_EC.RDS")

#upload object with all cell types (including contaminant PERI and MES)
integrated = readRDS("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/Subclustering_mm10_M16_TdTomato_cc_2/seurat_objects/all_cell_types_scBonanomi_full.Rds")

#extract cell names for each ident in the EC dataset integrated_EC
ARTERIAL <-WhichCells(object=integrated_EC, idents="ARTERIAL")
IMMATURE <-WhichCells(object=integrated_EC, idents="IMMATURE")
PROLIFERATING <-WhichCells(object=integrated_EC, idents="PROLIFERATING")
TIP <-WhichCells(object=integrated_EC, idents="TIP")
INTERMEDIATE <-WhichCells(object=integrated_EC, idents="INTERMEDIATE")
BARR_END_CAP <-WhichCells(object=integrated_EC, idents="BARR_END_CAP")
VENOUS_PLVAP_M <-WhichCells(object=integrated_EC, idents="VENOUS_PLVAP-")
VENOUS_PLVAP_P <-WhichCells(object=integrated_EC, idents="VENOUS_PLVAP+")
CAPILLARY_PLVAP_M <-WhichCells(object=integrated_EC, idents="CAPILLARY_PLVAP-")
CAPILLARY_PLVAP_P <-WhichCells(object=integrated_EC, idents="CAPILLARY_PLVAP+")

#highligh specific cells in the full dataset integrated
DimPlot(integrated, label=T, cells.highlight=INTERMEDIATE, cols.highlight = c("cyan"), cols= "grey", order = T)


# Assign cell ID identified after EC cell type subset reclustering to the full object. 
# The full object contains also contaminant cell types such as Fibroblasts, Pericytes and LEC.
Idents(object = integrated, cells = TIP) <- 'TIP'
Idents(object = integrated, cells = PROLIFERATING) <- 'PROLIFERATING'
Idents(object = integrated, cells = IMMATURE) <- 'IMMATURE'
Idents(object = integrated, cells = IMMATURE) <- 'INTERMEDIATE'


DefaultAssay(integrated) = "RNA"

#isolate intact sample 
object_clean_new <- subset(x = integrated, subset = stim == "1")
DefaultAssay(object_clean_new) = "RNA"

#Isolate cluster marker genes
thresh.use = 0.1
min.pct = 0
min.diff.pct = -Inf
mg_exp.use ="wilcox"

cluster.marker_genes = FindAllMarkers(object_clean_new, logfc.threshold = thresh.use, mg_exp.use=mg_exp.use, min.pct=min.pct, min.diff.pct=min.diff.pct, only.pos=TRUE)

#quick check
FeaturePlot(object_clean_new, "Ldlr", label = T)
cluster.marker_genes[grep("Ldlr", cluster.marker_genes$gene), ]


filename_xls <- 'FindAllMarkers_min.pct_0.001.xlsx'
write.xlsx(cluster.marker_genes,
           file= filename_xls, 
           row.names = F,
           asTable = T)


mg_1 = cluster.marker_genes

# filter the marker genes for expression and adj pval
mg_pvalue_1 = mg_1[mg_1$p_val_adj < 0.05, ] 
mg_avg_logFC_1 = mg_pvalue_1[mg_pvalue_1$avg_log2FC > 0, ] 

# create a 2 column dataframe c("gene, "celltype)
mg_exp_bis <- as.data.frame(mg_avg_logFC_1 %>%
  group_by(gene) %>%
  summarise(temp = toString(cluster)) %>%
  ungroup())

# for each gene report the number of celltypes it is associated with
mg_exp_bis$counts <- lengths(strsplit(mg_exp_bis$temp, ","))
head(mg_exp_bis)


# for each gene report the number of celltypes it is associated with
mg_exp_1 <- mg_avg_logFC_1 %>%
  group_by(gene) %>%
  summarise(temp = toString(cluster)) %>%
  ungroup()

mg_exp_df_1 = as.data.frame(mg_exp_1)

# merge the 2 dataframes mg_exp_df_1 and mg_exp_bis so that in the final dataframe are reported
# for each gene, both the celltype showing the highest expression, both the list and number 
# of celltypes in which the gene is expressed
mg_exp_df_1 <- merge(mg_exp_df_1, mg_exp_bis, by = "gene", all.x = TRUE)
rownames(mg_exp_df_1) = mg_exp_df_1$gene
colnames(mg_exp_df_1) =c("gene_in_sc","celltypes_highest_LR", "celltypes","number_of_diff_celltypes")

# add column celltype_raw where, if multiple celltypes express the same gene, I report only one representative macro cluster.
# some celltypes have priority over others.
# e.g. VENOUS_PLVAP+  becomes EC; 
# VENOUS_PLVAP+ + LEC becomes EC;
# FIBRO + PERI becomes FIBRO

mg_exp_df = (mg_exp_df_1 %>% 
    mutate(celltype_raw = ifelse(celltypes_highest_LR %in% c("FIBRO", grep("FIBRO", celltypes_highest_LR, value=T)), "FIBRO",
                  ifelse(celltypes_highest_LR %in% c(grep("BARR_END_CAP", celltypes_highest_LR, value=T),"IMMATURE", "VENOUS_PLVAP-", "VENOUS_PLVAP+", "PROLIFERATING", "ARTERIAL", "LEC", "VENOUS", "CAPILLARY_PLVAP+", "CAPILLARY_PLVAP-"), "EC",
                   ifelse(celltypes_highest_LR %in% c(grep("IMMATURE", celltypes_highest_LR, value=T)), "EC",    
                  ifelse(celltypes_highest_LR %in% c(grep("NK", celltypes_highest_LR, value=T)), "NK",
                  ifelse(celltypes_highest_LR %in% c(grep("IMMATURE", celltypes_highest_LR, value=T)), "EC",          
                  ifelse(celltypes_highest_LR %in% c(grep("PERI", celltypes_highest_LR, value=T)), "PERICYTES", NA))))))))

# add column celltype_raw where, if multiple celltypes express the same gene, I report only one representative macro cluster.
# an option of the above one
# if one gene is expressed in more than one celltype report MIX

mg_exp_df = (mg_exp_df_1 %>% 
    mutate(celltype_raw = ifelse(number_of_diff_celltypes > 1 , "MIX",
                  ifelse(celltypes_highest_LR %in% c(grep("BARR_END_CAP", celltypes_highest_LR, value=T),"IMMATURE", "VENOUS_PLVAP-", "VENOUS_PLVAP+", "PROLIFERATING", "ARTERIAL", "LEC", "VENOUS", "CAPILLARY_PLVAP+", "CAPILLARY_PLVAP-"), "EC",
                   ifelse(celltypes %in% c(grep("IMMATURE", celltypes_highest_LR, value=T)), "EC",
                  ifelse(celltypes %in% c(grep("FIBRO", celltypes_highest_LR, value=T)), "FIBRO", 
                  ifelse(celltypes %in% c(grep("NK", celltypes_highest_LR, value=T)), "NK",
                  ifelse(celltypes %in% c(grep("IMMATURE", celltypes_highest_LR, value=T)), "EC",  
                  ifelse(celltypes %in% c(grep("TIP", celltypes_highest_LR, value=T)), "EC",  
                  ifelse(celltypes %in% c(grep("PERI", celltypes_highest_LR, value=T)), "PERICYTES", NA))))))))))

tail(mg_exp_df, 40)

# check some genes 
mg_exp_df[grep("Dsn1", mg_exp_df$gene_in_sc), ]

###################

# Upload the Cherry samples bulk RNAseq dataset (DGE results table ) 

DGE_cherry = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/only_KDRCherry_Combat_no_outliers/DGE_results_Cherry/DGE_results_Cherry.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)


# merge bulk DGE results table with teh table created starting from sc data 
cherry <- merge(DGE_cherry, mg_exp_df, by = "row.names", all.x = TRUE)
head(mg_exp_df)

# sort by pvalue
cherry_sort <- cherry[order(cherry$pvalue, decreasing = F),]  

SAVE_variable <- list()
filename_xls <- "DGE_results_Cherry_celltype_0.001.xlsx"
variable2save_names <- 'KO_vs_WT'


# all_counts

SAVE_variable[[variable2save_names[1]]] <- cherry_sort

# expGENES_counts


write.xlsx(SAVE_variable,
           file = "DGE_results_Cherry_celltype_0.001.xlsx", 
           row.names = F,
           asTable = F, 
           sheetName =variable2save_names)



###############################

# Upload the Tomato samples bulk RNAseq dataset (DGE results table ) 

DGE_tomato_14_vs_7 = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)



DGE_tomato_14_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 2,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)


DGE_tomato_7_vs_intact = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/DGE_results_Tomato/DGE_results_Tomato.xlsx",
  sheet = 3,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)


merged_1 <- merge(DGE_tomato_14_vs_7, mg_exp_df, by = "row.names", all.x = TRUE)
merged_1_sort <- merged_1[order(merged_1$pvalue, decreasing = F),]  



merged_2 <- merge(DGE_tomato_14_vs_intact, mg_exp_df, by = "row.names", all.x = TRUE)
merged_2_sort <- merged_2[order(merged_2$pvalue, decreasing = F),]  


merged_3 <- merge(DGE_tomato_7_vs_intact, mg_exp_df, by = "row.names", all.x = TRUE)
merged_3_sort <- merged_3[order(merged_3$pvalue, decreasing = F),]  


SAVE_variable <- list()
filename_xls <- "DGE_results_Tomato_celltype_0.001.xlsx"
variable2save_names <- c('crush_d14_vs_crush_d7', 'crush_d14_vs_intact','crush_d7_vs_intact')


SAVE_variable[[variable2save_names[1]]] <- merged_1_sort

SAVE_variable[[variable2save_names[2]]] <- merged_2_sort

SAVE_variable[[variable2save_names[3]]] <- merged_3_sort

write.xlsx(SAVE_variable,
           file = filename_xls, 
           row.names = F,
           asTable = F, 
           sheetName =variable2save_names)

