library(monocle) # OK
suppressPackageStartupMessages(library(SeuratWrappers)) # NO
suppressPackageStartupMessages(library(Seurat)) # OK
#suppressPackageStartupMessages(library(SeuratData)) # NO
suppressPackageStartupMessages(library(ggplot2)) # OK
suppressPackageStartupMessages(library(patchwork)) # OK
suppressPackageStartupMessages(library(magrittr)) # OK
suppressPackageStartupMessages(library(future)) # OK
suppressPackageStartupMessages(library(cowplot)) # OK
suppressPackageStartupMessages(library(dplyr)) # OK

###################################################################
getwd()


integrated = readRDS(file ="/beegfs/scratch/ric.cosr/ric.brendolan/BrendolanA_1280_scRNAseq/dataset/090521/integrated_without_cluster8.rds")

#DimPlot(integrated, split.by = "stim")

#DefaultAssay(integrated) = "integrated"

#integrated$CellType = Idents(integrated)

#########################################
#########################################
#########################################

var_genes <- integrated[["integrated"]]@var.features

#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(integrated@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = integrated@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
i_monocle2 <- newCellDataSet(data,
                             phenoData = pd,
                             featureData = fd,
                             #lowerDetectionLimit = 0.5,
                             expressionFamily = uninormal())# since I have already normalized, thresholded and scaled data

var_genes <- integrated[["integrated"]]@var.features
ordering_genes <- var_genes

i_monocle2 <- setOrderingFilter(i_monocle2, ordering_genes)

i_monocle2 <- reduceDimension(i_monocle2, max_components = 2,
                              method = 'DDRTree', norm_method="none", pseudo_expr=0,scaling=TRUE)


## black and white heatmap State vs CellType
i_monocle2 <- orderCells(i_monocle2)
st_ce_table = pData(i_monocle2)[, c('State', 'CellType')]
st_ce_table$CellType <- factor(st_ce_table$CellType)
state_cluster_stat = table(st_ce_table)
state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))
state_cluster_stat_ordered <- t(state_cluster_stat)


PH <- pheatmap::pheatmap(state_cluster_stat_ordered, 
                         cluster_cols = F, 
                         cluster_rows = F, 
                         color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(100),
                         fontsize = 30, fontsize_row = 30, fontsize_col = 30,cellwidth=30,cellheight = 30 )


pdf("Bonanomi_trajectory_HeatMap_BW.pdf", 20,10)
PH
dev.off()


#saveRDS(i_monocle2, "Bonanomi_trajectory.rds")

#i_monocle2 = readRDS("Bonanomi_trajectory.rds")

#color by celltype

P1 = plot_cell_trajectory(i_monocle2, color_by = "CellType", cell_name_size = 8) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  theme(legend.text = element_text(size=10)) +
  theme(legend.title = element_text(size=10, face="bold")) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                                         'BARR_END_CAP' = '#336666', 
                                                                                         'CAPILLARY_PLVAP-' = '#399933',
                                                                                         'CAPILLARY_PLVAP+' = '#99CC33',
                                                                                         'TIP_1' = '#6600CC',
                                                                                         'TIP_2' = '#FF99CC',
                                                                                         'TIP_3' = '#FF00FF',
                                                                                         'eight' = 'grey',
                                                                                         'VENOUS_PLVAP-' = '#990000',
                                                                                         'VENOUS_PLVAP+' = '#FF6666'))


pdf("Bonanomi_trajectory_CellType.pdf")
P1
dev.off()

#color by State
    
P3 = plot_cell_trajectory(i_monocle2, color_by = "State")+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold")) + scale_color_manual(values=c('1' = '#535D87',
                                                                                        '2' = '#7C5D67', 
                                                                                        '3' = '#808B96'))

pdf("Bonanomi_trajectory_state_grey_new.pdf")
P3
dev.off()

#color and wrap by CellType     
    
P4 = plot_cell_trajectory(i_monocle2, color_by = "CellType") + facet_wrap(~CellType, nrow = 1)+
      guides(color = guide_legend(override.aes = list(size = 6))) +
      theme(legend.text = element_text(size=20)) +
      theme(text = element_text(size=20)) +
      theme(legend.title = element_text(size=20,face="bold")) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                                            'BARR_END_CAP' = '#336666', 
                                                                                            'CAPILLARY_PLVAP-' = '#399933',
                                                                                            'CAPILLARY_PLVAP+' = '#99CC33',
                                                                                            'TIP_1' = '#6600CC',
                                                                                            'TIP_2' = '#FF99CC',
                                                                                            'TIP_3' = '#FF00FF',
                                                                                            'eight' = 'grey',
                                                                                            'VENOUS_PLVAP-' = '#990000',
                                                                                            'VENOUS_PLVAP+' = '#FF6666'))
    
    
    
    
pdf("Bonanomi_trajectory_CellType_wrap.pdf", 30, 10) 
P4
dev.off()
    
    

# impose trajectory start in proliferating cells for pseudotime   
    
    GM_state <- function(cds){
      if (length(unique(cds@phenoData@data$State)) > 1){
        T0_counts <- table(cds@phenoData@data$State, cds@phenoData@data$CellType)[,"TIP_2"]
        return(as.numeric(names(T0_counts)[which
                                           (T0_counts == max(T0_counts))]))
      } else {
        return (1)
      }
    }
    
    
    i_monocle2 <- orderCells(i_monocle2, root_state = GM_state(i_monocle2))
    
# color by pseudotime
    
P6 = plot_cell_trajectory(i_monocle2, color_by = "Pseudotime") + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))


pdf("Bonanomi_trajectory_Pseudotime.pdf")
P6
dev.off()



my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Car4", "Mfsd2a", "Mki67", "Apln", "Itgb4", "Ackr1")))

cds_subset <- i_monocle2[my_genes,]



pdf("plot_genes_branched_pseudotime_Clust6.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               branch_labels = c("Arterial-Barrier","Venous"),
                               ncol = 1) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                       'BARR_END_CAP' = '#336666', 
                                                                       'CAPILLARY_PLVAP-' = '#399933',
                                                                       'CAPILLARY_PLVAP+' = '#99CC33',
                                                                       'TIP_1' = '#6600CC',
                                                                       'TIP_2' = '#FF99CC',
                                                                       'TIP_3' = '#FF00FF',
                                                                       'eight' = 'grey',
                                                                       'VENOUS_PLVAP-' = '#990000',
                                                                       'VENOUS_PLVAP+' = '#FF6666'))+ 
  
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))
dev.off() 


# use the beam function to generate a table of significance scores for each gene. 
# Genes that score significant are said to be branch-dependent in their expression.
    

i_monocle2_res <- BEAM(i_monocle2, branch_point = 1, cores = 32)
i_monocle2_res <- i_monocle2_res[order(i_monocle2_res$qval),]
i_monocle2_res <- i_monocle2_res[,c("gene_short_name", "pval", "qval")]

write.table(i_monocle2_res, "Bonanomi_BEAM.txt")

#save.image(file = "my_work_space_bonanomi.RData")
#load("my_work_space_bonanomi.RData")

# write a table with gene names reported in heatmap rows

plot_table <- plot_genes_branched_heatmap(i_monocle2[var_genes,], 
                                          branch_point = 1,
                                          cores = 32,
                                          num_clusters = 6,
                                          use_gene_short_name = T,
                                          return_heatmap = TRUE,
                                          show_rownames = T)
write.table(plot_table$annotation_row, 'Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_2000variablegenes_6clust.txt')


# generate heatmap

# small without gene names

pdf("Bonanomi_trajectory_Heatmap_branch_point_2pc_bp1_2000variablegenes_6clust.pdf", 20, 220)
plot_genes_branched_heatmap(i_monocle2[var_genes,],
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 5,
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()   

# big with gene names

pdf("Bonanomi_trajectory_Heatmap_branch_point_2pc_bp1_2000variablegenes_6clust.pdf", 20, 220)
plot_genes_branched_heatmap(i_monocle2[var_genes,],
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 3,
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()   



#saveRDS(br_monocle2, "br_monocle2.RDS")

#write.table(plot_table$annotation_row, 'plot_newheatmap_br_subset_2pc_bp2_Cd34.txt')





