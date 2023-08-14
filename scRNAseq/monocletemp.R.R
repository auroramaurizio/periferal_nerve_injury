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
suppressPackageStartupMessages(library(monocle)) 

integrated = readRDS(file ="/beegfs/scratch/ric.cosr/ric.brendolan/BrendolanA_1280_scRNAseq/dataset/090521/integrated_EC_3TIP.Rds")

#DimPlot(integrated, split.by = "stim")

#DefaultAssay(integrated) = "integrated"

#integrated$CellType = Idents(integrated)

#########################################
#########################################
#########################################


var_genes <- integrated[["integrated"]]@var.features


#Extract data, phenotype data, and feature data from the SeuratObject
#data <- as(as.matrix(integrated@assays$RNA@data), 'sparseMatrix')


#pd <- new('AnnotatedDataFrame', data = integrated@meta.data)


#fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
#fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
#i_monocle2 <- newCellDataSet(data,
#                             phenoData = pd,
#                             featureData = fd,
#                             #lowerDetectionLimit = 0.5,
#                             expressionFamily = uninormal())# since I have already normalized, thresholded and scaled data

#var_genes <- integrated[["integrated"]]@var.features
#ordering_genes <- var_genes

#i_monocle2 <- setOrderingFilter(i_monocle2, ordering_genes)

#i_monocle2 <- reduceDimension(i_monocle2, max_components = 2,
#                              method = 'DDRTree', norm_method="none", pseudo_expr=0,scaling=TRUE)

## order cells change colors and theta to match your plot
#i_monocle2 <- orderCells(i_monocle2)



#prova = pData(i_monocle2)[, c('State', 'CellType')]

#prova$CellType <- factor(prova$CellType)

#state_cluster_stat = table(prova)

#state_cluster_stat <- apply(state_cluster_stat, 2, function(x) x / sum(x))

#state_cluster_stat_ordered <- t(state_cluster_stat)


#detailed_cell_type_color <- c("Bas" = "#E088B8", "DC" = "#46C7EF", "Eos" = "#EFAD1E", "Ery" = "#8CB3DF", "M" = "#53C0AD", "MP/EP" = "#4EB859", "GMP" = "#D097C4", "MK" = "#ACC436", "Neu" = "#F5918A")
#annotation_colors = list(`Lineage` = detailed_cell_type_color)
#options(repr.plot.width=4, repr.plot.height=4)


#PH <- pheatmap::pheatmap(state_cluster_stat_ordered, 
#                         cluster_cols = F, 
#                         cluster_rows = F, 
#                         color = colorRampPalette(RColorBrewer::brewer.pal(n=9, name='Greys'))(100),
#                         fontsize = 30, fontsize_row = 30, fontsize_col = 30,cellwidth=30,cellheight = 30 )


#pdf("Bonanomi_trajectory_HeatMap_BW.pdf", 20,10)
#PH
#dev.off()


#saveRDS(i_monocle2, "Bonanomi_trajectory_assayRNA.rds")

i_monocle2 = readRDS("Bonanomi_trajectory_assayRNA.rds")

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









#pdf("Bonanomi_trajectory_CellType.pdf")
#P1
#dev.off()




#pdf("Bonanomi_trajectory_CellType.pdf")


P3 = plot_cell_trajectory(i_monocle2, color_by = "State")+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold")) + scale_color_manual(values=c('1' = '#828282',
                                                                                        '2' = '#B9BBB6', 
                                                                                        '3' = '#222021'))



P3 = plot_cell_trajectory(i_monocle2, color_by = "State")+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold")) + scale_color_manual(values=c('1' = '#535D87',
                                                                                        '2' = '#7C5D67', 
                                                                                        '3' = '#808B96'))

#questa va bene
pdf("Bonanomi_trajectory_state_grey_new.pdf")
P3
dev.off()

getwd()




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

#pdf("Bonanomi_trajectory_state.pdf")
#P3
#dev.off()





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

P6 = plot_cell_trajectory(i_monocle2, color_by = "Pseudotime") + 
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))


pdf("Bonanomi_trajectory_Pseudotime_.pdf")
P6
dev.off()




i_monocle2_res <- BEAM(i_monocle2, branch_point = 1, cores = 32)
i_monocle2_res <- i_monocle2_res[order(i_monocle2_res$qval),]
i_monocle2_res <- i_monocle2_res[,c("gene_short_name", "pval", "qval")]

write.table(i_monocle2_res, "Bonanomi_BEAM.txt")

head(i_monocle2_res)

sub = i_monocle2_res[i_monocle2_res$qval < 1e-4,]

Arterial
sub[grepl("Gja4", sub$gene_short_name),]
sub[grepl("Sema3g", sub$gene_short_name),]
sub[grepl("Sema3g", sub$gene_short_name),]

Venous
sub[grepl("Plvap", sub$gene_short_name),]
sub[grepl("Icam1", sub$gene_short_name),]
sub[grepl("Vwf", sub$gene_short_name),]   

Barrier
sub[grepl("Slc2a1", sub$gene_short_name),]
sub[grepl("Cxcl12", sub$gene_short_name),]
sub[grepl("Mal", sub$gene_short_name),]

my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Car4", "Mfsd2a", "Mki67","Apln","Itgb4","Ackr1")))

cds_subset <- i_monocle2[my_genes,]

pdf("genes_colorgradient_Barrier.pdf")
plot_cell_trajectory(i_monocle2, markers = c("Egr1", "Ret", "Ackr3"), 
                     use_color_gradient = TRUE)
dev.off()

pdf("genes_pseudotime_barrier.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)
dev.off()


pdf("new_genes_fig1.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               ncol = 1) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                       'BARR_END_CAP' = '#336666', 
                                                                       'CAPILLARY_PLVAP-' = '#399933',
                                                                       'CAPILLARY_PLVAP+' = '#99CC33',
                                                                       'TIP_1' = '#6600CC',
                                                                       'TIP_2' = '#FF99CC',
                                                                       'TIP_3' = '#FF00FF',
                                                                       'eight' = 'grey',
                                                                       'VENOUS_PLVAP-' = '#990000',
                                                                       'VENOUS_PLVAP+' = '#FF6666'))
dev.off() 

pdf("genes_state_barrier.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1) 
dev.off() 

##########################

my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Apln", "Cd276")))

cds_subset <- i_monocle2[my_genes,]

pdf("genes_colorgradient_TIP.pdf")
plot_cell_trajectory(i_monocle2, markers = my_genes, 
                     use_color_gradient = TRUE)
dev.off()

pdf("genes_pseudotime_TIP.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)
dev.off()


pdf("genes_celltype_TIP_APLN_Cd276.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               ncol = 1) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                       'BARR_END_CAP' = '#336666', 
                                                                       'CAPILLARY_PLVAP-' = '#399933',
                                                                       'CAPILLARY_PLVAP+' = '#99CC33',
                                                                       'TIP_1' = '#6600CC',
                                                                       'TIP_2' = '#FF99CC',
                                                                       'TIP_3' = '#FF00FF',
                                                                       'eight' = 'grey',
                                                                       'VENOUS_PLVAP-' = '#990000',
                                                                       'VENOUS_PLVAP+' = '#FF6666'))
dev.off() 

pdf("genes_state_TIP.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1) 
dev.off() 


##########################

save.image(file = "my_work_space_bonanomi.RData")
load("my_work_space_bonanomi.RData")


#Venous

my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Aurkb", "Mki67")))

cds_subset <- i_monocle2[my_genes,]



pdf("genes_pseudotime_Venous.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)
dev.off()

pdf("genes_colorgradient_Venous.pdf")
plot_cell_trajectory(i_monocle2, markers = c("Plvap", "Icam1", "Vwf"), 
                     use_color_gradient = TRUE)
dev.off()


pdf("genes_celltype_Venous.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               ncol = 1) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                       'BARR_END_CAP' = '#336666', 
                                                                       'CAPILLARY_PLVAP-' = '#399933',
                                                                       'CAPILLARY_PLVAP+' = '#99CC33',
                                                                       'TIP_1' = '#6600CC',
                                                                       'TIP_2' = '#FF99CC',
                                                                       'TIP_3' = '#FF00FF',
                                                                       'eight' = 'grey',
                                                                       'VENOUS_PLVAP-' = '#990000',
                                                                       'VENOUS_PLVAP+' = '#FF6666'))
dev.off() 



pdf("genes_celltype_Venous.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               ncol = 1) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                       'BARR_END_CAP' = '#336666', 
                                                                       'CAPILLARY_PLVAP-' = '#399933',
                                                                       'CAPILLARY_PLVAP+' = '#99CC33',
                                                                       'TIP_1' = '#6600CC',
                                                                       'TIP_2' = '#FF99CC',
                                                                       'TIP_3' = '#FF00FF',
                                                                       'eight' = 'grey',
                                                                       'VENOUS_PLVAP-' = '#990000',
                                                                       'VENOUS_PLVAP+' = '#FF6666'))
dev.off() 



pdf("genes_state_Venous.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1) 
dev.off() 



##########################

#Arterial


my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Gja4", "Sema3g", "Stmn2")))

cds_subset <- i_monocle2[my_genes,]



pdf("genes_pseudotime_Arterial.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)
dev.off()





pdf("genes_celltype_Arterial.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               ncol = 1) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
                                                                       'BARR_END_CAP' = '#336666', 
                                                                       'CAPILLARY_PLVAP-' = '#399933',
                                                                       'CAPILLARY_PLVAP+' = '#99CC33',
                                                                       'TIP_1' = '#6600CC',
                                                                       'TIP_2' = '#FF99CC',
                                                                       'TIP_3' = '#FF00FF',
                                                                       'eight' = 'grey',
                                                                       'VENOUS_PLVAP-' = '#990000',
                                                                       'VENOUS_PLVAP+' = '#FF6666'))
dev.off() 


#questi OK

my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Stat1", "Slc2a1", "Mfsd2a")))

cds_subset <- i_monocle2[my_genes,]



pdf("plot_genes_branched_pseudotime_Clust1.pdf")
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


########for paper#####


?plot_genes_branched_pseudotime

pdf("plot_genes_branched_pseudotime_in_column_ordered_fig1.pdf",6,8)
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               panel_order = top-to-bottom,
                               branch_labels = c("1","2"),
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
  
  theme(legend.text = element_text(size=15)) +
  theme(text = element_text(size=15)) +
  theme(legend.title = element_text(size=15,face="bold"))
dev.off() 



my_genes <- row.names(subset(fData(i_monocle2),
                             #gene_short_name %in% c("Sox17", "Mfsd2a", "Mki67", "Apln", "Itgb4", "Vcam1")))
                             gene_short_name %in% c("Car4", "Mfsd2a", "Mki67","Apln","Itgb4","Ackr1")))




cds_subset <- i_monocle2[my_genes,]



cell_size

pdf("plot_genes_branched_pseudotime_in_column_new.pdf", 6, 8)
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               branch_labels = c("1","2"),
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
  
  theme(legend.text = element_text(size=15)) +
  theme(text = element_text(size=15)) +
  theme(legend.title = element_text(size=15,face="bold"))
dev.off() 

#1 arterial
#2 venous

pdf("plot_genes_branched_pseudotime_in_2columns.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "CellType",
                               branch_labels = c("1","2"),
                               ncol = 2) + scale_color_manual(values=c('ARTERIAL' = '#0066FF',
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

pdf("plot_genes_branched_pseudotime_in_column_bigger_cell.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               cell_size = 5,
                               color_by = "CellType",
                               branch_labels = c("1","2"),
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




P3 = plot_cell_trajectory(i_monocle2, color_by = "State")+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold"))

####################

#questi OK

my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Vcam1", "Icam1", "Ackr1", "Ccl11" )))

cds_subset <- i_monocle2[my_genes,]


pdf("plot_genes_branched_pseudotime_Clust2.pdf")
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



Sox17
Car4


my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Sox17", "Car4")))

cds_subset <- i_monocle2[my_genes,]

#questi OK

pdf("plot_genes_branched_pseudotime_Clust3.pdf")
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

###########

Mki67
Top2a
Birc5
Plk1
Aurkb


my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Mki67", "Top2a", "Birc5", "Plk1", "Aurkb")))

cds_subset <- i_monocle2[my_genes,]


pdf("plot_genes_branched_pseudotime_Clust4.pdf")
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


##############

Mest 
Chga
Chst1
Apln
Cxcr4

my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Mest", "Chga", "Chst1", "Apln", "Cxcr4")))

cds_subset <- i_monocle2[my_genes,]


pdf("plot_genes_branched_pseudotime_Clust5.pdf")
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




my_genes <- row.names(subset(fData(i_monocle2),
                             gene_short_name %in% c("Vwf", "Itgb4", "Ackr3", "Lgals3")))

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


getwd()
















pdf("genes_state_Arterial.pdf")
plot_genes_branched_pseudotime(cds_subset,
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1) 
dev.off() 

pdf("genes_colorgradient_Arterial.pdf")
plot_cell_trajectory(i_monocle2, markers = c("Gja4", "Sema3g", "Stmn2"), 
                     use_color_gradient = TRUE)
dev.off()


pdf("genes_colorgradient.pdf", 5,5 )
plot_cell_trajectory(i_monocle2, markers = c("Cd276", "Mal", "Gja4, "), 
                     use_color_gradient = TRUE)
dev.off()









pdf("genes_colorgradient_TIP_Apln.pdf", 5,5 )
plot_cell_trajectory(i_monocle2, markers = c("Apln"), 
                     use_color_gradient = TRUE)
dev.off()


pdf("genes_colorgradient_Arterial_Gja4.pdf", 5,5)
plot_cell_trajectory(i_monocle2, markers = c("Gja4"), 
                     use_color_gradient = TRUE)
dev.off()


pdf("genes_colorgradient_TIP_Cd276.pdf", 5,5)
plot_cell_trajectory(i_monocle2, markers = c("Cd276"), 
                     use_color_gradient = TRUE)
dev.off()


pdf("genes_colorgradient_TIP_Apln.pdf", 5,5)
plot_cell_trajectory(i_monocle2, markers = c("Apln"), 
                     use_color_gradient = TRUE)
dev.off()  


plot_cell_trajectory(i_monocle2, markers = c("Stat1"), 
                     use_color_gradient = TRUE)


cl_1	cl_2	cl_3	cl_4	cl_5	cl_6
Stat1	VCAM1	Sox17	Mki67	Mest 	Vwf
Slc2a1	ICAM1	Car4	Top2a	Chga	Itgb4
Mfsd2a	Ackr1		Birc5	Chst1	ACKR3
CCL11		Plk1	Apln	Lgals3
Aurkb	Cxcr4	








plot_genes_branched_pseudotime(i_monocle2["Slc2a1",],
                               branch_point = 1,
                               color_by = "Time",
                               ncol = 1)
1e-4
10000
0.0001

plot_cell_trajectory_plot_heatmap_i_subset_2pc_bp1_Cd34 <- plot_genes_branched_heatmap(i_monocle2[row.names(subset(i_monocle2_res,
                                                                                                                   qval < 1e-4)),],
                                                                                       branch_point = 1,
                                                                                       cores = 32,
                                                                                       num_clusters = 6,
                                                                                       use_gene_short_name = T,
                                                                                       return_heatmap = TRUE,
                                                                                       show_rownames = T)


write.table(plot_cell_trajectory_plot_heatmap_i_subset_2pc_bp1_Cd34$annotation_row, 'Bonanomi_trajectory_Heatmap_2pc_bp1.txt')



pdf("Bonanomi_trajectory_Heatmap_branch_point_2pc_bp1_2000variablegenes_5clust.pdf", 20, 220)
plot_genes_branched_heatmap(i_monocle2[var_genes,],
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 5,
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()   


pdf("Bonanomi_trajectory_Heatmap_branch_point_2pc_bp1_2000variablegenes_4clust.pdf", 20, 220)
plot_genes_branched_heatmap(i_monocle2[var_genes,],
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 4,
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()   


pdf("Bonanomi_trajectory_Heatmap_branch_point_2pc_bp1_2000variablegenes_3clust.pdf", 20, 220)
plot_genes_branched_heatmap(i_monocle2[var_genes,],
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 3,
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()   

#questo qui

pdf("Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_2000variablegenes_6clust.pdf",)
plot_genes_branched_heatmap(i_monocle2[var_genes,],
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 6,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()   



pdf("Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_2000variablegenes_5clust.pdf",)
plot_genes_branched_heatmap(i_monocle2[var_genes,],
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 5,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()   


plot_table <- plot_genes_branched_heatmap(i_monocle2[var_genes,], 
                                          branch_point = 1,
                                          cores = 32,
                                          num_clusters = 6,
                                          use_gene_short_name = T,
                                          return_heatmap = TRUE,
                                          show_rownames = T)
write.table(plot_table$annotation_row, 'Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_2000variablegenes_6clust.txt')



plot_table <- plot_genes_branched_heatmap(i_monocle2[var_genes,], 
                                          branch_point = 1,
                                          cores = 32,
                                          num_clusters = 5,
                                          use_gene_short_name = T,
                                          return_heatmap = TRUE,
                                          show_rownames = T)
write.table(plot_table$annotation_row, 'Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_2000variablegenes_5clust.txt')


plot_table <- plot_genes_branched_heatmap(i_monocle2[var_genes,], 
                                          branch_point = 1,
                                          cores = 32,
                                          num_clusters = 4,
                                          use_gene_short_name = T,
                                          return_heatmap = TRUE,
                                          show_rownames = T)
write.table(plot_table$annotation_row, 'Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_2000variablegenes_4clust.txt')


plot_table <- plot_genes_branched_heatmap(i_monocle2[var_genes,], 
                                          branch_point = 1,
                                          cores = 32,
                                          num_clusters = 3,
                                          use_gene_short_name = T,
                                          return_heatmap = TRUE,
                                          show_rownames = T)
write.table(plot_table$annotation_row, 'Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_2000variablegenes_3clust.txt')





pdf("Bonanomi_trajectory_Heatmap_branch_point_2pc_bp1_gn_false.pdf", 20, 150)
plot_genes_branched_heatmap(i_monocle2[row.names(subset(i_monocle2_res,
                                                        qval < 1e-4)),],
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 6,
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()



pdf("Bonanomi_trajectory_Heatmap_branch_point_2pc_bp1_gn.pdf", 20, 150)
plot_genes_branched_heatmap(i_monocle2[row.names(subset(i_monocle2_res,
                                                        qval < 1e-4)),],
                            branch_point = 1,
                            cores = 32,
                            num_clusters = 6,
                            branch_colors = c("#619CFF","#F8766D","#00BA38"),
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()



P3 = plot_cell_trajectory(i_monocle2, color_by = "State")+
  guides(color = guide_legend(override.aes = list(size = 6))) +
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20,face="bold")) + scale_color_manual(values=c('1' = '#535D87',
                                                                                        '2' = '#7C5D67', 
                                                                                        '3' = '#808B96'))

pdf("Bonanomi_trajectory_Heatmap_branch_point_small_2pc_bp1_newcols_4_5.pdf", 4, 5)
plot_genes_branched_heatmap(i_monocle2[row.names(subset(i_monocle2_res,
                                                        qval < 1e-4)),],
                            branch_colors = c("#808B96","#535D87","#7C5D67"), #blue #orange #green ## # # # 
                            branch_point = 1,
                            cores = 16,
                            num_clusters = 6,
                            use_gene_short_name = T,
                            show_rownames = F)
dev.off()   




plot_table <- plot_multiple_branches_heatmap(i_monocle2,
                                             branches = c(1,2,3),
                                             branches_name = c("1","2","3"),
                                             cores = 16,
                                             return_heatmap = TRUE,
                                             use_gene_short_name = T,
                                             show_rownames = T)







