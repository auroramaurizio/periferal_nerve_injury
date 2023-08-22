library(CellChat)
library(patchwork)
library(Seurat)
library(ComplexHeatmap)
library(R.utils)
options(stringsAsFactors = FALSE)


# The signaling crosstalk between ECs and other cellTypes involved in nerve injury was investigated by merging 
# ECs from our intact/7 dpi scRNA-seq datasets, together with Schwann cells, macrophages, and Mesenchymal cells from other studies 
# {Toma, 2020 ;Ydens, 2020; Carr, 2019  }. 
# In parallel the same analysis was performed with MES cells from our study instead of Carr's. 
# In this script we use Carr's MES.

# create cellchat intact and injury objects as described in https://github.com/sqjin/CellChat

#load injury R object
injury = readRDS("Injury_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")

DefaultAssay(injury) = "RNA"

count_norm_injury = injury@assays$RNA@data # normalized data matrix

meta_data_injury <- cbind(rownames(injury@meta.data),as.data.frame(Idents(injury)))
colnames(meta_data_injury)<- c("Cell","labels")
cellchat_injury <- createCellChat(object = count_norm_injury, meta = meta_data_injury, group.by = "labels")
cellchat_injury <- addMeta(cellchat_injury, meta = meta_data_injury)
cellchat_injury <- setIdent(cellchat_injury, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_injury@idents) # show factor levels
groupSize <- as.numeric(table(cellchat_injury@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB 

# set the used database in the object
cellchat_injury@DB <- CellChatDB.use
cellchat_injury <- subsetData(cellchat_injury) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) 

cellchat_injury <- identifyOverExpressedGenes(cellchat_injury)
cellchat_injury <- identifyOverExpressedInteractions(cellchat_injury)
cellchat_injury <- projectData(cellchat_injury, PPI.mouse)
cellchat_injury <- computeCommunProb(cellchat_injury)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_injury <- filterCommunication(cellchat_injury, min.cells = 10)
cellchat_injury <- computeCommunProbPathway(cellchat_injury)
cellchat_injury <- aggregateNet(cellchat_injury)
groupSize <- as.numeric(table(cellchat_injury@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("injury_only_circleNofI.pdf")
netVisual_circle(cellchat_injury@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("injury_only_circleSofI.pdf")
netVisual_circle(cellchat_injury@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()


mat_injury <- cellchat_injury@net$weight
par(mfrow = c(3,4), xpd=TRUE)

pdf("injury_only.pdf")
for (i in 1:nrow(mat_injury)) {
  mat2 <- matrix(0, nrow = nrow(mat_injury), ncol = ncol(mat_injury), dimnames = dimnames(mat_injury))
  mat2[i, ] <- mat_injury[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_injury), title.name = rownames(mat_injury)[i])
}
dev.off()


#load intact R object

intact = readRDS("intact_EC_MES_MACS_SCHWANN_forCDB_090222_merged_seurat.Rds")
DimPlot(intact, split.by = "orig.ident", ncol = 2)
DefaultAssay(intact) = "RNA"
count_norm_intact = intact@assays$RNA@data # normalized data matrix
meta_data_intact <- cbind(rownames(intact@meta.data),as.data.frame(Idents(intact)))
colnames(meta_data_intact)<- c("Cell","labels")
cellchat_intact <- createCellChat(object = count_norm_intact, meta = meta_data_intact, group.by = "labels")
cellchat_intact <- addMeta(cellchat_intact, meta = meta_data_intact)
cellchat_intact <- setIdent(cellchat_intact, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat_intact@idents) # show factor levels
groupSize <- as.numeric(table(cellchat_intact@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse 
showDatabaseCategory(CellChatDB)

filename_xls <- 'CellChatDB_interactionDB.xlsx'
write.xlsx(CellChatDB$interaction,
           file= filename_xls, 
           row.names = T,
           asTable = T)

(CellChatDB$interaction)
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat_intact@DB <- CellChatDB.use
cellchat_intact <- subsetData(cellchat_intact) # This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) #

cellchat_intact <- identifyOverExpressedGenes(cellchat_intact)
cellchat_intact <- identifyOverExpressedInteractions(cellchat_intact)

cellchat_intact <- projectData(cellchat_intact, PPI.mouse)
cellchat_intact <- computeCommunProb(cellchat_intact)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat_intact <- filterCommunication(cellchat_intact, min.cells = 10)
cellchat_intact <- computeCommunProbPathway(cellchat_intact)
cellchat_intact <- aggregateNet(cellchat_intact)
groupSize <- as.numeric(table(cellchat_intact@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf("intact_only_circleNofI.pdf")
netVisual_circle(cellchat_intact@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

pdf("intact_only_circleSofI.pdf")
netVisual_circle(cellchat_intact@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
mat_intact <- cellchat_intact@net$weight
par(mfrow = c(3,4), xpd=TRUE)

pdf("intact_only.pdf")
for (i in 1:nrow(mat_intact)) {
  mat2 <- matrix(0, nrow = nrow(mat_intact), ncol = ncol(mat_intact), dimnames = dimnames(mat_intact))
  mat2[i, ] <- mat_intact[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat_intact), title.name = rownames(mat_intact)[i])
}
dev.off()

# save CellChat objects

saveRDS(cellchat_intact, "cellchat_intact_bis.Rds")
saveRDS(cellchat_injury, "cellchat_injury_bis.Rds")

# load cellchat intact and injury objects 

cellchat_intact = readRDS("/Users/maurizio.aurora/cellchat_intact_bis.Rds")
cellchat_injury = readRDS("/Users/maurizio.aurora/cellchat_injury_bis.Rds")

# merge object list

object.list <- list(intact = cellchat_intact, injury = cellchat_injury)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")


# plot heatmap on the merged object

gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg1 + gg2


# customize the above heatmap: adjusting the scale of the bars to be comparable among the 2 plots and the colors

#get_color function from https://rdrr.io/github/sqjin/CellChat/api/

# x: vector
# break1 single value
# break2 single value
# rgb1 vector with 3 elements
# rgb2 vector with 3 elements
.get_color = function(x, break1, break2, col1, col2, space) {
  
  col1 = colorspace::coords(as(colorspace::sRGB(col1[1], col1[2], col1[3]), space))
  col2 = colorspace::coords(as(colorspace::sRGB(col2[1], col2[2], col2[3]), space))
  
  res_col = matrix(ncol = 3, nrow = length(x))
  for(j in 1:3) {
    xx = (x - break2)*(col2[j] - col1[j]) / (break2 - break1) + col2[j]
    res_col[, j] = xx
  }
  
  res_col = get(space)(res_col)
  res_col = colorspace::coords(as(res_col, "sRGB"))
  res_col[, 1] = .restrict_in(res_col[,1], 0, 1)
  res_col[, 2] = .restrict_in(res_col[,2], 0, 1)
  res_col[, 3] = .restrict_in(res_col[,3], 0, 1)
  colorspace::hex(colorspace::sRGB(res_col))
}



gg1 <- netVisual_heatmap(cellchat_intact)
gg2 <- netVisual_heatmap(cellchat_injury)

gg2_orig = gg2
gg1_orig = gg1

#----------------------------------------
gg1@matrix[is.na(gg1@matrix)] <- 0
gg2@matrix[is.na(gg2@matrix)] <- 0
#---------------------------------------
matrix_color_mapping <- gg2@matrix_color_mapping
matrix_color_mapping@colors <- colorRampPalette(c('dodgerblue4','white','darkred'))(5)
gg1@matrix_color_mapping <- gg2@matrix_color_mapping <- matrix_color_mapping
#---------------------------------------
# gg1@right_annotation@anno_list$Strength@fun@data_scale <- gg2@right_annotation@anno_list$Strength@fun@data_scale
# gg1@top_annotation@anno_list$Strength@fun@data_scale <- gg2@top_annotation@anno_list$Strength@fun@data_scale
#---------------------------------------
gg1 + gg2


mat <- gg2_orig@matrix
mat[is.na(mat)] <- 0
legend.name <- 'Number of interactions'
color.heatmap.use <- colorRampPalette(c('dodgerblue4','snow','darkred'))(5)
font.size = 8
font.size.title = 10
title.name = NULL
cluster.rows = TRUE; cluster.cols = TRUE
color.use <- c(
  'VENOUS_PLVAP+' = '#F8766D', 
  'VENOUS_PLVAP-' = 'brown',
  'TIP_1' = '#CD9600',
  'TIP_2' = 'pink',
  'BARR_END_CAP' = '#7CAE00',
  'CAPILLARY_PLVAP-' = '#00BE67', 
  'ARTERIAL' = '#00BFC4',
  'MACROPHAGES' = '#00A9FF', 
  'MES_DIFF'='#C97398',
  'MES_ENDONEUR'='#7DF3FF',
  'MES_EPINEUR'='#3993D0',
  'MES_PERINEUR'='#4339D0',
  'LEC'='#FF61CC', 
  'CAPILLARY_PLVAP+' = 'grey', 
  'SCHWANN' = '#DCE961',
  'TIP_3' = "purple"
)[rownames(mat)]


df<- data.frame(group = colnames(mat)); rownames(df) <- colnames(mat)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))
row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))

ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use, col=color.use)), show_annotation_name = FALSE)



ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
              show_column_dend = F, show_row_dend = F,
              bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
              cluster_rows = cluster.rows,cluster_columns = cluster.rows,
              row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
              # width = unit(width, "cm"), height = unit(height, "cm"),
              column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
              row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                          border = NA, #at = colorbar.break,
                                          legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
)
(injury_heatmap <- ht1)


mat <- gg1_orig@matrix
mat[is.na(mat)] <- 0
legend.name <- 'Number of interactions'
color.heatmap.use <- colorRampPalette(c('dodgerblue4','white','darkred'))(5)
font.size = 8
font.size.title = 10
title.name = NULL
cluster.rows = FALSE; cluster.cols = FALSE
color.use <- c(
  'VENOUS_PLVAP+' = '#F8766D', 
  'VENOUS_PLVAP-' = 'brown',
  'TIP_1' = '#CD9600',
  'TIP_2' = 'pink',
  'BARR_END_CAP' = '#7CAE00',
  'CAPILLARY_PLVAP-' = '#00BE67', 
  'ARTERIAL' = '#00BFC4',
  'MACROPHAGES' = '#00A9FF', 
  'MES_DIFF'='#C97398',
  'MES_ENDONEUR'='#7DF3FF',
  'MES_EPINEUR'='#3993D0',
  'MES_PERINEUR'='#4339D0',
  'LEC'='#FF61CC', 
  'CAPILLARY_PLVAP+' = 'grey', 
  'SCHWANN' = '#DCE961',
  'TIP_3' = "purple"
)[rownames(mat)]



df<- data.frame(group = colnames(mat)); 
rownames(df) <- colnames(mat)
col_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use),which = "column",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))
row_annotation <- HeatmapAnnotation(df = df, col = list(group = color.use), which = "row",
                                    show_legend = FALSE, show_annotation_name = FALSE,
                                    simple_anno_size = grid::unit(0.2, "cm"))

ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), 
                                            ylim = injury_heatmap@right_annotation@anno_list$Strength@fun@data_scale, 
                                            border = FALSE, gp = gpar(fill = color.use, col=color.use)), 
                    show_annotation_name = FALSE)
ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), 
                                                ylim = injury_heatmap@top_annotation@anno_list$Strength@fun@data_scale, 
                                                border = FALSE,gp = gpar(fill = color.use, col=color.use)), 
                        show_annotation_name = FALSE)


ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
              bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
              cluster_rows = cluster.rows,cluster_columns = cluster.rows,
              row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
              # width = unit(width, "cm"), height = unit(height, "cm"),
              column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
              row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
              heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                          border = NA, #at = colorbar.break,
                                          legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
)

ht1@column_order <- c(12,13,1,3,14,11,15,7,9,4,10,2,8,6,5)
ht1@row_order <- c(7,9,6,12,13,3,1,11,14,15,8,4,2,10,5)
intact_heatmap <- ht1
intact_heatmap@matrix_color_mapping <- injury_heatmap@matrix_color_mapping


pdf('cellChat_interactions_heatmap_clustered.pdf', height = 4)
intact_heatmap + injury_heatmap
dev.off()


# make circos



# make differential

# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "injury"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "injury",ligand.logFC = 0.2, receptor.logFC = 0.2)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "intact",ligand.logFC = -0.2, receptor.logFC = -0.2)



filename_xls <- 'net.down_0.2.xlsx'

write.xlsx(net.down,
           file= filename_xls,
           row.names = F,
           asTable = T)



filename_xls <- 'net.up_0.2.xlsx'

write.xlsx(net.up,
           file= filename_xls,
           row.names = F,
           asTable = T)

#Identify dysfunctional signaling by using differential expression analysis

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)


# subsetCommunication returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors.

df.net.cellchat_intact_LR <- subsetCommunication(cellchat_intact)
df.net.cellchat_injury_LR <- subsetCommunication(cellchat_injury)


filename_xls <- 'df.net.cellchat_intact_LR.xlsx'

write.xlsx(df.net.cellchat_intact_LR,
           file= filename_xls,
           row.names = F,
           asTable = T)



filename_xls <- 'df.net.cellchat_injury_LR.xlsx'

write.xlsx(df.net.cellchat_injury_LR,
           file= filename_xls,
           row.names = F,
           asTable = T)


# set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

df.net.cellchat_intact_sigpat <- subsetCommunication(cellchat_intact, slot.name = "netP" )
df.net.cellchat_injury_sigpat <- subsetCommunication(cellchat_injury, slot.name = "netP" )
df.net.cellchat_sigpat <- subsetCommunication(cellchat, slot.name = "netP" )


filename_xls <- 'df.net.cellchat_intact_sigpat.xlsx'

write.xlsx(df.net.cellchat_intact_sigpat,
           file= filename_xls,
           row.names = F,
           asTable = T)



filename_xls <- 'df.net.cellchat_injury_sigpat.xlsx'

write.xlsx(df.net.cellchat_injury_sigpat,
           file= filename_xls,
           row.names = F,
           asTable = T)


# Identify signaling groups based on their functional similarity

cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")

# Compute and visualize the pathway distance in the learned joint manifold

pdf("rankSimilarity_MES_intact_inj_2.0.pdf")
rankSimilarity(cellchat, type = "functional")
dev.off()

# Compare the overall information flow of each signaling pathway

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
pdf("CompareInformationFlow_MES_intact_inj_2.0.pdf",6, 8)
gg1 + gg2
dev.off()


# Compute and visualize the network centrality scores

for (i in 1:length(object.list)) {
  object.list[[i]] = netAnalysis_computeCentrality(object.list[[i]], slot.name = "netP") 
}

cellchat_intact <- netAnalysis_computeCentrality(cellchat_intact, slot.name = "netP") 
cellchat_injury <- netAnalysis_computeCentrality(cellchat_injury, slot.name = "netP") 


#Visualize the dominant senders (sources) and receivers (targets) in a 2D space


for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]],  title = names(object.list)[i], weight.MinMax = weight.MinMax)
}


gg1 <- netAnalysis_signalingRole_scatter(cellchat_injury, color.use = c('VENOUS_PLVAP+' = '#F8766D',
                                                                        'VENOUS_PLVAP-' = 'brown',
                                                                        'ARTERIAL' = '#00BFC4',
                                                                        'TIP_1' = '#CD9600',
                                                                        'TIP_2' = 'pink',
                                                                        'BARR_END_CAP' = '#7CAE00',
                                                                        'CAPILLARY_PLVAP-' = '#00BE67',
                                                                        'MACROPHAGES' = '#00A9FF',
                                                                        'MES_DIFF'='#C97398',
                                                                        'MES_ENDONEUR'='#7DF3FF',
                                                                        'MES_EPINEUR'='#3993D0',
                                                                        'MES_PERINEUR'='#4339D0',
                                                                        'LEC'='#FF61CC',
                                                                        'CAPILLARY_PLVAP+' = 'grey',
                                                                        'SCHWANN' = '#DCE961',
                                                                        'TIP_3' = 'purple'))  + xlim(0, 30) + ylim(0, 30)



gg2 <- netAnalysis_signalingRole_scatter(cellchat_intact, color.use = c('VENOUS_PLVAP+' = '#F8766D',
                                                                        'VENOUS_PLVAP-' = 'brown',
                                                                        'ARTERIAL' = '#00BFC4',
                                                                        'TIP_1' = '#CD9600',
                                                                        'TIP_2' = 'pink',
                                                                        'BARR_END_CAP' = '#7CAE00',
                                                                        'CAPILLARY_PLVAP-' = '#00BE67',
                                                                        'MACROPHAGES' = '#00A9FF',
                                                                        'MES_DIFF'='#C97398',
                                                                        'MES_ENDONEUR'='#7DF3FF',
                                                                        'MES_EPINEUR'='#3993D0',
                                                                        'MES_PERINEUR'='#4339D0',
                                                                        'LEC'='#FF61CC',
                                                                        'CAPILLARY_PLVAP+' = 'grey',
                                                                        'SCHWANN' = '#DCE961',
                                                                        'TIP_3' = 'purple'))  + xlim(0, 30) + ylim(0, 30)


pdf("diagonal_injury.pdf")
gg1
dev.off()




pdf("diagonal_intact.pdf")
gg2
dev.off()



#Chord diagrams


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')

netVisual_chord_gene(cellchat, sources.use = "ARTERIAL", targets.use = c(1:15), lab.cex = 0.5,legend.pos.y = 30)

colors = c('VENOUS_PLVAP+' = '#F8766D',
           'VENOUS_PLVAP-' = 'brown',
           'ARTERIAL' = '#00BFC4',
           'TIP_1' = '#CD9600',
           'TIP_2' = 'pink',
           'BARR_END_CAP' = '#7CAE00',
           'CAPILLARY_PLVAP-' = '#00BE67',
           'MACROPHAGES' = '#00A9FF',
           'MES_DIFF'='#C97398',
           'MES_ENDONEUR'='#7DF3FF',
           'MES_EPINEUR'='#3993D0',
           'MES_PERINEUR'='#4339D0',
           'LEC'='#FF61CC',
           'CAPILLARY_PLVAP+' = 'grey',
           'SCHWANN' = '#DCE961',
           'TIP_3' = 'purple')


pathways.show <- c("SEMA3")

pdf(file ="SEMA3C_PLXND1_LR_Chorddiagram.pdf", width = 20, height =16)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]],  color.use = colors,
                       pairLR.use = data.frame(interaction_name='SEMA3C_PLXND1'))
}
dev.off()


pdf(file ="SEMA3_Chorddiagram_gene.pdf", width = 8, height =8)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], signaling = pathways.show, color.use = colors, title.name = paste(pathways.show, names(object.list)[i]),show.legend = FALSE)
}
dev.off()


pdf(file ="SEMA3_Chorddiagram_CellType.pdf", width = 8, height =8)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], pt.title	= 25, signaling = pathways.show, layout = "chord", vertex.label.cex = 1, signaling.name = paste(pathways.show, names(object.list)[i]), color.use = colors)
}
dev.off()



