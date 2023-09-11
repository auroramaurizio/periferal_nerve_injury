########## code for the interaction pairs heatmap (Suppl Fig 5 U and V)

library(clusterProfiler)
library(org.Mm.eg.db)
library(reshape2)
library(pheatmap)
library(openxlsx)
library(stringr)
library(Seurat)
library(ggplot2)
library(ggrepel)
library(patchwork)



view_pair <- function(mat_lr, source_cell, target_cell) {
  pair_mat <- acast(subset(mat_lr, source==source_cell & target==target_cell), 
                    ligand ~ receptor, value.var = 'prob')
  pair_mat_fig <- pair_mat
  pair_mat[is.na(pair_mat)] <- 1
  rowOrd <- hclust(dist(pair_mat))
  colOrd <- hclust(dist(t(pair_mat)))
  pheatmap(-log10(pair_mat_fig), cluster_rows = rowOrd, 
           cluster_cols = colOrd, 
           color = colorRampPalette(c('snow', 'pink','red','darkred','black'))(100))
}

view_pair_diff <- function(mat_lr_1, mat_lr_2, source_cell, target_cell) {
  pair_mat_1 <- acast(subset(mat_lr_1, source==source_cell & target==target_cell), 
                      ligand ~ receptor, value.var = 'prob')
  pair_mat_2 <- acast(subset(mat_lr_2, source==source_cell & target==target_cell), 
                      ligand ~ receptor, value.var = 'prob')
  new_row <- unique(c(rownames(pair_mat_1), rownames(pair_mat_2)))
  new_col <- unique(c(colnames(pair_mat_1), colnames(pair_mat_2)))
  diff_mat_1 <- diff_mat_2 <- matrix(1, nrow = length(new_row), ncol = length(new_col), 
                                     dimnames = list(new_row, new_col))
  pair_mat_1[is.na(pair_mat_1)] <- 1
  pair_mat_2[is.na(pair_mat_2)] <- 1
  diff_mat_1[rownames(pair_mat_1), colnames(pair_mat_1)] <- pair_mat_1
  diff_mat_2[rownames(pair_mat_2), colnames(pair_mat_2)] <- pair_mat_2
  values <- log10(diff_mat_2/diff_mat_1)
  breaks <- seq(min(values), max(values), length.out = 101)
  sym_id <- which.min((breaks + breaks[1])^2)
  palette <- c(
    colorRampPalette(c('blue','lightblue','snow', 'pink','red'))(sym_id),
    colorRampPalette(c('red','darkred','black'))(101 - sym_id + 1)
  )
  pheatmap(values, breaks = breaks, col = palette, 
           main = paste(source_cell,'->',target_cell))
}
#



## Load sc datasets


sc_ec_inj <- subset(readRDS('../integrated_EC_3TIP_clstrmvd.Rds'), stim=='2')
sc_mes_inj <- subset(readRDS('../integrated_MES.rds'), orig.ident=='INJURY_D7')
avg_ec_inj <- AverageExpression(sc_ec_inj)$RNA
avg_mes_inj <- AverageExpression(sc_mes_inj)$RNA


## Load interactions


carr_lr <- read.xlsx('../cellChat/raw/df.net.cellchat_injury_LR.xlsx', check.names = FALSE)
pdgfrb_lr <- read.xlsx('raw/df.net.cellchat_injury_LR_MES.xlsx', check.names = FALSE)
pdgfrb_intact_lr <- read.xlsx('raw/df.net.cellchat_intact_LR_MES.xlsx', check.names = FALSE)
head(pdgfrb_lr)


## obtain ligands and receptors from interaction names



interaction_names <- unique(c(carr_lr$interaction_name_2, 
                              pdgfrb_lr$interaction_name_2, 
                              pdgfrb_intact_lr$interaction_name_2))
inferred_genes <- do.call('rbind', 
                          strsplit(interaction_names, '[ ]*- '))
inferred_receptors <- strsplit(inferred_genes[,2], '\\+')
inferred_receptors <- do.call('rbind', 
                              lapply(inferred_receptors, function(x) 
                                if(length(x)==1) c(x, '') else x))
inferred_receptors[,1] <- sub('\\(', '', inferred_receptors[,1])
inferred_receptors[,2] <- sub('\\)', '', inferred_receptors[,2])
gene_ann <- data.frame(
  lig=inferred_genes[,1], 
  r1=inferred_receptors[,1],
  r2=inferred_receptors[,2], 
  row.names = interaction_names)


## remove cell lines not in common (MES_DIVIDING, MES_DIFF\*, UNDET)


setdiff(unique(c(pdgfrb_lr$source, pdgfrb_lr$target)), 
        unique(c(carr_lr$source, carr_lr$target)))



common_cell_types <- intersect(unique(c(pdgfrb_lr$source, pdgfrb_lr$target)), 
                               unique(c(carr_lr$source, carr_lr$target)))
pdgfrb_lr <- pdgfrb_lr[pdgfrb_lr$source %in% common_cell_types & pdgfrb_lr$target %in% common_cell_types,]
carr_lr <- carr_lr[carr_lr$source %in% common_cell_types & carr_lr$target %in% common_cell_types,]
pdgfrb_intact_lr <- pdgfrb_intact_lr[pdgfrb_intact_lr$source %in% common_cell_types & pdgfrb_intact_lr$target %in% common_cell_types,]


## number of interactions


pdf('figures/cellChat_interaction_numbers.pdf', height = 4, width = 4)
par(mar=c(5,4,1,1)+.1)
int_num <- c(intact=nrow(pdgfrb_intact_lr),carr=nrow(carr_lr), pdgfrb=nrow(pdgfrb_lr))
bpOut <- barplot(int_num, ylim = c(0,max(int_num)*1.15), las=2)
text(bpOut[,1], int_num, int_num, pos = 3)
dev.off()


## interactions by pairs


carr_int_mat <- acast(carr_lr, source ~ target)
pdgfrb_int_mat <- acast(pdgfrb_lr, source ~ target)[rownames(carr_int_mat), colnames(carr_int_mat)]

#### Difference 


values <- pdgfrb_int_mat-carr_int_mat
values[!is.finite(values)] <- NA
breaks <- seq(min(values, na.rm = T), max(values, na.rm = T), length.out = 101)
palette <- colorRampPalette(c('snow', 'snow','red','darkred','black'))(100)
pheatmap(values, breaks = breaks, col = palette, treeheight_col = 20, treeheight_row = 20,display_numbers = TRUE,number_format = "%.0f",
         filename = 'figures/interaction_analysis_diff_withVals.pdf', height = 5, width = 5.5)


## PDGFRb vs CARR

my_pairs <- c()
for(snd in c('TIP_1','TIP_2','TIP_3','VENOUS_PLVAP+')) 
  for(rcv in c('MES_DIFF','MES_EPINEUR','MES_PERINEUR','MES_ENDONEUR'))
    my_pairs <- c(my_pairs, paste(snd, '->', rcv))

pdgfrb_my_pairs <- pdgfrb_lr
pdgfrb_my_pairs$pair <- paste(pdgfrb_my_pairs$source, '->', 
                              pdgfrb_my_pairs$target)
pdgfrb_my_pairs <- pdgfrb_my_pairs[pdgfrb_my_pairs$pair %in% my_pairs,]
pdgfrb_my_mat <- acast(pdgfrb_my_pairs, interaction_name_2 ~ pair, value.var = 'prob')
pdgfrb_my_mat[is.na(pdgfrb_my_mat)] <- 1

carr_my_pairs <- carr_lr
carr_my_pairs$pair <- paste(carr_my_pairs$source, '->', 
                            carr_my_pairs$target)
carr_my_pairs <- carr_my_pairs[carr_my_pairs$pair %in% my_pairs,]
carr_my_mat <- acast(carr_my_pairs, interaction_name_2 ~ pair, value.var = 'prob')
carr_my_mat[is.na(carr_my_mat)] <- 1

intact_my_pairs <- pdgfrb_intact_lr
intact_my_pairs$pair <- paste(intact_my_pairs$source, '->', 
                              intact_my_pairs$target)
intact_my_pairs <- intact_my_pairs[intact_my_pairs$pair %in% my_pairs,]
intact_my_mat <- acast(intact_my_pairs, interaction_name_2 ~ pair, value.var = 'prob')
intact_my_mat[is.na(intact_my_mat)] <- 1

all_interactions <- unique(c(rownames(pdgfrb_my_mat),rownames(carr_my_mat), rownames(intact_my_mat)))
intact_all_mat <- carr_all_mat <- pdgfrb_all_mat <- 
  matrix(1, length(all_interactions), ncol(pdgfrb_my_mat),
         dimnames = list(all_interactions, colnames(carr_my_mat)))
intact_all_mat[rownames(intact_my_mat), colnames(intact_my_mat)] <- intact_my_mat
carr_all_mat[rownames(carr_my_mat), colnames(carr_my_mat)] <- carr_my_mat
pdgfrb_all_mat[rownames(pdgfrb_my_mat), colnames(pdgfrb_my_mat)] <- pdgfrb_my_mat

all_pairs <- rbind(carr_my_pairs, intact_my_pairs, pdgfrb_my_pairs)
ann_df <- unique(all_pairs[,c('interaction_name_2','pathway_name')])
rownames(ann_df) <- ann_df$interaction_name_2
ann_df$interaction_name_2 <- NULL

values <- log10(carr_all_mat/pdgfrb_all_mat)

ec_genes <- gene_ann[rownames(values),'lig']
ec_cells <- sapply(strsplit(colnames(values), ' -> '), '[', 1)
mes_genes_r1 <- gene_ann[rownames(values),'r1']
mes_genes_r2 <- gene_ann[rownames(values),'r2']
mes_genes_r2[mes_genes_r2==''] <- 'not_a_gene'
mes_cells <- sapply(strsplit(colnames(values), ' -> '), '[', 2)

if(!'not_a_gene' %in% rownames(avg_mes_inj)) {
  avg_mes_inj <- rbind(avg_mes_inj, 'not_a_gene' = rep(1, ncol(avg_mes_inj)))
}

expr_mask_ec <- avg_ec_inj[ec_genes, ec_cells]
expr_mask_mes_r1 <- avg_mes_inj[mes_genes_r1, mes_cells]
expr_mask_mes_r2 <- avg_mes_inj[mes_genes_r2, mes_cells]

# geni per cellula selezionati a diverse soglie di espressione
# (simile sia per EC che per MES)
###   0.3 ~ 4.5k
###   0.5 ~ 2.5k
###   1.0 ~ 1.5k

e_thr <- .5
e_idx <- apply(expr_mask_ec >= e_thr & 
                 expr_mask_mes_r1 >= e_thr & 
                 expr_mask_mes_r2 >= e_thr, 1, any)

values_filt <- values[e_idx,]

breaks <- seq(min(values_filt), max(values_filt), length.out = 101)
sym_id <- which.min((breaks + breaks[1])^2)
palette <- c(
  colorRampPalette(c('blue','lightblue','snow', 'pink','red'))(sym_id),
  colorRampPalette(c('red','darkred','black'))(101 - sym_id + 1)
)

phOut <- pheatmap(values_filt, breaks = breaks,
                  color = palette, annotation_row = ann_df,
                  filename = 'figures/PDGFRb_vs_CARR_interaction.pdf', 
                  main = 'PDGFRb vs CARR',
                  height = 12, width = 8)



## PDGFRb vs INTACT


my_cols_palette <- c(
  '#F0A3FF', #Amethyst
           '#0075DC', #Blue
           '#993F00', #Caramel
           '#4C005C', #Damson
           '#191919', #Ebony
           '#005C31', #Forest
           '#2BCE48', #Green
           '#FFCC99', #Honeydew
           '#808080', #Iron
           '#94FFB5', #Jade
           '#8F7C00', #Khaki
           '#9DCC00', #Lime
           '#C20088', #Mallow
           '#003380', #Navy
           '#FFA405', #Orpiment
           '#FFA8BB', #Pink
           '#426600', #Quagmire
           '#FF0010', #Red
           '#5EF1F2', #Sky
           '#00998F', #Turquoise
           '#E0FF66', #Uranium
           '#740AFF', #Violet
           '#990000', #Wine
           '#FFFF80', #Xanthin
           '#FFFF00', #yellow
           '#FF5005' #Zinnia
)


my_pairs <- c()
for(snd in c('TIP_1','TIP_2','TIP_3')) 
  for(rcv in c('MES_DIFF','MES_EPINEUR','MES_PERINEUR','MES_ENDONEUR'))
    my_pairs <- c(my_pairs, paste(snd, '->', rcv))

pdgfrb_my_pairs <- pdgfrb_lr
pdgfrb_my_pairs$pair <- paste(pdgfrb_my_pairs$source, '->', 
                              pdgfrb_my_pairs$target)
pdgfrb_my_pairs <- pdgfrb_my_pairs[pdgfrb_my_pairs$pair %in% my_pairs,]
pdgfrb_my_mat <- acast(pdgfrb_my_pairs, interaction_name_2 ~ pair, value.var = 'prob')
pdgfrb_my_mat[is.na(pdgfrb_my_mat)] <- 1

intact_my_pairs <- pdgfrb_intact_lr
intact_my_pairs$pair <- paste(intact_my_pairs$source, '->', 
                              intact_my_pairs$target)
intact_my_pairs <- intact_my_pairs[intact_my_pairs$pair %in% my_pairs,]
intact_my_mat <- acast(intact_my_pairs, interaction_name_2 ~ pair, value.var = 'prob')
intact_my_mat[is.na(intact_my_mat)] <- 1

all_interactions <- unique(c(rownames(pdgfrb_my_mat), rownames(intact_my_mat)))
intact_all_mat <- pdgfrb_all_mat <- 
  matrix(1, length(all_interactions), ncol(pdgfrb_my_mat),
         dimnames = list(all_interactions, colnames(pdgfrb_my_mat)))
intact_all_mat[rownames(intact_my_mat), colnames(intact_my_mat)] <- intact_my_mat
pdgfrb_all_mat[rownames(pdgfrb_my_mat), colnames(pdgfrb_my_mat)] <- pdgfrb_my_mat


values <- log10(intact_all_mat/pdgfrb_all_mat)

ec_genes <- gene_ann[rownames(values),'lig']
ec_cells <- sapply(strsplit(colnames(values), ' -> '), '[', 1)
mes_genes_r1 <- gene_ann[rownames(values),'r1']
mes_genes_r2 <- gene_ann[rownames(values),'r2']
mes_genes_r2[mes_genes_r2==''] <- 'not_a_gene'
mes_cells <- sapply(strsplit(colnames(values), ' -> '), '[', 2)

if(!'not_a_gene' %in% rownames(avg_mes_inj)) {
  avg_mes_inj <- rbind(avg_mes_inj, 'not_a_gene' = rep(1, ncol(avg_mes_inj)))
}

expr_mask_ec <- avg_ec_inj[ec_genes, ec_cells]
expr_mask_mes_r1 <- avg_mes_inj[mes_genes_r1, mes_cells]
expr_mask_mes_r2 <- avg_mes_inj[mes_genes_r2, mes_cells]

# geni per cellula selezionati a diverse soglie di espressione
# (simile sia per EC che per MES)
###   0.3 ~ 4.5k
###   0.5 ~ 2.5k
###   1.0 ~ 1.5k

e_thr <- .5
e_idx <- apply(expr_mask_ec >= e_thr & 
                 expr_mask_mes_r1 >= e_thr & 
                 expr_mask_mes_r2 >= e_thr, 1, any)

values_filt <- values[e_idx,]
rn <- rownames(values_filt)

ann_names <- unique(ann_df[rn,])
ann_cols <- my_cols_palette[1:length(ann_names)]
names(ann_cols) <- ann_names
ann_cols <- list(pathway_name=ann_cols)

breaks <- seq(min(values_filt), max(values_filt), length.out = 101)
sym_id <- which.min((breaks + breaks[1])^2)
palette <- c(
  colorRampPalette(c('blue','lightblue','snow', 'pink','red'))(sym_id),
  colorRampPalette(c('red','darkred','black'))(101 - sym_id + 1)
)

pheatmap(values_filt, color = palette, breaks=breaks,
         annotation_row = ann_df[rn,,drop=F],
         annotation_colors = ann_cols,
         main = 'PDGFRb vs INTACT',
         filename = 'figures/InjuryD7_vs_intact_interaction.pdf', 
         height = 12, width = 8)