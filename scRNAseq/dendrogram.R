setwd("/beegfs/scratch/ric.cosr/ric.bonanomi")
library(Seurat)
library(ggtree)
library(tidyverse)
library(ggplot2)




integrated_3TIP_new = readRDS("integrated_3TIP_new.RDS")

levels(integrated_3TIP_new)

new.cluster.ids <- c('VEN_PLVAP+', 
                     'IMMAT',
                     'BARR',
                     'ARTERIAL',
                     'CAP_PLVAP-',
                     'PROLIF',
                     'VEN_PLVAP-',
                     'TIP',
                     'CAP_PLVAP+'
)


names(new.cluster.ids) <- levels(integrated_3TIP_new)
integrated_3TIP_new <- RenameIdents(integrated_3TIP_new, new.cluster.ids)

DefaultAssay(integrated_3TIP_new) = "integrated"


seurat <- BuildClusterTree(
  integrated_3TIP_new,
  dims = 1:25,
  reorder = FALSE,
  reorder.numeric = FALSE
)


tree <- seurat@tools$BuildClusterTree

dataset1<-data.frame("name" = c("VEN_PLVAP+","IMMATURE", "BARR_END_CAP","ARTERIAL","CAP_PLVAP-","PROLIF","VEN_PLVAP-","TIP","CAP_PLVAP+"), 
                     "colour" = c("#FF6666","#6600CC","#336666","#0066FF","#399933","#FF99CC","#990000","#FF00FF","#99CC33"))


p <- ggtree::ggtree(tree, aes(x, y)) +
  # scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(color=dataset1$colour, shape = 16, size = 8) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,3.5,1,0), 'cm'))


pdf("nice_dendrogram_integrated_EC_25PC.pdf")
p
dev.off()



#################################################################################
sample1 <- subset(x = integrated_3TIP_new, subset = stim == "1")


names(new.cluster.ids) <- levels(sample1)
sample1 <- RenameIdents(sample1, new.cluster.ids)

DefaultAssay(sample1) = "integrated"


seurat1 <- BuildClusterTree(
  sample1,
  dims = 1:25,
  reorder = FALSE,
  reorder.numeric = FALSE
)

tree1 <- seurat1@tools$BuildClusterTree


p <- ggtree::ggtree(tree1, aes(x, y)) +
  # scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(color=dataset1$colour, shape = 16, size = 8) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,3.5,1,0), 'cm'))

pdf("nice_dendrogram_Intact_EC_25PC.pdf")
p
dev.off()


###############################################################################
sample2 <- subset(x = integrated_3TIP_new, subset = stim == "2")


names(new.cluster.ids) <- levels(sample2)
sample2 <- RenameIdents(sample2, new.cluster.ids)

DefaultAssay(sample2) = "integrated"


seurat2 <- BuildClusterTree(
  sample2,
  dims = 1:25,
  reorder = FALSE,
  reorder.numeric = FALSE
)

tree2 <- seurat2@tools$BuildClusterTree


p <- ggtree::ggtree(tree2, aes(x, y)) +
  # scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(color=dataset1$colour, shape = 16, size = 8) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,3.5,1,0), 'cm'))

pdf("nice_dendrogram_InjuryD7_EC_25PC.pdf")
p
dev.off()











