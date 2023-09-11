# suppressMessages(library(ggbeeswarm))
suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
suppressMessages(library(gplots))
suppressMessages(library(ggplot2))
suppressMessages(library(openxlsx))
suppressMessages(library(cowplot))
suppressMessages(library(patchwork))

# nerve EC dataset

nerve_tissue <- readRDS('../../integrated_EC_3TIP.Rds')
nerve_tissue[['integrated']] <- NULL
nerve_tissue_uninj <- nerve_tissue[,nerve_tissue$stim == '1']
nerve_tissue_uninj$orig.ident <- 'nerve'


# following raw scRNA datasets were downloaed from
# https://endotheliomics.shinyapps.io/ec_atlas/:
# - raw_count_matrix_brain.txt.gz
# - raw_count_matrix_colon.txt.gz
# - raw_count_matrix_heart.txt.gz
# - raw_count_matrix_kidney.txt.gz
# - raw_count_matrix_liver.txt.gz
# - raw_count_matrix_lung.txt.gz
# - raw_count_matrix_muscle_EDL.txt.gz
# - raw_count_matrix_muscle_soleus.txt.gz
# - raw_count_matrix_small_intestine.txt.gz
# - raw_count_matrix_spleen.txt.gz
# - raw_count_matrix_testis.txt.gz

tissues_list <- lapply(
  list.files('.','txt.gz',full.names=T),
  function(f) {
    pname <- sub('\\.txt\\.gz','',
                 sub('raw_count_matrix_','',basename(f)))
    print(paste(f, pname))
    counts <- read.table(gzfile(f), header = T)
    rownames(counts) <- counts$Feature
    counts$Feature <- NULL
    CreateSeuratObject(
      counts = counts, 
      project = pname,
      min.cells = 3, min.features = 200)
  })
names(tissues_list) <- sapply(tissues_list, Project)



tissues_merged <- tissues_list[[1]]
for(i in 2:length(tissues_list)) 
  tissues_merged <- merge(tissues_merged, tissues_list[[i]])
tissues_merged



Idents(tissues_merged) <- 'all_tissues'
tissues_merged[["percent.mt"]] <- PercentageFeatureSet(tissues_merged, pattern = "^mt-")
VlnPlot(tissues_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# (i) genes expressed by less than 10 cells or with a row average of < 0.002 were not considered; 
# (ii) cells that expressed fewer than 300 genes (low quality), and cells that expressed over 4,000 genes (potential doublets) were excluded from further analysis;
# (iii) cells in which over 10% of unique molecular identifiers (UMIs) were derived from the mitochondrial genome were removed.


# For in silico EC selection, we first identified highly variable genes using the Seurat FindVariableGenes function (mean lower threshold = 0.0125, mean higher threshold = 8, dispersion threshold = 0.5. Data 

tissues_merged <- NormalizeData(tissues_merged)
tissues_merged <- FindVariableFeatures(tissues_merged, mean.cutoff = c(0.0125, 8), dispersion.cutoff = c(.5, Inf))
tissues_merged <- ScaleData(tissues_merged, features = rownames(tissues_merged))
tissues_merged <- RunPCA(tissues_merged, npcs = 20, features = VariableFeatures(object = tissues_merged))
ElbowPlot(tissues_merged, ndims = 20)



# (t-SNE, Rtsne package; top 8 principal components (PCs))
# FindClusters function in Seurat (clus- tering resolution = 1, k-nearest neighbors = 10). 

tissues_merged <- FindNeighbors(tissues_merged, dims = 1:10, k.param = 10)
tissues_merged <- FindClusters(tissues_merged, resolution = 1.2)
tissues_merged <- RunUMAP(tissues_merged, dims = 1:10)

# EC clusters were annotated based on the expression of known EC and 
# non-EC marker genes, including Pecam1 and Cdh5 (vascular ECs), 
# Prox1 and Lyve-1 (lymphatic ECs), Col1a1 (fibroblasts), 
# Hba-a1, Hba- a2, Hbb-bs (red blood cells) and 
# Acta2 (smooth muscle cells). 
# Contaminating cell clusters (non-ECs) were removed, and all 
# down- stream analysis was performed on ECs only.

# FeaturePlot(tissues_merged, c('Pecam1','Cdh5','Prox1','Lyve1','Col1a1','Acta2','Hba-a1','Hba-a2','Hbb-bs'), ncol=2)
# VlnPlot(tissues_merged, c('Pecam1','Cdh5','Prox1','Lyve1','Col1a1','Acta2','Hba-a1','Hba-a2','Hbb-bs'), stack = TRUE, flip=TRUE, same.y.lims=TRUE)
to_remove <- c(8 ,17, 21, 25, 28, 29, 34, 36)
to_keep <- setdiff(0:36,to_remove)
sum(table(tissues_merged$seurat_clusters)[as.character(to_remove)])
sum(table(tissues_merged$seurat_clusters)[as.character(to_keep)])

tissues_merged <- tissues_merged[
  ,tissues_merged$seurat_clusters %in% as.character(to_keep)]

# keep names of retained cells
filtered_ECcells <- tissues_merged$orig.ident

# For each tissue independently, we identified highly variable genes using the Seurat FindVariableGenes on the in silico selected ECs (mean lower threshold = 0.0125, mean higher threshold = 8, dispersion threshold = 0.5). 
# We performed PCA on highly variable genes, followed by t-SNE or Uniform Manifold Approximation and Projection (UMAP) (umap package) (McInnes et al., 2018) visualization (top 8 PCs). 
# We used an in-house developed tool to color-code t-SNE plots for all detected genes, to empirically define the clustering resolution for each tissue. 
# Next, we applied graph-based clustering using the FindClusters function in Seurat (clustering resolution = 1, k-nearest neighbors = 10) and verified by t-SNE visualization that all expected clusters were captured. 
# Clusters with highly similar expression patterns indicative to underlie the same EC phenotype were merged into the same cluster. 

filtered_ECcells_bytissue <- split(names(filtered_ECcells), filtered_ECcells)
tissues_list <- lapply(
  names(filtered_ECcells_bytissue),
  function(pname) {
    f <- paste0('raw_count_matrix_',pname,'.txt.gz')
    print(paste(f, pname))
    counts <- read.table(gzfile(f), header = T)
    rownames(counts) <- counts$Feature
    counts$Feature <- NULL
    CreateSeuratObject(
      counts = counts[,filtered_ECcells_bytissue[[pname]]], 
      project = pname)
  })
names(tissues_list) <- sapply(tissues_list, Project)
names(tissues_list)[names(tissues_list)=='muscle_EDL'] <- 'EDL'
names(tissues_list)[names(tissues_list)=='muscle_soleus'] <- 'soleus'
names(tissues_list)[names(tissues_list)=='small_intestine'] <- 'small intestine'
#
names(filtered_ECcells_bytissue)[names(filtered_ECcells_bytissue)=='muscle_EDL'] <- 'EDL'
names(filtered_ECcells_bytissue)[names(filtered_ECcells_bytissue)=='muscle_soleus'] <- 'soleus'
names(filtered_ECcells_bytissue)[names(filtered_ECcells_bytissue)=='small_intestine'] <- 'small intestine'

kalucka_pipeline <- function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, mean.cutoff = c(0.0125, 8), dispersion.cutoff = c(.5, Inf))
  x <- ScaleData(x, features = rownames(x))
  x <- RunPCA(x, npcs = 20, features = VariableFeatures(object = x))
  x <- FindNeighbors(x, dims = 1:10, k.param = 10)
  x <- FindClusters(x, resolution = 1)
  x <- RunUMAP(x, dims = 1:10)  
}
tissues_list <- lapply(tissues_list, kalucka_pipeline)

# downscale every experiment to 1k
tissues_list_subset <- lapply(tissues_list, function(x) {
  x[,sample(1:ncol(x), 999)]
})
# add nerve
tissues_list_subset$nerve <- nerve_tissue_uninj

# merge and analyse jointly

tissues_merged <- tissues_list_subset[[1]]
for(i in 2:length(tissues_list_subset)) 
  tissues_merged <- merge(tissues_merged, tissues_list_subset[[i]])
tissues_merged

tissues_merged <- NormalizeData(tissues_merged)
tissues_merged <- FindVariableFeatures(tissues_merged, mean.cutoff = c(0.0125, 8), dispersion.cutoff = c(.5, Inf))
tissues_merged <- ScaleData(tissues_merged, features = rownames(tissues_merged))
tissues_merged <- RunPCA(tissues_merged, npcs = 20, features = VariableFeatures(object = tissues_merged))
ElbowPlot(tissues_merged, ndims = 20)
tissues_merged <- FindNeighbors(tissues_merged, dims = 1:10, k.param = 10)
tissues_merged <- FindClusters(tissues_merged, resolution = 1.2)
tissues_merged <- RunUMAP(tissues_merged, dims = 1:10)



umap_merged_withnerve <- DimPlot(tissues_merged, 
                                 group.by = 'orig.ident', 
                                 label = TRUE)

pdf('umap_merged_withnerve.pdf', height = 5, width=6); 
umap_merged_withnerve; 
dev.off()


######## EC dendrogram

library(pheatmap)
library(pvclust)

Idents(tissues_merged) <-tissues_merged$orig.ident 
tissue_marks <- FindAllMarkers(tissues_merged, 
                               max.cells.per.ident = 1e3)

best_markers <- with(tissue_marks[order(tissue_marks$avg_log2FC),], sapply(split(gene, cluster), head, n=200))
str(best_marker_genes <- unique(c(best_markers)))
avgexpr_merged <- AverageExpression(tissues_merged, 
                                    group.by = 'orig.ident', slot='scale.data',
                                    features = best_marker_genes)

pvclust_merged_RNA <- pvclust(avgexpr_merged$RNA, 
                              method.hclust = 'complete', method.dist = 'correlation')
pvclust_merged_RNA_nonerve <- pvclust(avgexpr_merged$RNA[,colnames(avgexpr_merged$RNA) != "nerve"],
                                      method.hclust = 'complete', method.dist = 'correlation')


pdf('tissues_merged_1k_dendrogram_N20.pdf', height = 4, width = 3)
par(mar=c(1,1,1,1))
plot(pvclust_merged_RNA$hclust, hang = 0.01, frame.plot = F, ylab='', main ='', yaxt = 'n', 
     xlab='', sub='')
dev.off()

######## label transfer from nerve to kalucka

library(Seurat)
library(CelliD)

integrated_EC_3TIP <- readRDS('../../integrated_EC_3TIP.Rds')
Idents(integrated_EC_3TIP, 
       cells = colnames(integrated_EC_3TIP)[
         integrated_EC_3TIP$seurat_clusters == 8]) <- 
  'INTERMEDIATE_BEC'
palette <- c('#0066FF','#336666','#99CC33','#399933','#E0E0E0',
                      '#6600CC','#FF99CC','#FF00FF','#FF6666','#990000')
names(palette) <- sort(levels(Idents(integrated_EC_3TIP)))
integrated_EC_3TIP$cell_types <- factor(Idents(integrated_EC_3TIP), 
                                        levels = names(palette))
DefaultAssay(integrated_EC_3TIP) <- 'RNA'
integrated_EC_3TIP <- RunMCA(integrated_EC_3TIP)
EC_all_1_gs <- GetGroupGeneSet(subset(integrated_EC_3TIP, stim == 1), 
                               dims = 1:50, n.features = 200)


tissues_merged <- merge(tissues_list[[1]], tissues_list[-1])
tissues_merged <- FindNeighbors(tissues_merged, 
                                    dims = 1:30, k.param = 10)
tissues_merged <- RunUMAP(tissues_merged, 
                              dims = 1:30)
tissues_merged <- RunMCA(tissues_merged)
HGTkm_EC_all_1 <- RunCellHGT(tissues_merged, 
                             pathways = EC_all_1_gs, dims = 1:50)

mapHGT <- function(obj, HGT, name, conf=2) {
  HGT_prediction <- rownames(HGT)[apply(HGT, 2, which.max)]
  HGT_conf <- apply(HGT, 2, max)
  HGT_signif <- ifelse(HGT_conf>conf, yes = HGT_prediction, "unassigned")
  obj[[name]] <- HGT_signif
  obj[[paste0(name,'_conf')]] <- HGT_conf
  return(obj)
}
HGTtable <- function(obj, name, lev='orig.ident', norm=TRUE) {
  bp_data <- table(tissues_merged[[name]][,1], tissues_merged[[lev]][,1])
  if(norm) bp_data <- t(t(bp_data)/colSums(bp_data))
  return(bp_data)
}

tissues_merged <- mapHGT(tissues_merged, 
                             HGTkm_EC_all_1, 
                             'EC_all_1')
tab_EC_all_1 <- HGTtable(tissues_merged, 'EC_all_1')


pdf('kalucka_label_transfer.pdf', height = 4, width = 6)
par(mar=c(6,3,1,0)+.1)
barplot(100*tab_EC_all_1, col = c(palette, unassigned='black')[rownames(tab_EC_all_1)],
        angle = 90, cex.names = .7, las=2)
dev.off()
