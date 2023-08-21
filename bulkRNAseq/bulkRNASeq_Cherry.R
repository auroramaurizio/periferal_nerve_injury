---
  title: "Bonanomi_Bulk"
output: html_document
---

  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# set working diretory
prj = 'BonanomiD_1292_RNASeq'
PI = 'Bonanomi'

#setwd("/Users/maurizio.aurora/Documents/BonanomiD_1287_scRNA_injury_TdT/CHERRY_NEW")


# load libraries
suppressMessages(library("dplyr"))
suppressMessages(library("edgeR"))
suppressMessages(library("data.table"))
suppressMessages(library("DESeq2"))
suppressMessages(library("openxlsx"))
suppressMessages(library("dplyr"))
suppressMessages(library("ggplot2"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("enrichR"))
suppressMessages(library("grid"))
suppressMessages(library("gridExtra"))
suppressMessages(library('VennDiagram'))
suppressMessages(library("venn"))
suppressMessages(library("cowplot"))
suppressMessages(library("ggpubr"))
suppressMessages(library("viridis"))
suppressMessages(library('pals'))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("wesanderson"))
suppressMessages(library("patchwork"))
suppressMessages(library("magick"))
suppressMessages(library('philentropy'))
suppressMessages(library('IntClust'))
suppressMessages(library('pheatmap'))
suppressMessages(library("assertr"))
suppressMessages(library("remotes"))
suppressMessages(library("GeneOverlap"))
suppressMessages(library('stringr'))
suppressMessages(library("ggrepel"))
#install.packages("multcompView", repos="http://R-Forge.R-project.org")
suppressMessages(library("multcompView"))
suppressMessages(library("sva"))
suppressMessages(library("ComplexHeatmap"))


# path to the count files (2 seq runs containg both Cherry and Tomato samples)
filecount_1 = '/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/raw_not_corrected_counts_tables_1133_1292/1292/all.counts.gz'
filecount_2='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/raw_not_corrected_counts_tables_1133_1292/1133/all.counts.gz'

#define the number of replicates
Nreplica= 3

# import metadata
metadata = read.xlsx("/Users/maurizio.aurora/Documents/metadata.xlsx")
#remove outliers from metadata
metadata <- subset(metadata, sample!= c("72","DB_2_130"))

# import counts
annotation <- c('GeneID','Chr','Start','End','Strand','Length')
fCounts_1 <- read.delim(file=filecount_1, header=TRUE, check.names = F)
fCounts_2 <- read.delim(file=filecount_2, header=TRUE, check.names = F)

fCounts <- merge(fCounts_1, fCounts_2, by=c("Geneid","Chr","Start","End","Length","Strand"))


fCountsData_not_adj <- fCounts[
  ,
  -which(
    tolower(names(fCounts))
    %in%
      tolower(annotation))]


fCountsAnnotation <- fCounts[
  ,
  which(
    tolower(names(fCounts))
    %in%
      tolower(annotation))]


fCounts_bis= select(fCounts, "Geneid","Chr","Start","End","Length","Strand" )

geneidColname <- 'Geneid'
geneidIdx <- which(tolower(annotation) %in% tolower(geneidColname))
rownames(fCountsData_not_adj) <- fCounts[[geneidIdx]]

# Reordering counts matrix to have samples ordered as in metadata
# it assumes that the first column in metadata is sample name
fCountsData_not_adj <- fCountsData_not_adj[,match(metadata[,1], colnames(fCountsData_not_adj))]
rownames(fCounts_bis) = fCounts_bis$Geneid
fCountsData_not_adj_temp = fCountsData_not_adj
fCountsData_not_adj_temp$Geneid = rownames(fCountsData_not_adj)

fCounts_cherry_tomato <- merge(fCounts_bis, fCountsData_not_adj_temp, by="Geneid")

fCounts_onlycherry = fCounts_cherry_tomato[,c("Geneid","Chr","Start","End","Length","Strand","DB_1_127","DB_5_148","DB_6_151","101","104","74")]

fCounts_onlytomato = fCounts_cherry_tomato[,c("Geneid","Chr","Start","End","Length","Strand","DB_3_131","DB_4_133", "DB_7_152", "DB_8_153",
                                                            "DB_9_154","DB_10_189","DB_11_190", "DB_12_191", "77", "80", "83", "92", "95", "98")]


write.table(fCounts_onlycherry, "Raw_counts_Cherry.tsv", sep = "\t", row.names = T, col.names = T )
write.table(fCounts_onlytomato, "Raw_counts_Tomato.tsv", sep = "\t", row.names = T, col.names = T )




#####################
# ONLY CHERRY SAMPLES
# 7DPI
# KO vs WT
#####################


metadata_d = metadata[metadata$reporter == 'KDR_cherry',]
#removed samples 72 and 130

getwd()

row.names(metadata_d) = metadata_d$sample
cov1=as.factor(metadata_d$nerve)
cov2=factor(metadata_d$Plxnd1, levels = c('WT', 'KO'))
cov3=as.factor(metadata_d$reporter)
exp_batch=as.factor(metadata_d$exp)
covar_mat <- cbind(cov2, cov3)

fCountsData_notadjusted_p = fCountsData_not_adj[,row.names(metadata_d)]

# correct batch effect with ComBatSeq
# Using DEseq2 covariates in exp design formula is not enough.

fCountsData_adjusted_d <- ComBat_seq(as.matrix(fCountsData_notadjusted_p), batch=exp_batch, group=cov2 ) #combined


# save tables about Cherry exp in multisheet Excel file.

SAVE_variable <- list()
filename_xls <- paste('COUNTS_cherry',prj,'.xlsx', sep='')
variable2save_names <- c('all_counts', 'expGenes_counts','expGenes_LogCPM', 'expGenes_LogRPKM', 'expGenes_RPKM')

# all_counts
y <- DGEList(counts=fCountsData_adjusted_d, genes = fCountsAnnotation)
SAVE_variable[[variable2save_names[1]]] <- as.data.frame(y$counts)

# expGENES_counts
keep <- rowSums(cpm(y)>1)>=Nreplica
table(keep)
yf <- y[keep,]
SAVE_variable[[variable2save_names[2]]] <- as.data.frame(yf$counts)

#CPM
SAVE_variable[[variable2save_names[3]]] <- as.data.frame(cpm(yf, log=T))

#RPKM log
SAVE_variable[[variable2save_names[4]]] <- as.data.frame(rpkm(yf, log=T, gene.length =yf$genes$Length))

#RPKM not log
SAVE_variable[[variable2save_names[5]]] <- as.data.frame(rpkm(yf, log=F, gene.length =yf$genes$Length))

SAVE_variable <- list()
prj = "logRpkm_all_genes_cherry"
filename_xls <- paste('COUNTS_combatSeq',prj,'.xlsx', sep='')
variable2save_names <- c('all_genes_log_rpkm')

#RPKM log all genes
SAVE_variable[[variable2save_names[1]]] <- as.data.frame(rpkm(y, log=T, gene.length =y$genes$Length))

write.xlsx(SAVE_variable,
           file = filename_xls,
           row.names = T,
           asTable = T,
           sheetName =variable2save_names)



fCountsData_adjusted_df = as.data.frame(fCountsData_adjusted_d)
fCountsData_adjusted_df$Geneid = rownames(fCountsData_adjusted_df)
fCountsData_adjusted_df =  merge(fCountsData_adjusted_df, fCounts_bis, by="Geneid")
fCountsData_adjusted_df2 <- fCountsData_adjusted_df[, c("Geneid", "Chr","Start","End","Length","Strand","DB_1_127","DB_5_148","DB_6_151","101","104","74")]

write.table(fCountsData_adjusted_d, "Counts_Cherry_adjusted_Combat.tsv", sep = "\t", row.names = T, col.names = T )
write.table(fCountsData_adjusted_df2, "Raw_counts_Cherry_corrected_ComBat_Seq.tsv", sep = "\t", row.names = T, col.names = T )

y <- DGEList(counts=fCountsData_adjusted_d, genes = fCountsAnnotation)
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$Length)

provaA=rowSums(cpm(y)>1)>=Nreplica
head(provaA)
keep <- rowSums(cpm(y)>1)>=Nreplica
length(provaA)
yf <- y[keep,]
nrow(yf)
N=500
vary <- apply(fCountsRPKM[keep,],1,var)
vary_s <- sort(vary, decreasing = T)
TOP_N <- names(vary_s[1:N])
yTOP <-  y[TOP_N,]
fCountsRPKMTOP <- fCountsRPKM[TOP_N,]
#PCA parameters
pcx = 1
pcy = 2
centering = TRUE
scaling = TRUE
# PCA
pca = prcomp(t(fCountsRPKMTOP), center=centering, scale=scaling)
var = round(matrix(((pca$sdev^2)/(sum(pca$sdev^2))), ncol=1)*100,1)
score = as.data.frame(pca$x)
# plot paramters
xlab = paste("PC", pcx, " (",var[pcx],"%)", sep="")
ylab = paste("PC", pcy, " (",var[pcy],"%)", sep="")
cum = var[pcx]+var[pcy]
names = rownames(pca$x)

score$exp = metadata_d$exp
score$reporter = metadata_d$reporter
score$nerve = metadata_d$nerve
score$Plxnd1 = metadata_d$Plxnd1
score$sampleID = metadata_d$sample


#top 500 most variable genes PCA

pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy],
                        color=Plxnd1, shape = exp))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy],
                                    color=Plxnd1, shape=exp, label = str_sub(x <- sampleID,-5,-1)),
                   size = 5,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
                   segment.color = 'grey50') +
  geom_point(size= 7)+
  labs(x=xlab, y=ylab, title=paste("PC",pcx," vs PC",pcy," scoreplot",sep="")) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
  geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
  theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
        axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5),
        axis.title.x = element_text(face = "bold", color = "black", size = 24),
        axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
        axis.title.y = element_text(face = "bold", color = "black", size = 24),
        legend.text = element_text(face = "bold", color = "black", size = 18),
        legend.position="right",
        panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
  scale_color_manual(values = brewer.pal(n = 4, name = "Spectral"))

options(repr.plot.width=13, repr.plot.height=7)
pdf('PCA_top500rpkm_Cherry_test_geo.pdf', width = 14, height = 10)
pca
dev.off()

annotation_column <- metadata_d[,2:(dim(metadata_d)[2])]

score$reporter = metadata_d$reporter
score$nerve = metadata_d$nerve
score$Plxnd1 = metadata_d$Plxnd1
score$sampleID = metadata_d$sample

row.names(annotation_column) <- metadata_d[,1]
annotation_column = annotation_column[,c('nerve','Plxnd1')]
options(repr.plot.width=12, repr.plot.height=10)


# top 500 most variable genes heatmap

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)

HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         cluster_rows = T,
                         cluster_cols = T,
                         show_rownames = F,
                         show_colnames = F,
                         width = 5,
                         height = 5,
                         col= colors,
                         treeheight_col  = 0,
                         treeheight_row  = 0,
                         annotation_colors = ann_colors,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                         display_numbers = F,
                         filename = 'Heatmap_500rpkm_Cherry_white.pdf')


# filtered genes heatmap

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
ann_colors = list(nerve = c( 'crush_d7' = 'purple'),
                  Plxnd1 = c('WT' = '#B1624EFF','KO' = '#5CC8D7FF'))

toMatch = readLines("/Users/maurizio.aurora/Documents/updw_annotation_cherry.txt")
subset = readLines("/Users/maurizio.aurora/Documents/heatmap_cherry.txt")
filtered = fCountsRPKM[rownames(fCountsRPKM) %in% subset ,]
refiltered <- filtered[match(subset, rownames(filtered)), ]

refiltered <- refiltered[,c("DB_1_127","DB_5_148","101","DB_6_151","104","74")]
annotation_column_ = annotation_column[c("DB_1_127","DB_5_148","101","DB_6_151","104","74"),]

matches <- grep(paste(toMatch,collapse="|"),
                rownames(refiltered))

matched <- grep(paste(toMatch,collapse="|"),
                rownames(refiltered), value = T)

ha = rowAnnotation(foo = anno_mark(at = matches, labels = matched))

# FigS3 A

pdf('Heatmap_logrpkm_Cherry_white_filtered_annotated_smallish_test.pdf',10, 12)
ComplexHeatmap::pheatmap(refiltered,
                         scale = 'row',
                         annotation_col = annotation_column_,
                         cluster_rows = F,
                         cluster_cols = F,
                         show_rownames = F,
                         show_colnames = F,
                         #height = 11,
                         #width = 10,
                         col= colors,
                         treeheight_col  = 0,
                         treeheight_row  = 0,
                         annotation_colors = ann_colors,
                         fontsize = 18, fontsize_row = 15, fontsize_col = 14,
                         display_numbers = F,
                         right_annotation = ha)
dev.off()

#####MiniHeatmaps#######
# Fig3 F
annotation_column = as.data.frame(colnames(KOWTac))
colnames(annotation_column)= "Plxnd1"
rownames(annotation_column) <- colnames(KOWTac)
annotation_column$nerve = "crush_d7"
head(annotation_column)


comparison = list(KO = c("KO"),
                  WT = c("WT"))


metadata_crush_KO = metadata_d[metadata_d$Plxnd1 %in% comparison[[1]],]
KO = as.data.frame(fCountsRPKM[,row.names(metadata_crush_KO)])
KO$KO = rowMeans(KO)


metadata_crush_WT = metadata_d[metadata_d$Plxnd1 %in% comparison[[2]],]
WT = as.data.frame(fCountsRPKM[,row.names(metadata_crush_WT)])
WT$WT = rowMeans(WT)


angio_cherry = readLines("/Users/maurizio.aurora/Documents/angio_cherry")
interferon_cherry = readLines("/Users/maurizio.aurora/Documents/interferon_cherry")
TNFalpha_cherry = readLines("/Users/maurizio.aurora/Documents/TNFalpha_cherry")
TGFbeta_cherry = readLines("/Users/maurizio.aurora/Documents/TGFbeta_cherry")
inflammation_cherry = readLines("/Users/maurizio.aurora/Documents/inflammation_cherry")


WTac = WT[rownames(WT) %in% angio_cherry,]
KOac = KO[rownames(KO) %in% angio_cherry,]


KOWTac = cbind(KOac, WTac)
head(KOWTac)

KOWTac <- KOWTac[,c("DB_1_127","DB_5_148","101","DB_6_151","104","74")]


# Fig3 F
pdf('MiniHeatmap_Angiogenesis_Cherry_not_avg.pdf')
pheatmap(as.matrix(KOWTac),
         main = "Angiogenesis",
         show_rownames = T,
         show_colnames = F,
         cellheight = 35,
         cellwidth = 35,
         cluster_rows = F,
         cluster_cols = F,
         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                               "RdYlBu")))(100),
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column_)
dev.off()

# Fig3 F
pdf('MiniHeatmap_Angiogenesis_Cherry.pdf')
pheatmap(as.matrix(KOWTac),
         main = "Angiogenesis",
         show_rownames = T,
         show_colnames = F,
         cellheight = 45,
         cellwidth = 45,
         cluster_rows = F,
         cluster_cols = F,
         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                               "RdYlBu")))(100),
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()

WTac = WT[rownames(WT) %in% interferon_cherry,]
KOac = KO[rownames(KO) %in% interferon_cherry,]

KOWTac = cbind(KOac, WTac)
head(KOWTac)

KOWTac <- KOWTac[,c("DB_1_127","DB_5_148","101","DB_6_151","104","74")]

# Fig3 F
pdf('MiniHeatmap_Interferon_Cherry_not_avg.pdf')
pheatmap(as.matrix(KOWTac),
         main = "Interferon",
         show_rownames = T,
         show_colnames = F,
         cellheight = 35,
         cellwidth = 35,
         cluster_rows = F,
         cluster_cols = F,
         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                               "RdYlBu")))(100),
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column_)
dev.off()

WTac = WT[rownames(WT) %in% TNFalpha_cherry,]
KOac = KO[rownames(KO) %in% TNFalpha_cherry,]


KOWTac = cbind(KOac, WTac)
head(KOWTac)

KOWTac <- KOWTac[,c("DB_1_127","DB_5_148","101","DB_6_151","104","74")]

# Fig3 F
pdf('MiniHeatmap_TNFalpha_Cherry_not_avg.pdf')
pheatmap(as.matrix(KOWTac),
         main = "TNFalpha",
         show_rownames = T,
         show_colnames = F,
         cellheight = 35,
         cellwidth = 35,
         cluster_rows = F,
         cluster_cols = F,
         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                               "RdYlBu")))(100),
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column_)
dev.off()


WTac = WT[rownames(WT) %in% inflammation_cherry,]
KOac = KO[rownames(KO) %in% inflammation_cherry,]

KOWTac = cbind(KOac, WTac)

KOWTac <- KOWTac[,c("DB_1_127","DB_5_148","101","DB_6_151","104","74")]

# Fig3 F
pdf('MiniHeatmap_inflammation_Cherry_not_avg.pdf')
pheatmap(as.matrix(KOWTac),
         main = "Inflammation",
         show_rownames = T,
         show_colnames = F,
         cellheight = 35,
         cellwidth = 35,
         cluster_rows = F,
         cluster_cols = F,
         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                               "RdYlBu")))(100),
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column_)
dev.off()



WTac = WT[rownames(WT) %in% TGFbeta_cherry,]
KOac = KO[rownames(KO) %in% TGFbeta_cherry,]


KOWTac = cbind(KOac, WTac)
head(KOWTac)

KOWTac <- KOWTac[,c("DB_1_127","DB_5_148","101","DB_6_151","104","74")]

# Fig3 F

pdf('MiniHeatmap_TGFbeta_Cherry_not_avg.pdf')
pheatmap(as.matrix(KOWTac),
         main = "TGFbeta",
         show_rownames = T,
         show_colnames = F,
         cellheight = 35,
         cellwidth = 35,
         cluster_rows = F,
         cluster_cols = F,
         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                               "RdYlBu")))(100),
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column_)
dev.off()




f = "Plxnd1"
seqc_pvalue = 0.01
dgeResults = list()

metadata_d[metadata_d$Plxnd1 %in% comparison[[1]],]
fCountsData_adjusted_s = fCountsData_adjusted_d[,row.names(metadata_d)]
colnames(fCountsData_adjusted_d)
row.names(metadata_d) = metadata_d$sample

comparison = list(KO_vs_WT = c("KO","WT"))

for (i in names(comparison)) {

  metadata_s = metadata_d[metadata_d$Plxnd1 %in% comparison[[i]],]
  fCountsData_adjusted_s = fCountsData_adjusted_d[,row.names(metadata_s)]
  dds <- DESeqDataSetFromMatrix(
    countData = fCountsData_adjusted_s,
    colData  = metadata_s,
    #design   = as.formula('~ exp + Plxnd1'))
    design   = as.formula('~Plxnd1'))
  filter <- rowSums(cpm(counts(dds)) >= 1) >= Nreplica
  table(filter)
  ddsFiltered <- dds[filter,]
  dga <- DESeq(
    object = ddsFiltered,
    test = "Wald",
    fitType = "parametric",
    betaPrior = FALSE,
    minReplicatesForReplace = Inf)

  alpha = 0.05
  print(paste(comparison[[i]][1],"_vs_",comparison[[i]][2],sep=''))
  dgeResults.tmp <- results(dga,
                            contrast             = c(f,comparison[[i]][1],comparison[[i]][2]),
                            cooksCutoff          = Inf,
                            independentFiltering = TRUE,
                            alpha                = alpha,
                            pAdjustMethod        = "BH")
  summary(dgeResults.tmp)
  head(dgeResults.tmp$pvalue)
  dgeResults[[i]] <- dgeResults.tmp[order(dgeResults.tmp$pvalue, decreasing = F),]

  #PCA
  vsd <- vst(dga, blind=FALSE)
  main.factor = "Plxnd1"



  pcaData <- plotPCA(vsd, intgroup=c(main.factor),returnData=TRUE)
  pcaData$exp=metadata_s$exp
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  #ggplot(pcaData, aes(PC1, PC2, color=Plxnd1, shape=exp))
  PCA = ggplot(pcaData, aes(PC1, PC2, color=Plxnd1, shape=exp)) +
    geom_label_repel(data= pcaData, aes(PC1, PC2, color=Plxnd1, label = str_sub(name,-5,-1)),
                     size = 6,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
                     segment.color = 'grey50') +
    geom_point(size=6) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    ggtitle(paste("PCA",i)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 12),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    scale_color_manual(values = c('salmon','grey'))
  pdf(paste('pca_',i,'.pdf',sep=''),width=8, height=6)
  print(PCA)
  dev.off()

  seqcUP = row.names(dgeResults[[i]])[dgeResults[[i]]$pvalue <= seqc_pvalue &
                                        !is.na(dgeResults[[i]]$padj)&
                                        dgeResults[[i]]$log2FoldChange > 1]

  if (length(seqcUP) > 0) {
    # print heatmap
    annotation_column <- metadata_s[,2:(dim(metadata_s)[2])]
    row.names(annotation_column) <- metadata_s[,1]
    options(repr.plot.width=12, repr.plot.height=10)
    crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
    colors = crp(255)
    head(cpm(counts(dga)))
    length(seqcUP)
    HP <- pheatmap::pheatmap(cpm(counts(dga))[seqcUP,],
                             scale = 'row',
                             annotation_col = annotation_column,
                             annotation_colors = ann_colors,
                             cluster_rows = T,
                             cluster_cols = T,
                             show_rownames = F,
                             cutree_cols = 2,
                             fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                             display_numbers = F,
                             col=colors,
                             filename = paste('Heatmap_seqcUP',i,'.pdf'),
                             width = 10, height = 11 )
  }

}

getwd()
f = 'DGE_results_Cherry'
dir.create(f, showWarnings=TRUE, recursive=TRUE)

lapply(
  names(dgeResults),
  function(x) write.table(
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname),
    file.path(f, paste(x, ".tsv", sep="")),
    append=F,
    row.names=F,
    col.names=T,
    quote=F,
    sep="\t"))


dgeResults_table = list()
dgeResults_table = lapply(
  names(dgeResults),
  function(x)
    data.table(
      data.frame(dgeResults[[x]]),
      keep.rownames=geneidColname))

names(dgeResults_table) = names(dgeResults)


write.xlsx(dgeResults_table,
           file = 'DGE_results_Cherry.xlsx',
           row.names = F,
           asTable = T,
           startRow = 1,
           sheetName = str_sub(names(dgeResults),1,31))






#plots
n.label = 20
FDR = T
pvalue = 0.01
for (Plxnd1 in names(dgeResults)) {
  results = as.data.frame(dgeResults[[Plxnd1]])
  results$DE = 'unm'
  if (!FDR) {
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'up'}
    if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
      results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < .1,]$DE = 'down' }
  } else {
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'up'}
    if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
      results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'down' }
  }
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange > 1,]$DE = 'SEQCup'}
  if (length(rownames(results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]))>0) {
    results[results$pvalue < pvalue & !is.na(results$padj) & results$log2FoldChange < -1,]$DE = 'SEQCdown' }
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange > 0,]$DE = 'FDRup'}
  if (length(rownames(results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]))>0) {
    results[results$padj < alpha & !is.na(results$padj) & results$log2FoldChange < 0,]$DE = 'FDRdown' }

  results$DE <- factor(x = results$DE, levels = c("unm", "FDRdown","FDRup", 'SEQCdown','SEQCup'))


  mycolors = c('grey','dodgerblue4','darkred','dodgerblue2','coral'); names(mycolors) = levels(results$DE)
  results$DE2 = 'unm'; results[results$DE!='unm',]$DE2 = 'mod'
  results$DE2 <- factor(x = results$DE2, levels = c("mod","unm"))
  mysize = c(3,2); names(mysize) = levels(results$DE2)
  myalpha = c(1,0.2); names(mysize) = unique(results$DE2)

  # label N genes
  N = min(n.label, length(rownames(results[results$DE == 'FDRup',])))
  up_label = rownames(results[results$DE == 'FDRup',])[1:N]
  N = min(n.label, length(rownames(results[results$DE == 'FDRdown',])))
  down_label = rownames(results[results$DE == 'FDRdown',])[1:N]

}



  results1 = results
  pvaluenew = 0.05
  results1$DE1 <- "NO"
  # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
  results1$DE1[results1$log2FoldChange > 0.58 & results1$pvalue < 0.05] <- "UP"
  # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
  results1$DE1[results1$log2FoldChange < -0.58 & results1$pvalue < 0.05] <- "DOWN"



  MAplot = ggplot(results) +
    geom_point(aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    xlim(c(0,1.e5)) +
    scale_x_continuous(trans='log10') +
    ggtitle(paste("MAPlot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "Mean expression", y = "log2 fold change")

  print(MAplot)
  pdf(paste(f,'/','MAplot_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(MAplot)
  dev.off()

  # Vulcano plot
  VP = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[down_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results[down_label,]), size = 3,
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")

  pdf(paste('VulcanoPlot_Dw_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(VP)
  dev.off()

  # Vulcano plot
  VP0 = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[up_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results[up_label,]), size = 3,
                     box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")

  pdf(paste('VulcanoPlot_Up_',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(VP0)
  dev.off()



  VP1 = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =2) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")

  print(VP1)
  pdf('VolcanoPlot_FDR_SEQC.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()

  VP1 = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =2) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")  +
    xlim(-5, 5)+
    ylim(0, 150)

  print(VP1)
  pdf('VolcanoPlot_FDR_SEQC_zoom.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()

  levels(results1$DE) = c("NO","DOWN","UP","NO","NO")
  results1$DE[results$DE == "SEQCdown"] <-"NO"
  results1$DE[results$DE == "SEQCup"] <- "NO"
  results1$DE[results$DE == "FDRup"] <- "UP"
  results1$DE[results$DE == "FDRdown"] <- "DOWN"
  levels(results1$DE) = c("NO","DOWN","UP","NO","NO")

  head(results1)

  mynewcolors = c('grey','dodgerblue4','darkred'); names(mynewcolors) = levels(results1$DE)

  unique(results1$DE)
  mynewcolors

  results1[grepl("Plxnd1", rownames(results1),]

  table((results1$DE))
  # Vulcano plot
  VP1 = ggplot(results1) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results1, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mynewcolors) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results1[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results1[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")

  print(VP1)
  pdf('VolcanoPlot_new_onlyFDR_pvalue0.05___.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()


  VP1 = ggplot(results1) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE), size =3) +
    geom_point(data = subset(results1, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE), size =3) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mynewcolors) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results1[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results1[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")  +
    xlim(-5, 5)+
    ylim(0, 150)

  print(VP1)
  pdf('VolcanoPlot_new_onlyFDR_pvalue0.05_zoom__.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()

  levels(results$DE)
  unique(results$DE)
  unique(results1$DE)

  results1 = results
  levels(results1$DE)
  levels(results1$DE) = c("NO","NO","NO","DOWN","UP")
  results1$DE[results$DE == "SEQCdown"] <-"DOWN"
  results1$DE[results$DE == "SEQCup"] <- "UP"
  results1$DE[results$DE == "FDRup"] <- "NO"

  results1$DE
  levels(results1$DE) = c("NO","NO","NO","DOWN","UP")


  # Vulcano plot
  VP1 = ggplot(results1) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results1, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mynewcolors) +
    #scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results1[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results1[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")

  print(VP1)
  pdf('VolcanoPlot_new_onlySEQC_pvalue0.01_log2FC1.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()


  tail(results1, 30)


  head(results1$DE)
  VP1 = ggplot(results1) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE), size =3) +
    geom_point(data = subset(results1, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE), size =3) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mynewcolors) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results1[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results1[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")  +
    xlim(-5, 5)+
    ylim(0, 150)

  print(VP1)
  pdf('VolcanoPlot_new_onlySEQC_pvalue0.01_log2FC1_zoom.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()





  # Vulcano plot
  VP1 = ggplot(results1) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE1), size =2) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mynewcolors) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results1[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE1),
                     label = row.names(results1[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")

  print(VP1)
  pdf('VolcanoPlot_new_pvalue0.05_log2FC0.58.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()


  tail(results1, 30)


  head(results1$DE)
  VP1 = ggplot(results1) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE1), size =2) +
    #geom_point(data = subset(results1, DE2 == 'mod'),
    #           aes(x=log2FoldChange, y=-log10(pvalue), color = DE), size =3) +
    ggtitle(paste("Volcano Plot,", Plxnd1)) +
    scale_color_manual(values = mynewcolors) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results1[genes_dario,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE1),
                     label = row.names(results1[genes_dario, ]), size = 3,
                     box.padding = unit(0.5, "lines"), max.overlaps = Inf, point.padding = unit(0.5, "lines"),segment.color = 'grey50') +
    theme(plot.title = element_text(color="black", size=10, face="bold.italic"),
          axis.text.x = element_text(angle = 90, face = "bold", color = "black", size=16, hjust =1),
          axis.title.x = element_text(face = "bold", color = "black", size = 16),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=16),
          axis.title.y = element_text(face = "bold", color = "black", size = 16),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.title = element_text(face = "bold", color = "black", size = 0),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid")) +
    labs(x = "log2 fold change", y = "-log10 p-value")  +
    xlim(-5, 5)+
    ylim(0, 150)

  print(VP1)
  pdf('VolcanoPlot_new_pvalue0.05_log2FC0.58_zoom.pdf',width=8, height=6.5)
  plot(VP1)
  dev.off()

  
  
  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +
          plot_layout(guides = "collect"))

  pdf(paste(f,'/','MA_VP_',Plxnd1,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +
          plot_layout(guides = "collect"))
  dev.off()



}



####### FDR ########

fdrUP = list()

alpha = 0.05
fdrUP = lapply(names(dgeResults),
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha &
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange > 0])
names(fdrUP)= names(dgeResults)

fdrDW = list()
fdrDW = lapply(names(dgeResults),
               function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$padj <= alpha &
                                                        !is.na(dgeResults[[x]]$padj)&
                                                        dgeResults[[x]]$log2FoldChange < 0])
names(fdrDW)= names(dgeResults)

####### SEQC ########
seqcUP = list()

pvalue = 0.01
seqcUP = lapply(names(dgeResults),
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue &
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange > 1])
names(seqcUP)= names(dgeResults)

seqcDW = list()
seqcDW = lapply(names(dgeResults),
                function(x) row.names(dgeResults[[x]])[dgeResults[[x]]$pvalue <= pvalue &
                                                         !is.na(dgeResults[[x]]$padj)&
                                                         dgeResults[[x]]$log2FoldChange < -1])
names(seqcDW)= names(dgeResults)

print('FDRup');lengths(fdrUP); print('FDRdw');lengths(fdrDW)
print('SEQCup');lengths(seqcUP); print('SEQCdw');lengths(seqcDW)

print('FDR_tot');lengths(fdrUP) + lengths(fdrDW)
print('SEQC_tot'); lengths(seqcUP) + lengths(seqcDW)


###########################
### ENRICHR enrichment ####
###########################


databases <- listEnrichrDbs()

dir.create('enrichR/', showWarnings=FALSE, recursive=TRUE)
# -------------------------
# enrichment Parameters
# -------------------------
# databases to make the enrichment of
enrich.databases <- c("GO_Biological_Process_2018",
                      "GO_Cellular_Component_2018",
                      "GO_Molecular_Function_2018",
                      "Reactome_2016",
                      "KEGG_2016",
                      "WikiPathways_2016",
                      "BioCarta_2016",
                      "Jensen_TISSUES",
                      "Jensen_COMPARTMENTS",
                      "Jensen_DISEASES")
# alpha used in DGE
padj.cutoff = alpha;

dgeResults
# -------------------------
# Perform Enrichment
# -------------------------
enrichr.list <- list()
for (i in 1:length(dgeResults)){
  print(names(dgeResults)[i])
  .res <- dgeResults[[i]]
  up.genes   <- fdrUP[[i]]
  down.genes <- fdrDW[[i]]
  both.genes <- c(up.genes, down.genes)
  write.table(up.genes, paste('./enrichR/FDRup_',names(dgeResults)[i],
                              '.txt', sep =''), quote = F,
              row.names = F, col.names = F)
  write.table(down.genes, paste('./enrichR/FDRdw_',names(dgeResults)[i],
                                '.txt', sep =''),
              quote = F, row.names = F, col.names = F)
  write.table(both.genes, paste('./enrichR/FDRboth_',names(dgeResults)[i],
                                '.txt', sep =''), quote = F,
              row.names = F, col.names = F)


  enrichr.list[[i]] <- lapply(list(up.genes,down.genes,both.genes),function(x) {
    enrichR::enrichr(genes = x, databases = enrich.databases)
  })
  names(enrichr.list[[i]]) <-  c("fdr_up","fdr_down","fdr_both")
}
names(enrichr.list) <- names(dgeResults)

# -----------------------------
# Write excels files
# -----------------------------


for (i in 1:length(dgeResults)){
  for (j in c("fdr_up","fdr_down","fdr_both")){
    filename = paste(
      file.path('enrichR',
                names(dgeResults)[[i]]),
      j,
      ".xlsx",
      sep="_")
    write.xlsx(x = enrichr.list[[i]][[j]], file = filename)
  }
}

write.xlsx(x = metadata, file= "metadata.xlsx")




###### enrichment #######

require("biomaRt")

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  return(genesV2)
}

# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  return(genesV2)
}


require("biomaRt")

# GSEA list
outdir = 'GSEA/'
dir.create(outdir)
for (c in names(dgeResults)) {
  l_ranked = data.frame(GeneID= row.names(dgeResults[[c]]), LogFC = dgeResults[[c]][,'log2FoldChange'])
  l_ranked = l_ranked[order(l_ranked$LogFC, decreasing = T),]
  musGenes <-  l_ranked$GeneID
  converted=convertMouseGeneList(musGenes)
  result = merge(l_ranked, converted, by.x = "GeneID", by.y = "MGI.symbol")
  result1 <- result[c(3,2)]
  colnames(result1)<- c("GeneID","LogFC")
  result1 = result1[order(result1$LogFC, decreasing = T),]
  write.table(result1, file = paste(outdir, c,'_ranked_list.rnk', sep =''),
              quote = F, row.names= F, col.names = F, sep ='\t')

}



#########################

### Jaccard distances ###
# -------------------------
# Jaccard distances
# -------------------------
dir = paste("./enrichR/JaccardPlots/", sep = '')
dir.create(dir)
#'WikiPathways_2016'
#database = c('KEGG_2016','Reactome_2016','GO_Biological_Process_2018')
p.value.thr = 0.01

breaksList = seq(0, 1, by = 0.001)

#cutree_rows_values = c(5,12,3,10,5,1)



lf = list.files(paste("./enrichR/", sep=''), pattern=glob2rx("*.xlsx"))
print(lf)
for (file in lf) {
  enrichR.file = paste("./enrichR/",file,sep='')
  s = openxlsx::getSheetNames(enrichR.file)
  Pathways.Table = data.frame()
  for (dat in s[1:length(s)]) {
    Table <- read.xlsx(xlsxFile = enrichR.file,
                       sheet = dat,
                       startRow = 1,
                       colNames = TRUE,
                       rowNames = TRUE,
                       detectDates = FALSE,
                       skipEmptyRows = TRUE,
                       skipEmptyCols = TRUE,
                       na.strings = "NA",
                       fillMergedCells = FALSE)

    Pathways.Table = rbind(Pathways.Table, Table)
  }

  gene.list = list()
  gene.all = character()
  pathways = unlist(row.names(Pathways.Table[Pathways.Table$Adjusted.P.value < p.value.thr,]))
  length(pathways)
  for (p in pathways) {
    gene.list[[p]] <- unlist(strsplit(Pathways.Table[p,]$Genes, ';'))
    gene.all = c(gene.all, unlist(strsplit(Pathways.Table[p,]$Genes, ';')))
  }
  gene.all = unique(gene.all)

  if (length(gene.all) !=0) {
    # MAtrix
    M = matrix(0, nrow = length(pathways), ncol = length(gene.all))
    row.names(M) = pathways
    colnames(M) = gene.all

    for (pat in pathways) {
      for (gene in gene.all) {
        if (gene %in% gene.list[[pat]]) {
          M[pat,gene] <- 1
        }
      }
    }

    if (length(pathways) >1) {
      # Jaccard dist
      Jacard.Matrix <- distance(M, method = "jaccard")
      if (length(pathways)==2) {
        Jacard.Matrix_new = as.matrix(rbind(c(0,Jacard.Matrix),c(Jacard.Matrix,0)))
        Jacard.Matrix = Jacard.Matrix_new
      }

      row.names(Jacard.Matrix) <- pathways
      colnames(Jacard.Matrix) <- pathways
      w=30; h=80; fs = 5; #cutree_rows_N = 10
      myb = seq(0,1,by = 0.01)
      myc = colorRampPalette(brewer.pal(n = 7, name ="RdYlBu"))(length(myb))
      pheatmap(Jacard.Matrix,
               border_color = 'darkgrey',
               color = myc,
               breaks = myb,
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               cutree_rows = cutree_rows_N,
               show_colnames = FALSE,
               main = paste(file,'- Jaccard distance heatmap'),
               fontsize = 8,
               fontsize_row = fs,
               filename = paste(dir,file,'cherry_JaccardDist.pdf', sep=''),
               width=w, height=h)
    }
  }
}




