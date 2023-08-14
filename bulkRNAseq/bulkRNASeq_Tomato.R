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
#setwd(paste("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/",
#            PI,"/",prj,
#            "/7_bioinfo/",
#            sep=''))



getwd()
```{r pressure, echo=FALSE}
plot(pressure)
```


# load libraries
suppressMessages(library("edgeR"))
suppressMessages(library(data.table))
suppressMessages(library("DESeq2"))
suppressMessages(library(openxlsx))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library("RColorBrewer"))
suppressMessages(library(enrichR))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library('VennDiagram'))
suppressMessages(library(venn))
suppressMessages(library(cowplot))
suppressMessages(library(ggpubr))
suppressMessages(library("viridis"))
suppressMessages(library('pals'))
suppressMessages(library(RColorBrewer))
suppressMessages(library(wesanderson))
suppressMessages(library(patchwork))
suppressMessages(library(magick))
suppressMessages(library('philentropy'))
suppressMessages(library('IntClust'))
suppressMessages(library('pheatmap'))
suppressMessages(library(assertr))
suppressMessages(library("remotes"))
suppressMessages(library(GeneOverlap))
suppressMessages(library('stringr'))
suppressMessages(library(ggrepel))



# path to the count file generated with pypette
filecount_1 = '/Users/maurizio.aurora/Documents/filecount_1/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged/all.counts.gz'
#giulys_filecount='/Users/maurizio.aurora/Documents/BOOB/BonanomiD_1133_RNAseq/dataset/20191211/counts.gz'
filecount_2 ='/Users/maurizio.aurora/Documents/filecount_2/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged/all.counts.gz'


# Min number of replica in each group
Nreplica= 4

# import metadata

metadata = read.xlsx("/Users/maurizio.aurora/Documents/BonanomiD_1287_scRNA_injury_TdT/metadata.xlsx")


# import counts
annotation <- c('GeneID','Chr','Start','End','Strand','Length')
fCounts_1 <- read.delim(file=filecount_1, header=TRUE, check.names = F)
fCounts_2 <- read.delim(file=filecount_2, header=TRUE, check.names = F)


fCounts <- merge(fCounts_1, fCounts_2, by=c("Geneid","Chr","Start","End","Length","Strand"))


fCountsData <- fCounts[
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



geneidColname <- 'Geneid'
geneidIdx <- which(tolower(annotation) %in% tolower(geneidColname))
rownames(fCountsData) <- fCounts[[geneidIdx]]

# Reordering counts matrix to have samples ordered as in metadata
# it assumes that the first column in metadata is sample name
fCountsData <- fCountsData[,match(metadata[,1], colnames(fCountsData))]


# biotypes
#BiocManager::install("GEOquery")
#load biotypes.R
source('/Users/maurizio.aurora/Downloads/plot-biotypes.R')

biotypesFile = '/Users/maurizio.aurora/Documents/dataset/181020/gencode.vM22.basic.annotation.BIOTYPES.DICTIONARY.gz'
biotypes_function(
  countsFile    = filecount,
  biotypesFile  = biotypesFile,
  pngFolder     = 'Biotypes/',
  minSamples    = 3,
  filterExp     = TRUE,
  useRpkm       = TRUE,
  plotPie       = TRUE,
  sglSamplePlot = TRUE,
  writeTable    = TRUE,
  perc2plot     = 0.001,
  useGgplot     = TRUE)


SAVE_variable <- list()
filename_xls <- paste('COUNTS_',prj,'.xlsx', sep='')
variable2save_names <- c('all_counts', 'expGenes_counts','expGenes_LogCPM', 'expGenes_LogRPKM', 'expGenes_RPKM')


# all_counts
y <- DGEList(counts=fCountsData, genes = fCountsAnnotation)
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


write.xlsx(SAVE_variable,
           file = filename_xls,
           row.names = T,
           asTable = T,
           sheetName =variable2save_names)



y <- DGEList(counts=fCountsData, genes = fCountsAnnotation)
# rpkm
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$Length)
# filter for expression
keep <- rowSums(cpm(y)>1)>=Nreplica
yf <- y[keep,]


# calculate 500 most variant genes and save their counts in fCountsRPKMTOP
N=500
vary <- apply(fCountsRPKM[keep,],1,var)
vary_s <- sort(vary, decreasing = T)
TOP_N <- names(vary_s[1:N])
yTOP <-  y[TOP_N,]
fCountsRPKMTOP <- fCountsRPKM[TOP_N,]


gene_control = c('Xist','Tsix','Actb','Gapdh','Gusb','B2m','Rps27a','Sry','Ddx3y','Eif2s3y','Eif2s3x')

rpkm.control = fCountsRPKM[which(rownames(fCountsRPKM) %in% gene_control),]

#### heatmap control genes
colors <- colorRampPalette( rev(brewer.pal(11, "RdYlBu")) )(255)
crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)


pheatmap::pheatmap(rpkm.control,
                   cluster_rows = F,
                   cluster_cols = F,
                   main = 'Heatmap of housekeeping genes - RPKM',
                   show_rownames = T,
                   show_columnames = F,
                   fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                   display_numbers = F,
                   col=colors,
                   filename = 'Heatmap_controlGenes.pdf',
                   width = 10, height = 7)




row.names(annotation_column) <- metadata[,1]
options(repr.plot.width=12, repr.plot.height=10)
row.names(annotation_column) <- metadata[,1]

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)

getwd()


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


#plot
pca_plot <- list()

i=0
for (col in 2:dim(metadata)[2]) {
  score$factor <-  metadata[,col]
  i=i+1
  pca_plot[[i]] <- ggplot(score, aes(x=score[,pcx], y=score[,pcy], color=factor))+
    geom_point(size=5)+
    labs(x=xlab, y=ylab, title=paste("PC",pcx," vs PC",pcy," scoreplot",sep="")) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    theme(plot.title = element_text(color="black", size=26, face="bold.italic"),
          axis.text.x = element_text(angle = 0, face = "bold", color = "black", size=22, hjust =.5),
          axis.title.x = element_text(face = "bold", color = "black", size = 24),
          axis.text.y = element_text(angle = 0, face = "bold", color = "black", size=22),
          axis.title.y = element_text(face = "bold", color = "black", size = 24),
          legend.text = element_text(face = "bold", color = "black", size = 16),
          legend.position="right",
          panel.background = element_rect(fill = "white",colour = "black", size = 1, linetype = "solid"))
}


png("pca1.png")
pca_plot[[1]]
dev.off()

png("pca2.png")
pca_plot[[2]]
dev.off()

png("pca3.png")
pca_plot[[3]]
dev.off()

png("pca4.png")
pca_plot[[4]]
dev.off()

metadata$reporter = as.factor(metadata$reporter)
metadata$reporter = factor(metadata$reporter, levels = c('KDR_cherry','Tomato'))

metadata$nerve = as.factor(metadata$nerve)
metadata$nerve = factor(metadata$nerve, levels = c('intact','crush_d7','crush_d14') )

metadata$Plxnd1 = as.factor(metadata$Plxnd1)
metadata$Plxnd1 = factor(metadata$Plxnd1, levels = c('WT','KO','na') )

metadata$exp = as.factor(metadata$exp)
metadata$exp = factor(metadata$exp, levels = c('1292','1133') )


score$exp = metadata$exp
score$Plxnd1 = metadata$Plxnd1
score$nerve = metadata$nerve
score$reporter = metadata$reporter
score$sampleID = metadata$sample



pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy],
                        color=reporter, shape = Plxnd1))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy],
                                    color=reporter, label = str_sub(x <- sampleID)),
                   #color=reporter, label = str_sub(x <- sampleID,7,-1)),
                   size = 6,  box.padding = unit(1, "lines"), point.padding = unit(0.1, "lines"),
                   segment.color = 'grey50') +
  geom_point(size= 8)+
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
  scale_color_manual(values = brewer.pal(n = 3, name = "Dark2"))

options(repr.plot.width=10, repr.plot.height=8.5)

pdf('PCA_top500rpkm_Plxnd1_bis.pdf', width = 10, height = 8.5)
pca
dev.off()



pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy],
                        color=nerve, shape = reporter))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy],
                                    #color=nerve, label = str_sub(x <- sampleID)),
                                    color=nerve, label = sampleID),
                   #color=reporter, label = str_sub(x <- sampleID,7,-1)),
                   size = 6,  box.padding = unit(1, "lines"), point.padding = unit(0.1, "lines"),
                   segment.color = 'grey50') +
  geom_point(size= 8)+
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
  scale_color_manual(values = brewer.pal(n = 3, name = "Dark2"))

options(repr.plot.width=10, repr.plot.height=8.5)

pdf('PCA_top500rpkm_reporter_bis.pdf', width = 10, height = 8.5)
pca
dev.off()



#####################################
####### only Tomato sample ##########
#####################################


metadata_d = metadata[metadata$reporter == 'Tomato',]
metadata$nerve = as.factor(metadata$nerve)

row.names(metadata_d) = metadata_d$sample

fCountsData_d = fCountsData[,row.names(metadata_d)]

y <- DGEList(counts=fCountsData_d, genes = fCountsAnnotation)
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$Length)
keep <- rowSums(cpm(y)>1)>=Nreplica
yf <- y[keep,]
nrow(yf)
N=500
vary <- apply(fCountsRPKM[keep,],1,var)
vary_s <- sort(vary, decreasing = T)
provaA <- names(vary)

TOP_N <- names(vary_s[1:N])
nrow(TOP_N)
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

score$Plxnd1 = metadata_d$Plxnd1
score$nerve = metadata_d$nerve
score$reporter = metadata_d$reporter
score$sampleID = metadata_d$sample

pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy],
                        color=nerve, shape = nerve))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy],
                                    #color=nerve, label = str_sub(x <- sampleID,-5,-1)),
                                    color=nerve, label = sampleID),
                   size = 6,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
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
  scale_color_manual(values = watlington(12))

options(repr.plot.width=13, repr.plot.height=7)
#pca
pdf('PCA_top500rpkm_Tomato.pdf', width = 14, height = 10)
pca
dev.off()

annotation_column <- metadata_d[,2:(dim(metadata_d)[2])]


row.names(annotation_column) <- metadata_d[,1]
options(repr.plot.width=12, repr.plot.height=10)


crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         annotation_colors = ann_colors,
                         cluster_rows = T,
                         cluster_cols = T,
                         show_rownames = F,
                         cutree_cols = 5,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                         display_numbers = F,
                         col=colors,
                         filename = 'Heatmap_500rpkm_Tomato.pdf',
                         width = 14, height = 11)


metadata_d = metadata[metadata$reporter == 'Tomato',]
row.names(metadata_d) = metadata_d$sample
fCountsData_d = fCountsData[,row.names(metadata_d)]


y <- DGEList(counts=fCountsData_d, genes = fCountsAnnotation)
fCountsRPKM = rpkm(y, log=T, gene.length =y$genes$Length)
keep <- rowSums(cpm(y)>1)>=Nreplica
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


score$reporter = metadata_d$reporter
score$nerve = metadata_d$nerve
score$Plxnd1 = metadata_d$Plxnd1
score$sampleID = metadata_d$sample


pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy],
                        color=nerve, shape = nerve))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy],
                                    color=nerve, label = str_sub(x <- sampleID,-5,-1)),
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
pdf('PCA_top500rpkm_Tomato.pdf', width = 14, height = 10)
pca
dev.off()

annotation_column <- metadata_d[,2:(dim(metadata_d)[2])]


score$reporter = metadata_d$reporter
score$nerve = metadata_d$nerve
score$Plxnd1 = metadata_d$Plxnd1
score$sampleID = metadata_d$sample


row.names(annotation_column) <- metadata_d[,1]
options(repr.plot.width=12, repr.plot.height=10)


crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
HP <- pheatmap::pheatmap(fCountsRPKMTOP,
                         scale = 'row',
                         annotation_col = annotation_column,
                         cluster_rows = T,
                         cluster_cols = T,
                         show_rownames = F,
                         fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                         display_numbers = F,
                         col=colors,
                         filename = 'Heatmap_500rpkm_Tomato.pdf',
                         width = 14, height = 11)




f = "nerve"
seqc_pvalue = 0.01
dgeResults = list()


comparison = list(crush_d14_vs_crush_d7 = c("crush_d14","crush_d7"),
                  crush_d14_vs_intact = c("crush_d14","intact"),
                  crush_d7_vs_intact = c("crush_d7","intact"))


metadata_d[metadata_d$nerve %in% comparison[[1]],]
fCountsData_s = fCountsData_d[,row.names(metadata_d)]
colnames(fCountsData_d)
row.names(metadata_d) = metadata_d$sample


for (i in names(comparison)) {
  print(i)

  metadata_s = metadata_d[metadata_d$nerve %in% comparison[[i]],]
  fCountsData_s = fCountsData_d[,row.names(metadata_s)]
  dds <- DESeqDataSetFromMatrix(
    countData = fCountsData_s,
    colData  = metadata_s,
    design   = as.formula('~ exp + nerve'))
  filter <- rowSums(cpm(counts(dds)) >= 1) >= Nreplica
  table(filter)
  
  ddsFiltered <- dds[filter,]
  
  dga <- DESeq(
    object = ddsFiltered,
    test = "Wald",
    fitType = "parametric",
    betaPrior = FALSE,
    minReplicatesForReplace = Inf)

  pdf(paste("Disp_cherry.pdf",i))
  plotDispEsts(dga)
  dev.off()

  alpha = 0.05
  
  print(paste(comparison[[i]][1],"_vs_",comparison[[i]][2],sep=''))
  
  dgeResults.tmp <- results(dga,
                            contrast             = c(f,comparison[[i]][1],comparison[[i]][2]),
                            cooksCutoff          = Inf,
                            independentFiltering = TRUE,
                            alpha                = alpha,
                            pAdjustMethod        = "BH")
  
  summary(dgeResults.tmp)
  
  dgeResults[[i]] <- dgeResults.tmp[order(dgeResults.tmp$pvalue, decreasing = F),]

  #PCA
  vsd <- vst(dga, blind=FALSE)
  main.factor = "nerve"
  pcaData <- plotPCA(vsd, intgroup=c(main.factor),returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  #ggplot(pcaData, aes(PC1, PC2, color=nerve, shape=sex))
  PCA = ggplot(pcaData, aes(PC1, PC2, color=nerve)) +
    geom_label_repel(data= pcaData, aes(PC1, PC2, color=nerve, label = name),
    #geom_label_repel(data= pcaData, aes(PC1, PC2, color=nerve, label = str_sub(name,-5,-1)),
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
    HP <- pheatmap::pheatmap(cpm(counts(dga))[seqcUP,],
                             scale = 'row',
                             annotation_col = annotation_column,
                             annotation_colors = ann_colors,
                             cluster_rows = T,
                             cluster_cols = T,
                             #main = 'Heatmap: 500 most variable genes - RPKM',
                             show_rownames = F,
                             #cutree_rows = 3,
                             cutree_cols = 2,
                             fontsize = 12, fontsize_row = 10, fontsize_col = 14,
                             display_numbers = F,
                             col=colors,
                             filename = paste('Heatmap_seqcUP',i,'.pdf'),
                             width = 10, height = 11 )
  }

}



f = 'DGE_results_Tomato'
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
           file = paste(f,'/DGE_results_Tomato.xlsx', sep=''),
           row.names = F,
           asTable = T,
           startRow = 1,
           sheetName = str_sub(names(dgeResults),1,31))






#plots


n.label = 20
FDR = T
pvalue = 0.01
for (nerve in names(dgeResults)) {
  results = as.data.frame(dgeResults[[nerve]])
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

  NSEQC = min(n.label, length(rownames(results[results$DE == 'SEQCup',])))
  SEQCup_label = rownames(results[results$DE == 'SEQCup',])[1:N]
  NSEQC = min(n.label, length(rownames(results[results$DE == 'SEQdown',])))
  SEQCdown_label = rownames(results[results$DE == 'SEQdown',])[1:N]




  MAplot = ggplot(results) +
    geom_point(aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=baseMean, y=log2FoldChange, color = DE, alpha = DE2), size = 3) +
    xlim(c(0,1.e5)) +
    scale_x_continuous(trans='log10') +
    ggtitle(paste("MAPlot,", nerve)) +
    scale_color_manual(values = mycolors) +
    #scale_size_manual(values = mysize) +
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
  pdf(paste('MAplot_',nerve,'.pdf',sep=''),width=8, height=6.5)
  plot(MAplot)
  dev.off()

  # Vulcano plot
  VP = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", nerve)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    #geom_label_repel(data= results[up_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
    #geom_label_repel(data= results[down_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
    geom_label_repel(data= results["Plxnd1",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results["Plxnd1",]), size = 3,
                     #label = row.names(results[up_label,]), size = 3,
                     #label = row.names(results[down_label,]), size = 3,
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

  print(VP)
  pdf(paste(f,'/','VulcanoPlot_Plxnd1',nerve,'.pdf',sep=''),width=8, height=6.5)
  plot(VP)
  dev.off()

  options(repr.plot.width=14, repr.plot.height=6.5)
  p1 = MAplot ; p2 = VP;
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +
          plot_layout(guides = "collect"))

  pdf(paste(f,'/','MA_VP_',nerve,'.pdf',sep=''),width=16, height=6.5)
  print((p1 + theme(plot.margin = unit(c(0,30,0,0), "pt"))) +
          (p2 + theme(plot.margin = unit(c(0,0,0,30), "pt"))) +
          plot_layout(guides = "collect"))
  dev.off()


}

getwd()

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


#######################
###### ENRICHR ########
#######################


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



#########
require("biomaRt")

# Basic function to convert mouse to human gene names
convertMouseGeneList <- function(x){
  #require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  #print(head(genesV2))
  #print(length(genesV2))
  #humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  #print(head(humanx))
  #return(humanx)
  return(genesV2)
}




###### enrichment #######
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

cutree_rows_values = c(5,12,3,10,5,1)

lf = list.files(paste("./enrichR/", sep=''), pattern=glob2rx("*.xlsx"))
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

      w=30; h=80; fs = 5; cutree_rows_N = 10
      if (file=='b.resting_IL7_vs_b.resting_seqc_up_.xlsx'){w=14; h=10; fs = 4; cutree_rows_N = 3}
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
               fontsize = 12,
               fontsize_row = fs,
               filename = paste(dir,file,'_JaccardDist.pdf', sep=''),
               width=w, height=h)
    }
  }
}






