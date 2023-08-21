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

getwd()
```{r pressure, echo=FALSE}
plot(pressure)
```


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


# path to the count files containing Cherry and Tomato data
filecount_1 = '/Users/maurizio.aurora/Documents/filecount_1/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged/all.counts.gz'
filecount_2 ='/Users/maurizio.aurora/Documents/filecount_2/samples/all/runs/all/fastq/merge-by-read/trimmed/trimmomatic/mapped/STAR/merged/featureCounts/merged/all.counts.gz'

# Min number of replica in each group
Nreplica= 4

# import metadata for Cherry and Tomato samples
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

# save tables as a multisheet excel (cherry and tomato)
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

#### check control genes to be sure that sample's sex matches what expected
#### first good practice check to spot sample switch (at macro level)

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



#####################################
####### only Tomato samples #########
### intact s.nerve, D7 PI, D14 PI ###
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
#top 500 most variable genes PCA
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


pdf('PCA_top500rpkm_Tomato.pdf', width = 14, height = 10)
pca
dev.off()



annotation_column <- metadata_d[,2:(dim(metadata_d)[2])]

row.names(annotation_column) <- metadata_d[,1]
options(repr.plot.width=12, repr.plot.height=10)


crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)

# Fig2 A Heatmap of top 500 most variable genes
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



########################

# MiniHeatmaps of Fig2 E

angiogenesis = readLines("/Users/maurizio.aurora/Downloads/angiogenesis.txt")
barrier = readLines("/Users/maurizio.aurora/Downloads/barrier_go.txt")
endmt = readLines("/Users/maurizio.aurora/Downloads/endmt.txt")
inflammation = readLines("/Users/maurizio.aurora/Downloads/inflammation.txt")
glycolisis = readLines("/Users/maurizio.aurora/Downloads/glycolisis.txt")
lipidmetabolism = readLines("/Users/maurizio.aurora/Downloads/lipidmetabolism.txt")
transporters = readLines("/Users/maurizio.aurora/Downloads/transporters.txt")


logrpkm = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/COUNTS_BonanomiD_1292_RNASeq.xlsx",
  sheet = 4,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)


metadata = read.xlsx(
  "/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/metadata.xlsx",
  sheet = 1,
  startRow = 1,
  colNames = TRUE,
  rowNames = TRUE,
)


metadata_d = metadata[metadata$reporter == 'Tomato',]
metadata$nerve = as.factor(metadata$nerve)

logrpkm = logrpkm[,row.names(metadata_d)]

comparison = list(crush_d14 = c("crush_d14"),
                  intact = c("intact"),
                  crush_d7 = c("crush_d7"))


metadata_crush_d14 = metadata_d[metadata_d$nerve %in% comparison[[1]],]
crush_d14 = logrpkm[,row.names(metadata_crush_d14)]
crush_d14$average = rowMeans(crush_d14)

metadata_crush_intact = metadata_d[metadata_d$nerve %in% comparison[[2]],]
intact = logrpkm[,row.names(metadata_crush_intact)]
intact$average = rowMeans(intact)

metadata_crush_d7 = metadata_d[metadata_d$nerve %in% comparison[[3]],]
crush_d7 = logrpkm[,row.names(metadata_crush_d7)]
crush_d7$average = rowMeans(crush_d7)


#############

# barrier D14

crush_d14_subset_barrier = crush_d14[row.names(crush_d14) %in% barrier,]
crush_d14_subset_barrier$annotation = "barrier"
crush_d14_subset_barrier$sample = "crush_d14"
crush_d14_subset_barrier =  crush_d14_subset_barrier[c("average","annotation","sample")]
crush_d14_subset_barrier$gene = row.names(crush_d14_subset_barrier)

#############

# barrier intact

intact_subset_barrier = intact[row.names(intact) %in% barrier,]
intact_subset_barrier$annotation = "barrier"
intact_subset_barrier$sample = "intact"

intact_subset_barrier =  intact_subset_barrier[c("average","annotation","sample")]
intact_subset_barrier$gene = row.names(intact_subset_barrier)

##############

# barrier D7

crush_d7_subset_barrier = crush_d7[row.names(crush_d7) %in% barrier,]
crush_d7_subset_barrier$annotation = "barrier"
crush_d7_subset_barrier$sample = "crush_d7"
crush_d7_subset_barrier =  crush_d7_subset_barrier[c("average","annotation","sample")]
crush_d7_subset_barrier$gene = row.names(crush_d7_subset_barrier)

###############
###############

# angiogenesis D14

crush_d14_subset_angiogenesis = crush_d14[row.names(crush_d14) %in% angiogenesis,]
crush_d14_subset_angiogenesis$annotation = "angiogenesis"
crush_d14_subset_angiogenesis$sample = "crush_d14"
crush_d14_subset_angiogenesis =  crush_d14_subset_angiogenesis[c("average","annotation","sample")]
crush_d14_subset_angiogenesis$gene = row.names(crush_d14_subset_angiogenesis)

############

# angiogenesis intact

intact_subset_angiogenesis = intact[row.names(intact) %in% angiogenesis,]
intact_subset_angiogenesis$annotation = "angiogenesis"
intact_subset_angiogenesis$sample = "intact"
intact_subset_angiogenesis =  intact_subset_angiogenesis[c("average","annotation","sample")]
intact_subset_angiogenesis$gene = row.names(intact_subset_angiogenesis)

##############

# angiogenesis D7

crush_d7_subset_angiogenesis = crush_d7[row.names(crush_d7) %in% angiogenesis,]
crush_d7_subset_angiogenesis$annotation = "angiogenesis"
crush_d7_subset_angiogenesis$sample = "crush_d7"
crush_d7_subset_angiogenesis =  crush_d7_subset_angiogenesis[c("average","annotation","sample")]
crush_d7_subset_angiogenesis$gene = row.names(crush_d7_subset_angiogenesis)

###############
###############

# endmt D14

crush_d14_subset_endmt = crush_d14[row.names(crush_d14) %in% endmt,]
crush_d14_subset_endmt$annotation = "endmt"
crush_d14_subset_endmt$sample = "crush_d14"
crush_d14_subset_endmt =  crush_d14_subset_endmt[c("average","annotation","sample")]
crush_d14_subset_endmt$gene = row.names(crush_d14_subset_endmt)

############

# endmt intact

intact_subset_endmt = intact[row.names(intact) %in% endmt,]
intact_subset_endmt$annotation = "endmt"
intact_subset_endmt$sample = "intact"
intact_subset_endmt =  intact_subset_endmt[c("average","annotation","sample")]
intact_subset_endmt$gene = row.names(intact_subset_endmt)

##############

# endmt D7

crush_d7_subset_endmt = crush_d7[row.names(crush_d7) %in% endmt,]
crush_d7_subset_endmt$annotation = "endmt"
crush_d7_subset_endmt$sample = "crush_d7"
crush_d7_subset_endmt =  crush_d7_subset_endmt[c("average","annotation","sample")]
crush_d7_subset_endmt$gene = row.names(crush_d7_subset_endmt)

###############
###############

# inflammation D14

crush_d14_subset_inflammation = crush_d14[row.names(crush_d14) %in% inflammation,]
crush_d14_subset_inflammation$annotation = "inflammation"
crush_d14_subset_inflammation$sample = "crush_d14"
crush_d14_subset_inflammation =  crush_d14_subset_inflammation[c("average","annotation","sample")]
crush_d14_subset_inflammation$gene = row.names(crush_d14_subset_inflammation)

############

# inflammation intact

intact_subset_inflammation = intact[row.names(intact) %in% inflammation,]
intact_subset_inflammation$annotation = "inflammation"
intact_subset_inflammation$sample = "intact"
intact_subset_inflammation =  intact_subset_inflammation[c("average","annotation","sample")]
intact_subset_inflammation$gene = row.names(intact_subset_inflammation)

##############

# inflammation D7

crush_d7_subset_inflammation = crush_d7[row.names(crush_d7) %in% inflammation,]
crush_d7_subset_inflammation$annotation = "inflammation"
crush_d7_subset_inflammation$sample = "crush_d7"
crush_d7_subset_inflammation =  crush_d7_subset_inflammation[c("average","annotation","sample")]
crush_d7_subset_inflammation$gene = row.names(crush_d7_subset_inflammation)


#############
#############

# glycolisis D14

crush_d14_subset_glycolisis = crush_d14[row.names(crush_d14) %in% glycolisis,]
crush_d14_subset_glycolisis$annotation = "glycolisis"
crush_d14_subset_glycolisis$sample = "crush_d14"
crush_d14_subset_glycolisis =  crush_d14_subset_glycolisis[c("average","annotation","sample")]
crush_d14_subset_glycolisis$gene = row.names(crush_d14_subset_glycolisis)

############

# glycolisis intact

intact_subset_glycolisis = intact[row.names(intact) %in% glycolisis,]
intact_subset_glycolisis$annotation = "glycolisis"
intact_subset_glycolisis$sample = "intact"
intact_subset_glycolisis =  intact_subset_glycolisis[c("average","annotation","sample")]
intact_subset_glycolisis$gene = row.names(intact_subset_glycolisis)

##############

# glycolisis D7

crush_d7_subset_glycolisis = crush_d7[row.names(crush_d7) %in% glycolisis,]
crush_d7_subset_glycolisis$annotation = "glycolisis"
crush_d7_subset_glycolisis$sample = "crush_d7"
crush_d7_subset_glycolisis =  crush_d7_subset_glycolisis[c("average","annotation","sample")]
crush_d7_subset_glycolisis$gene = row.names(crush_d7_subset_glycolisis)

##############
##############

# lipidmetabolism D14

crush_d14_subset_lipidmetabolism = crush_d14[row.names(crush_d14) %in% lipidmetabolism,]
crush_d14_subset_lipidmetabolism$annotation = "lipidmetabolism"
crush_d14_subset_lipidmetabolism$sample = "crush_d14"
crush_d14_subset_lipidmetabolism =  crush_d14_subset_lipidmetabolism[c("average","annotation","sample")]
crush_d14_subset_lipidmetabolism$gene = row.names(crush_d14_subset_lipidmetabolism)

############

# lipidmetabolism intact

intact_subset_lipidmetabolism = intact[row.names(intact) %in% lipidmetabolism,]
intact_subset_lipidmetabolism$annotation = "lipidmetabolism"
intact_subset_lipidmetabolism$sample = "intact"
intact_subset_lipidmetabolism =  intact_subset_lipidmetabolism[c("average","annotation","sample")]
intact_subset_lipidmetabolism$gene = row.names(intact_subset_lipidmetabolism)

##############

# lipidmetabolism D7

crush_d7_subset_lipidmetabolism = crush_d7[row.names(crush_d7) %in% lipidmetabolism,]
crush_d7_subset_lipidmetabolism$annotation = "lipidmetabolism"
crush_d7_subset_lipidmetabolism$sample = "crush_d7"
crush_d7_subset_lipidmetabolism =  crush_d7_subset_lipidmetabolism[c("average","annotation","sample")]
crush_d7_subset_lipidmetabolism$gene = row.names(crush_d7_subset_lipidmetabolism)

##############
##############

# transporters D14

crush_d14_subset_transporters = crush_d14[row.names(crush_d14) %in% transporters,]
crush_d14_subset_transporters$annotation = "transporters"
crush_d14_subset_transporters$sample = "crush_d14"
crush_d14_subset_transporters =  crush_d14_subset_transporters[c("average","annotation","sample")]
crush_d14_subset_transporters$gene = row.names(crush_d14_subset_transporters)

############

# transporters intact

intact_subset_transporters = intact[row.names(intact) %in% transporters,]
intact_subset_transporters$annotation = "transporters"
intact_subset_transporters$sample = "intact"
intact_subset_transporters =  intact_subset_transporters[c("average","annotation","sample")]
intact_subset_transporters$gene = row.names(intact_subset_transporters)

##############

# transporters D7

crush_d7_subset_transporters = crush_d7[row.names(crush_d7) %in% transporters,]
crush_d7_subset_transporters$annotation = "transporters"
crush_d7_subset_transporters$sample = "crush_d7"
crush_d7_subset_transporters =  crush_d7_subset_transporters[c("average","annotation","sample")]
crush_d7_subset_transporters$gene = row.names(crush_d7_subset_transporters)


#############


ann_colors = list(
  sample = c("Intact" = 'light salmon',
             "D7" = 'tomato',
             "D14" = 'maroon'))



crp <- colorRampPalette(c('dodgerblue4','white','darkred'))
colors = crp(255)
minH = -2; maxH=2
myb = seq(minH, maxH, by = 0.01)
myc <- crp(length(myb))


draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)


crush_d7_subset_transporters$D7 = crush_d7_subset_transporters$average

crush_d14_subset_transporters$D14 = crush_d14_subset_transporters$average

intact_subset_transporters$Intact = intact_subset_transporters$average


test_trans = cbind(intact_subset_transporters, crush_d7_subset_transporters, crush_d14_subset_transporters)
head(test_trans)
test_trans = test_trans[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_trans))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_trans)
head(annotation_column)

#Fig2 E

pdf('MiniHeatmap_Transporters.pdf')
pheatmap(as.matrix(test_trans),
         main = "Transporters",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()

#############

crush_d7_subset_barrier$D7 = crush_d7_subset_barrier$average

crush_d14_subset_barrier$D14 = crush_d14_subset_barrier$average

intact_subset_barrier$Intact = intact_subset_barrier$average


test_barrier = cbind(intact_subset_barrier, crush_d7_subset_barrier, crush_d14_subset_barrier)
head(test_barrier)
test_barr = test_barrier[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_barr))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_barr)
head(annotation_column)


#Fig2 E

pdf('MiniHeatmap_Barrier.pdf')
pheatmap(as.matrix(test_barr),
         main = "Barrier",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()


#############


crush_d7_subset_angiogenesis$D7 = crush_d7_subset_angiogenesis$average

crush_d14_subset_angiogenesis$D14 = crush_d14_subset_angiogenesis$average

intact_subset_angiogenesis$Intact = intact_subset_angiogenesis$average


test_angiogenesis = cbind(intact_subset_angiogenesis, crush_d7_subset_angiogenesis, crush_d14_subset_angiogenesis)
head(test_angiogenesis)
test_ang = test_angiogenesis[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_ang))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_ang)
head(annotation_column)


#Fig2 E

pdf('MiniHeatmap_angiogenesis.pdf')
pheatmap(as.matrix(test_ang),
         main = "Angiogenesis",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()


##########


crush_d7_subset_endmt$D7 = crush_d7_subset_endmt$average

crush_d14_subset_endmt$D14 = crush_d14_subset_endmt$average

intact_subset_endmt$Intact = intact_subset_endmt$average


test_endmt = cbind(intact_subset_endmt, crush_d7_subset_endmt, crush_d14_subset_endmt)
head(test_endmt)
test_endmt = test_endmt[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_endmt))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_endmt)
head(annotation_column)

#Fig2 E

pdf('MiniHeatmap_endmt.pdf')
pheatmap(as.matrix(test_endmt),
         main = "EMT",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()


##########

crush_d7_subset_inflammation$D7 = crush_d7_subset_inflammation$average

crush_d14_subset_inflammation$D14 = crush_d14_subset_inflammation$average

intact_subset_inflammation$Intact = intact_subset_inflammation$average


test_inflammation = cbind(intact_subset_inflammation, crush_d7_subset_inflammation, crush_d14_subset_inflammation)
head(test_inflammation)
test_inflammation = test_inflammation[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_inflammation))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_inflammation)
head(annotation_column)

#Fig2 E

pdf('MiniHeatmap_inflammation.pdf')
pheatmap(as.matrix(test_inflammation),
         main = "Inflammation",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()

############

crush_d7_subset_glycolisis$D7 = crush_d7_subset_glycolisis$average

crush_d14_subset_glycolisis$D14 = crush_d14_subset_glycolisis$average

intact_subset_glycolisis$Intact = intact_subset_glycolisis$average


test_glycolisis = cbind(intact_subset_glycolisis, crush_d7_subset_glycolisis, crush_d14_subset_glycolisis)
head(test_glycolisis)
test_glycolisis = test_glycolisis[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_glycolisis))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_glycolisis)
head(annotation_column)


#Fig2 E

pdf('MiniHeatmap_glycolisis.pdf')
pheatmap(as.matrix(test_glycolisis),
         main = "Glycolisis",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()


#############


crush_d7_subset_lipidmetabolism$D7 = crush_d7_subset_lipidmetabolism$average

crush_d14_subset_lipidmetabolism$D14 = crush_d14_subset_lipidmetabolism$average

intact_subset_lipidmetabolism$Intact = intact_subset_lipidmetabolism$average


test_lipidmetabolism = cbind(intact_subset_lipidmetabolism, crush_d7_subset_lipidmetabolism, crush_d14_subset_lipidmetabolism)
head(test_lipidmetabolism)
test_lipidmetabolism = test_lipidmetabolism[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_lipidmetabolism))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_lipidmetabolism)
head(annotation_column)


#Fig2 E

pdf('MiniHeatmap_lipidmetabolism.pdf')
pheatmap(as.matrix(test_lipidmetabolism),
         main = "Lipid metabolism",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()

########à



crush_d7_subset_lipidmetabolism$D7 = crush_d7_subset_lipidmetabolism$average

crush_d14_subset_lipidmetabolism$D14 = crush_d14_subset_lipidmetabolism$average

intact_subset_lipidmetabolism$Intact = intact_subset_lipidmetabolism$average


test_lipidmetabolism = cbind(intact_subset_lipidmetabolism, crush_d7_subset_lipidmetabolism, crush_d14_subset_lipidmetabolism)
head(test_lipidmetabolism)
test_lipidmetabolism = test_lipidmetabolism[,c("Intact","D7","D14")]

annotation_column = as.data.frame(colnames(test_lipidmetabolism))
colnames(annotation_column)= "sample"
rownames(annotation_column) <- colnames(test_lipidmetabolism)
head(annotation_column)

#Fig2 E

pdf('MiniHeatmap_lipidmetabolism.pdf')
pheatmap(as.matrix(test_lipidmetabolism),
         main = "lipidmetabolism",
         show_rownames = T,
         show_colnames = F,
         cellheight = 50,
         cellwidth = 50,
         cluster_rows = F,
         cluster_cols = F,
         #left_annotation = row_ha,
         scale = 'row',
         fontsize = 15,
         annotation_colors = ann_colors,
         annotation_col = annotation_column)
dev.off()





######### DGE #####################################################

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

