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
#install.packages('VennDiagram')
suppressMessages(library('VennDiagram'))
#install.packages("venn")
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
suppressMessages(library(sva))

setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/literature/Corada_embryo")

test=read.table("GSM3473961_EC-E15.5_WT_1_ReadsPerGene.norm.txt", header = T)
test1=read.table("GSM3473962_EC-E15.5_WT_2_ReadsPerGene.norm.txt", header = T)
test2=read.table("GSM3473963_EC-E15.5_WT_3_ReadsPerGene.norm.txt", header = T)
test3=read.table("GSM3473979_EC-Adult_WT_1_ReadsPerGene.norm.txt", header = T)
test4=read.table("GSM3473980_EC-Adult_WT_2_ReadsPerGene.norm.txt", header = T)
test5=read.table("GSM3473981_EC-Adult_WT_3_ReadsPerGene.norm.txt", header = T)

test$E15.5_WT1 = round(test$TPM)
counts = test[,c("gene_name","gene_id","gene_biotype", "chrom", "strand","txstart","txend","length", "E15.5_WT1")]
counts$E15.5_WT2 = round(test1$TPM)
counts$E15.5_WT3 = round(test2$TPM)

counts$Adult_WT1 = round(test3$TPM)
counts$Adult_WT2 = round(test4$TPM)
counts$Adult_WT3 = round(test5$TPM)

filename_xls <- 'counts_corada.xlsx'
write.xlsx(counts,
           file= filename_xls, 
           row.names = F,
           asTable = T)


setwd("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Corada_embryo")


# path to the count file generated with pypette
filecount_1 = '/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/raw_not_corrected_counts_tables_1133_1292/1292/all.counts.gz'
#giulys_filecount='/Users/maurizio.aurora/Documents/BOOB/BonanomiD_1133_RNAseq/dataset/20191211/counts.gz'
filecount_2 ='/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/BonanomiD_1292_RNASeq/7_bioinfo/combined_samples_1133_1292/raw_not_corrected_counts_tables_1133_1292/1133/all.counts.gz'


#head(metadata_file_auri)
# Min number of replica in each group
#Nreplica= 4
Nreplica= 3

# import metadata

metadata_d = read.xlsx("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Corada_embryo/metadata_corada.xlsx")
#remove outliers from metadata

coradas_filecount = read.xlsx("/Users/maurizio.aurora/Dropbox (HSR Global)/WORKSPACE/Bonanomi/Bonanomi_1287_scRNA_injury/7_bioinfo/public_data/Corada_embryo/counts_corada.xlsx")

coradas_filecount$Geneid = coradas_filecount$gene_name

fCounts_corada = coradas_filecount[,c("Geneid", "E15.5_WT1", "E15.5_WT2","E15.5_WT3","Adult_WT1","Adult_WT2","Adult_WT3")]

head(fCounts_corada)
# import counts
annotation <- c('GeneID','Chr','Start','End','Strand','Length')
fCounts_1 <- read.delim(file=filecount_1, header=TRUE, check.names = F)
fCounts_2 <- read.delim(file=filecount_2, header=TRUE, check.names = F)

fCounts_temp <- merge(fCounts_1, fCounts_2, by=c("Geneid","Chr","Start","End","Length","Strand"))

fCounts <- merge(fCounts_temp, fCounts_corada, by=c("Geneid"))
nrow(fCounts) 

colnames(fCounts)

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
colnames(fCounts_bis)


geneidColname <- 'Geneid'
geneidIdx <- which(tolower(annotation) %in% tolower(geneidColname))
rownames(fCountsData_not_adj) <- fCounts[[geneidIdx]]

# Reordering counts matrix to have samples ordered as in metadata
# it assumes that the first column in metadata is sample name
fCountsData_not_adj <- fCountsData_not_adj[,match(metadata[,1], colnames(fCountsData_not_adj))] 

#metadata$Plxnd1 = as.factor(metadata$Plxnd1)
row.names(metadata_d) = metadata_d$sample


cov1=as.factor(metadata_d$nerve)
cov2=as.factor(metadata_d$dataset)
exp_batch=as.factor(metadata_d$exp)
#covar_mat <- cbind(cov1, cov2, cov3)
covar_mat <- cbind(cov2, cov2)

fCountsData_adjusted_p = fCountsData_not_adj[,row.names(metadata_d)]

help(ComBat_seq)
fCountsData_adjusted_adjusted <- ComBat_seq(as.matrix(fCountsData_adjusted_p), batch=exp_batch, group=NULL, full_mod=FALSE) #combined

fCountsData_adjusted_d=fCountsData_adjusted_adjusted

colnames(fCountsData_adjusted_d)

SAVE_variable <- list()
filename_xls <- paste('COUNTS_bonanomi_corada',prj,'.xlsx', sep='')
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


write.xlsx(SAVE_variable,
           file = filename_xls, 
           row.names = T,
           asTable = T, 
           sheetName =variable2save_names)



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
score$nerve = metadata_d$nerve
score$sampleID = metadata_d$sample

pca<- ggplot(score, aes(x=score[,pcx], y=score[,pcy], 
                        color=nerve))+
  geom_label_repel(data= score, aes(x=score[,pcx], y=score[,pcy], 
                                    #color=nerve, label = str_sub(x <- sampleID)),
                                    color=nerve, label = ""),
                   size = 5,  box.padding = unit(0.55, "lines"), point.padding = unit(0.55, "lines"),
                   segment.color = 'grey50') +
  #geom_point(aes(size= IL7))+
  geom_point(size= 12)+
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
  #scale_color_manual(values = brewer.pal(n = 5, name = "Spectral"))  
  scale_color_manual(values=c("#33691E", "#FFAB91", 'purple',"#B388FF",'#4C9900'))


options(repr.plot.width=13, repr.plot.height=7)
pdf('PCA_top500rpkm_Corada_newcolors_squeezed.pdf', width = 8, height = 10)
pca
dev.off()




subset = readLines("/Users/maurizio.aurora/Documents/corada_genes.txt")


filtered = fCountsRPKMTOP[rownames(fCountsRPKMTOP) %in% subset ,]

annotation_column <- metadata_d[,2:(dim(metadata_d)[2])]
row.names(annotation_column) <- metadata_d[,1]
annotation_column_ = annotation_column[, c('dataset', 'nerve')]

ann_colors = list(nerve = c( 'adult' = '#33691E',
                             'crush_d14' = '#FFAB91', 
                             'crush_d7' = 'purple',
                             'embryo' = '#B388FF',
                             'intact' = '#4C9900'),
                  dataset = c('Bonanomi' = '#FFD600','Corada' = '#FF6F00'))

crp <- colorRampPalette(c('dodgerblue4','white','darkred'))

df1 <- annotation_column_ %>% 
  arrange(factor(nerve, levels = c("intact", "crush_d7", "crush_d14", "embryo", "adult")))

re_filtered <- filtered[, rownames(df1)]



HP <- pheatmap::pheatmap(re_filtered,
                         scale = 'row',
                         annotation_col = annotation_column_,
                         cluster_rows = T, 
                         cluster_cols = F, 
                         annotation_colors = ann_colors, 
                         show_rownames=F,
                         show_colnames=F,
                         treeheight_col  = 0,
                         treeheight_row  = 0,
                         display_numbers = F, 
                         border_color=NA,
                         col=colorRampPalette(rev(brewer.pal(n = 11, name =
                                                               "RdYlBu")))(100),
                         filename = 'Heatmap_rpkm_Corada_filtered.pdf',
                         width = 5, height = 5)



#metadata

f = "nerve"
seqc_pvalue = 0.01
dgeResults = list()



comparison = list(intact_vs_embryo = c("intact","embryo"), 
                  crush_d14_vs_embryo = c("crush_d14","embryo"),
                  crush_d7_vs_embryo = c("crush_d7","embryo"),
                  intact_vs_adult = c("intact","adult"), 
                  crush_d14_vs_adult = c("crush_d14","adult"),
                  crush_d7_vs_adult = c("crush_d7","adult"),
                  crush_d14_vs_crush_d7 = c("crush_d14","crush_d7"), 
                  crush_d14_vs_intact = c("crush_d14","intact"),
                  crush_d7_vs_intact = c("crush_d7","intact"),
                  embryo_vs_adult = c("embryo","adult")
                  )

metadata_d[metadata_d$nerve %in% comparison[[1]],]
fCountsData_adjusted_s = fCountsData_adjusted_d[,row.names(metadata_d)]
colnames(fCountsData_adjusted_d)
row.names(metadata_d) = metadata_d$sample


for (i in names(comparison)) {
  
  metadata_s = metadata_d[metadata_d$nerve %in% comparison[[i]],]
  fCountsData_adjusted_s = fCountsData_adjusted_d[,row.names(metadata_s)]
  dds <- DESeqDataSetFromMatrix(
    countData = fCountsData_adjusted_s,
    colData  = metadata_s,
    design   = as.formula('~nerve'))
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
  main.factor = "nerve"
  
  
  
  pcaData <- plotPCA(vsd, intgroup=c(main.factor),returnData=TRUE)
  pcaData$exp=metadata_s$exp
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  #ggplot(pcaData, aes(PC1, PC2, color=Plxnd1, shape=exp))
  PCA = ggplot(pcaData, aes(PC1, PC2, color=nerve, shape=exp)) +
    geom_label_repel(data= pcaData, aes(PC1, PC2, color=nerve, label = str_sub(name,-5,-1)),
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


f = 'DGE_results_Corada'
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
           file = 'DGE_results_Corada.xlsx', 
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
  pdf(paste(f,'/','MAplot_',nerve,'.pdf',sep=''),width=8, height=6.5)
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
    geom_label_repel(data= results[down_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #geom_label_repel(data= results["nerve",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     label = row.names(results[down_label,]), size = 3,
                     #label = row.names(results[up_label,]), size = 3,
                     #label = row.names(results["nerve",]), size = 3,
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
  
  #pdf(paste('VulcanoPlot_Up',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  pdf(paste('VulcanoPlot_Dw_',nerve,'.pdf',sep=''),width=8, height=6.5)
  plot(VP)
  dev.off()  
  
  
  
  
  # Vulcano plot 
  VP0 = ggplot(results) +
    geom_point(aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    geom_point(data = subset(results, DE2 == 'mod'),
               aes(x=log2FoldChange, y=-log10(pvalue), color = DE, alpha = DE2), size =3) +
    ggtitle(paste("Vulcano Plot,", nerve)) +
    scale_color_manual(values = mycolors) +
    scale_size_manual(values = mysize) +
    scale_alpha_manual(values = myalpha) +
    geom_hline(yintercept=0, linetype="dashed", color = "darkgrey") +
    geom_vline(xintercept=0, linetype="dashed", color = "darkgrey") +
    geom_label_repel(data= results[up_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #geom_label_repel(data= results[down_label,], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #geom_label_repel(data= results["nerve",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #label = row.names(results[down_label,]), size = 3,
                     label = row.names(results[up_label,]), size = 3,
                     #label = row.names(results["nerve",]), size = 3,
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
  
  print(VP0)
  #pdf(paste(f,'/','VulcanoPlot_',nerve,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Plnd1',nerve,'.pdf',sep=''),width=8, height=6.5)
  pdf(paste('VulcanoPlot_Up_',nerve,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Dw',nerve,'.pdf',sep=''),width=8, height=6.5)
  plot(VP0)
  dev.off()
  
  
  
  
  # Vulcano plot 
  VP1 = ggplot(results) +
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
    geom_label_repel(data= results["nerve",], aes(x = log2FoldChange, y = -log10(pvalue), color = DE),
                     #label = row.names(results[down_label,]), size = 3,
                     #label = row.names(results[up_label,]), size = 3,
                     label = row.names(results["nerve",]), size = 3,
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
  
  print(VP1)
  #pdf(paste(f,'/','VulcanoPlot_',nerve,'.pdf',sep=''),width=8, height=6.5)
  pdf(paste('VulcanoPlot_Plnd1_',nerve,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Up',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  #pdf(paste('VulcanoPlot_Dw',Plxnd1,'.pdf',sep=''),width=8, height=6.5)
  plot(VP1)
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
  
  #A1 <- image_read_pdf(paste(f,'/','MA_VP_',Plxnd1,'.pdf',sep=''), density = 70)
  #image_write(A1, path = paste(f,'/','MA_VP_',Plxnd1,'.tiff',sep=''), format = "tiff")
  
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


