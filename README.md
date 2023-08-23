# Peripheral nerve injury
scRNAseq and BulkRNAseq scripts for "*Structured wound angiogenesis instructs mesenchymal barrier compartments in the regenerating nerve*" publication 

### Raw counts data in GEO. 

## Directory structure


 * [scRNAseq](./scRNAseq)
   * [Tomato_scRNAseq_full_ds.R](./scRNAseq/Tomato_scRNAseq_full_ds.R) [TdTomato scRNAseq full dataset (containing also some fibroblasts, pericites, and LEC cells) prepropressing.] 
   * [Tomato_scRNAseq_EC_only_ds.R](./scRNAseq/Tomato_scRNAseq_EC_only_ds.R) [Analysis and carachterization of EC cell only clusters, isolated from the TdTomato scRNAseq full dataset]
   * [Carr_reanalysis.R](./scRNAseq/Carr_reanalysis.R) [Re-analysis of Carr et al. 2019 scRNAseq dataset. GSE120678]
   * [Toma_reanalysis.R](./scRNAseq/Toma_reanalysis.R) [Re-analysis of Toma et al. 2020 dataset. GSE147285]
   * [Kalinski_reanalysis.R](./scRNASeq/Kalinski_reanalysis.R) [Re-analysis of Kalinski et al. 2020 dataset. GSE153762]
   * [Toma_Carr_Kalinski_Bonanomi_alluvial.R](./scRNAseq/Toma_Carr_Kalinski_Bonanomi_alluvial.R) [Alluvial plot of the relative abundance of EC subtypes in the intact nerve or at different days post injury.]
   * [Milich_EC_label_transfer_alluvial.R](./scRNAseq/Milich_EC_label_transfer_alluvial.R ) [Re-analysis of Milich et al. GSE162610 spinal cord injury dataset]
   * [Ydens_reanalysis.R](./scRNAseq/Ydens_reanalysis.R) [Re-analyisis of Ydens et al. 2020]
   * [Merge_EC_MACs_SCHW_MES_Carr.R](./scRNAseq/Merge_EC_MACs_SCHW_MES_Carr.R) [Cells from different tissues and datasets were merged and clustered]
   * [CellChat_EC_MACs_SCHW_MES_Carr.R](./scRNAseq/CellChat_EC_MACs_SCHW_MES_Carr.R) [Cell-Cell communication analysis with CellChat. EC, FIBRO, MACS, SCHWANN]
   * [monocle_pseudotime.R](./scRNAseq/monocle_pseudotime.R) [EC pseudotime analysis with monocle2]
   * [dendrogram.R](./scRNAseq/dendrogram.R) [EC dendrogram representation]

 * [bulkRNAseq](./bulkRNAseq) 
   * [bulkRNASeq_Cherry.R](./bulkRNAseq/bulkRNASeq_Cherry.R) [DGE analysis of our 7DPI WT and KO Cherry EC samples]
   * [bulkRNASeq_Tomato.R](./bulkRNAseq/bulkRNASeq_Tomato.R) [DGE analysis of our intact, 7 DPI, 14 DPI TdTomato EC samples]
   * [metadata_cherry_tomato.xlsx](./bulkRNAseq/metadata_cherry_tomato.xlsx) [sample metadatda for our Cherry and TdTomato EC samples ]
   * [bulkRNASeq_Tomato_Corada.R](./bulkRNAseq/bulkRNASeq_Tomato_Corada.R) [DGE analysis of our intact, 7 DPI, 14 DPI TdTomato samples + Corada et al. 2019 embryo and Adult samples ]
   * [metadata_tomato_corada.xlsx](./bulkRNAseq/metadata_tomato_corada.xlsx) [sample metadatda for our TdTomato samples and Corada et al. 2019 samples]
   * [EWCE.R](./bulkRNAseq/EWCE.R) [Expression Weighted Cell Type Enrichment analysis to check the distribution of DEGs among EC subtypes]
   * [Associate_celltype_with_gene_name.R](./bulkRNAseq/Associate_celltype_with_gene_name.R) [DEG gene annotation according to scRNAseq expression]





