# Peripheral nerve injury
scRNAseq and BulkRNAseq scripts for "Structured wound angiogenesis instructs mesenchymal barrier compartments in the regenerating nerve" publication 

### Raw counts data in GEO. 

## Directory structure


 * [scRNAseq](./scRNAseq)
   * [](./qc_and_somatic_mutation_calling) 
   * [](./Snakefile) [QC, alignment, and somatic calling Snakefile]
   * [](./Plot_MAF.R) [Rscript to summarize, analyze, annotate, and visualize MAF files]
 * [bulkRNAseq](./bulkRNAseq) 
   * [bulkRNASeq_Cherry.R](./bulkRNAseq/bulkRNASeq_Cherry.R) [DGE analysis of our 7DPI WT and KO Cherry EC samples]
   * [bulkRNASeq_Tomato.R](./bulkRNAseq/bulkRNASeq_Tomato.R) [DGE analysis of our intact, 7 DPI, 14 DPI TdTomato EC samples]
   * [metadata_cherry_tomato.xlsx](./bulkRNAseq/metadata_cherry_tomato.xlsx) [sample metadatda for our Cherry and Tomato EC samples ]
   * [bulkRNASeq_Tomato_Corada.R](./bulkRNAseq/bulkRNASeq_Tomato_Corada.R) [DGE analysis of our intact, 7 DPI, 14 DPI samples + Corada et al. 2019 embryo and Adult samples ]
   * [metadata_tomato_corada.xlsx](./bulkRNAseq/metadata_tomato_corada.xlsx) [sample metadatda for our Tomato samples and Corada et al. 2019 samples]
   * [EWCE.R](./bulkRNAseq/EWCE.R)[Expression Weighted Cell Type Enrichment analysis to check the distribution of DEGs among EC subtypes]
   * [Associate_celltype_with_gene_name.R](./bulkRNAseq/Associate_celltype_with_gene_name.R)[DEG gene annotation according to scRNAseq expression]





