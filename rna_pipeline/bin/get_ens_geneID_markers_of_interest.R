setwd("./metadata")
library(biomaRt)
library(dplyr)

# Alveolar: CLDN18 Endothelial: CLDN5 Fibroblast: COL1A1 B cell: CD79A Myeloid: LYZ T cell: CD3D Cancer: EPCAM
gene_sym <- c('CLDN18','CLDN5','CAPS','COL1A1','CD79A', 'LYZ','CD3D','EPCAM')

# Read in the additional genes for which the DGE is associated with lung cancer
ref_DGE <- read.csv('D:/OneDrive - KU Leuven/phd/NMTseq/metadata/ref_DGE/ref_DGE.tsv', sep = '\t')[,-1]

# Get the genes with log2 fold change > 5
ref_DGE_genes <- ref_DGE[which(ref_DGE$Fold.Change..log2.>5),'Gene']
gene_sym <- c(gene_sym, ref_DGE_genes)

# Get ensembl gene names
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
mapping <- getBM(
  attributes = c('ensembl_gene_id', "external_gene_name"), 
  filters = "external_gene_name",
  values = gene_sym,
  mart = ensembl
)

# Reorder rows in mapping based on gene_sym
mapping <- mapping[match(gene_sym, mapping$external_gene_name), ]
cell_type <- c('Alveolar','Endothelial','Epithelial','Fibroblast',
              'B cell', 'Myeloid','T cell','Cancer')
cell_type <- c(cell_type, rep('Cancer', length(ref_DGE_genes)))
mapping$cell_type = paste(cell_type, ": ", mapping$external_gene_name, sep = "")
mapping
write.table(mapping, file = "markers_of_interest.tsv", row.names=FALSE, sep="\t")
