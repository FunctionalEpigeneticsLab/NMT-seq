library(SeuratDisk)
library(Seurat)

convert_srat_h5 <- function(rds, h5Seuratfile){
  tpm_ana <- readRDS(rds)
  # Access raw counts from the "RNA" assay
  raw_counts <- GetAssayData(object = tpm_ana, assay = "RNA", layer = "counts")
  raw_counts[1:5,1:5]
  # Access metadata
  meta <- tpm_ana@meta.data
  # Create new object
  options(Seurat.object.assay.version = "v3")
  srat <- CreateSeuratObject(counts = raw_counts, project = "NMT")
  # Add back metadata
  srat@meta.data <- meta
  print(all(Cells(srat)==rownames(srat@meta.data)))
  print(srat)
  # Save file
  SaveH5Seurat(srat, filename = h5Seuratfile, overwrite = TRUE)
  # Convert
  Convert(h5Seuratfile, dest = "h5ad", overwrite = TRUE)
}

setwd('')

convert_srat_h5(rds="tpm.rds",h5Seuratfile="tpm.h5Seurat")

convert_srat_h5(rds="tpm_ana.rds",h5Seuratfile="tpm_ana.h5Seurat")
