setwd("D:/OneDrive - KU Leuven/phd/NMTseq/230914_countmatrix/VariantCalling/4.1_Variant")
library(stringr)
library(dplyr)

# filename="GC135739_AAGAGGCA-AAGGAGTA_VariantCalling.tsv"
read.vc <- function(filename){
  df <- read.table(filename, header = F)
  colnames(df) <- c('chr','pos','allele_number','alternate_allele_count','frac',
                    'REF','ALT','depth','DP4','TGT','PL','QUAL')
  print(paste('Number of rows (loci): ',nrow(df)))
  # Split the counts for reference and variant allele
  df[,c('alt1','alt2','ref1','ref2')] <- as.integer(str_split_fixed(df$DP4, ',', 4))
  # Compute the nr of reads for the reference and variant allele
  df$ref_allele_count <- df$ref1 + df$ref2
  df$alt_allele_count <- df$alt1 + df$alt2
  # Generate the column for the identity of each allele
  df$ident <- paste(df$chr, df$pos, df$REF, df$ALT, sep='_')
  rownames(df) <- df$ident
  #df <- select(df, -chr, -pos, -tf, -TGT, -GT)
  # Extract the cell name from the file name
  cellname <- strsplit(filename, split='_V')[[1]][1]
  cellname <- strsplit(cellname, split='./')[[1]][2]
  # Split into two dataframes
  df_ref <- select(df, ident, ref_allele_count)
  colnames(df_ref) <- c('allele',cellname)
  df_alt <- select(df, ident, alt_allele_count)
  colnames(df_alt) <- c('allele',cellname)
  return(list(df_ref, df_alt))
}

getVC.df <- function(allele='ref'){
  index=1
  if(allele=='alt'){
    index=2
  }
  
  # Read the QC result
  qc <- read.csv('D:/OneDrive - KU Leuven/phd/NMTseq/metadata/cell_type_results/qc_202310.csv')
  pass_qc <- qc[which(qc$QC=='Pass'),'X']
  
  # List all the files in the folder
  list_files <- list.files(path = ".", full.names = TRUE, include.dirs=FALSE)
  list_files <- grep("VariantCalling.tsv", list_files, value = TRUE)
  nrfiles <- length(list_files)
  print(paste('Number of VariantCalling.tsv files:',nrfiles))
  
  # Extract the cell name from the file name
  list_cellname <- sapply(strsplit(list_files, split='_V'), function(x) x[[1]])
  list_cellname <- sapply(strsplit(list_cellname, split='./'), function(x) x[[2]])
  
  # Keep only the files for cells that passed QC
  idx_pass_qc <- which(list_cellname %in% pass_qc)
  list_files <- list_files[idx_pass_qc]
  print(paste('Number of cells that passed QC in these cells:',length(idx_pass_qc)))
  
  # Initialize
  df <- read.vc(list_files[1])[[index]]
  # Loop through the feature count files and read them into data frames
  # Perform outer joins on the rownames (first columns)
  for (file in list_files[-1]) {
    df_next <- read.vc(file)[[index]]
    df <- merge(x = df, y = df_next, by = "allele", all = TRUE)
  }
  without_missing = length(which(complete.cases(df)))
  print(paste('Number of alleles without missing values:',without_missing))
  return(df)
}

merged_ref <- getVC.df('ref')
merged_alt <- getVC.df('alt')

# Check the alleles without missing value
head(merged_ref[which(complete.cases(merged_ref)),])
# Check the alleles without missing value
head(merged_alt[which(complete.cases(merged_alt)),])

# Save the R image
save.image("20231003_vcf_merge.RData")
