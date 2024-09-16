setwd("D:/OneDrive - KU Leuven/phd/NMTseq/variant/VariantCalling")
load("vcf_20231017.RData")

library(stringr)
library(dplyr)
library(jsonlite)
library(tidyr)

TP53=c('chr17',7661779,7687538)
KRAS=c('chr12',25205246,25250936)
EGFR=c('chr7',55019017,55211628)
ERBB2=c('chr17',39687914,39730426)
BRAF=c('chr7',140719327,140924929)

read.vc.df <- function(filename){
  df <- read.table(filename, header = F)
  colnames(df) <- c('cell', 'chr','pos','allele_number','alternate_allele_count','frac',
                    'ALT','REF','depth','DP4','TGT','PL','QUAL')
  df[,c('a1','a2')] <- str_split_fixed(df$TGT, '/', 2)
  return(df)
}

filter.vc <- function(df, gene){
  df_gene <- df[which(df$chr==gene[1]),]
  df_gene <- df_gene[which(df_gene$pos>=as.integer(gene[2])),]
  df_gene <- df_gene[which(df_gene$pos<=as.integer(gene[3])),]
  return(df_gene)
}

check.all.vc <- function(){
  # List all the files in the folder
  list_files <- list.files(path = ".", full.names = TRUE, include.dirs=FALSE)
  list_files <- grep("VariantCalling.tsv", list_files, value = TRUE)
  nrfiles <- length(list_files)
  print(paste('Number of VariantCalling.tsv files:',nrfiles))

  ## Keep only the files for cells that passed QC
  # Extract the cell name from the file name
  #list_cellname <- sapply(strsplit(list_files, split='_V'), function(x) x[[1]])
  #list_cellname <- sapply(strsplit(list_cellname, split='./'), function(x) x[[2]])
  # Read the QC result
  #qc <- read.csv('D:/OneDrive - KU Leuven/phd/NMTseq/metadata/cell_type_results/qc_202310.csv')
  #pass_qc <- qc[which(qc$QC=='Pass'),'X']
  #idx_pass_qc <- which(list_cellname %in% pass_qc)
  #rm(list_cellname)
  #list_files <- list_files[idx_pass_qc]
  #print(paste('Number of cells that passed QC in these cells:',length(idx_pass_qc)))
  
  tp53_list = list()
  kras_list = list()
  egfr_list = list()
  erbb2_list = list()
  braf_list = list()
  
  i=1
  for (file in 1:length(list_files)) {
    file = list_files[i]
    GC_code = strsplit(file, split='_V')[[1]][1]
    GC_code = strsplit(GC_code, split='./')[[1]][2]
    print(GC_code)
    df <- read.vc.df(file)
    
    # Keep only the files for cells that passed QC
    #df <- df[which(df$cell %in% pass_qc),]
    
    # Get the dfs corresponding to the 3 genes
    tp53 <- filter.vc(df, TP53)
    kras <- filter.vc(df, KRAS)
    egfr <- filter.vc(df, EGFR)
    erbb2 <- filter.vc(df, ERBB2)
    braf <- filter.vc(df, BRAF)
    
    if (nrow(tp53)>0){
      tp53[,'gene'] <- 'TP53'
      tp53_list[[GC_code]] <- as.data.frame(tp53)
    }
    if (nrow(kras)>0){
      kras[,'gene'] <- 'KRAS'
      kras_list[[GC_code]] <- as.data.frame(kras)
    }
    if (nrow(egfr)>0){
      egfr[,'gene'] <- 'EGFR'
      egfr_list[[GC_code]] <- as.data.frame(egfr)
    }
    if (nrow(erbb2)>0){
      erbb2[,'gene'] <- 'ERBB2'
      erbb2_list[[GC_code]] <- as.data.frame(erbb2)
    }
    if (nrow(braf)>0){
      braf[,'gene'] <- 'BRAF'
      braf_list[[GC_code]] <- as.data.frame(braf)
    }
    i=i+1
  }

  df_tp53 <- bind_rows(tp53_list)
  df_kras <- bind_rows(kras_list)
  df_egfr <- bind_rows(egfr_list)
  df_erbb2 <- bind_rows(erbb2_list)
  df_braf <- bind_rows(braf_list)
  df_all <- rbind(df_tp53,df_kras,df_egfr,df_erbb2,df_braf)
  return(df_all)
}

all.vc <- check.all.vc()
all.vc$mutation <- paste(all.vc$gene,all.vc$chr,all.vc$pos,all.vc$TGT,all.vc$PL,all.vc$DP4, sep=':')
all.vc <- all.vc[,-c(4,5,6)]
var_sites <- paste(all.vc$chr, all.vc$pos, sep = '_')
all.vc$mutsites <- var_sites
length(unique(var_sites))
head(all.vc)
unique(all.vc$gene)
grep('397247',all.vc[which(all.vc$gene=='ERBB2'),'pos'],value=TRUE)
grep('140753355',all.vc[which(all.vc$gene=='BRAF'),'pos'],value=TRUE)
grep('25245350',all.vc[which(all.vc$gene=='KRAS'),'pos'],value=TRUE) #
grep('2522734',all.vc[which(all.vc$gene=='KRAS'),'pos'],value=TRUE)

all.vc[which(all.vc$pos=='25245350'),]

### Compute the total number of variants depth (number of read-covered bases that show the variant alleles)
all.vc$coverage <- rowSums(matrix(as.numeric(unlist(strsplit(str_split_fixed(all.vc$mutation, ":", 6)[,6], ","))), ncol = 4, byrow = TRUE))
all.vc$coverage_var <- rowSums(matrix(as.numeric(unlist(strsplit(str_split_fixed(all.vc$mutation, ":", 6)[,6], ","))), ncol = 4, byrow = TRUE)[,c(3,4)])
all.vc$mut <- paste(all.vc$mutsites, all.vc$ALT, sep='_')

# Select unique rows
all.vc <- all.vc %>% distinct()
dim(all.vc)


#=====================================================================
#### Reference SNVs from databases
combine_files <- function(dir){
  list_files <- list.files(path = dir, full.names = TRUE, include.dirs=FALSE)
  nrfiles <- length(list_files)
  print(paste('Number of files:',nrfiles))
  
  df_list = list()
  
  i=1
  for (file in 1:length(list_files)) {
    file = list_files[i]
    df <- read.csv(file, sep='\t')
    df_list[[i]] <- as.data.frame(df)
    i=i+1
  }
  df_full <- bind_rows(df_list)
  return(df_full)
}

### COSMIC, the Catalogue Of Somatic Mutations In Cancer
# Nonsense substitution	681 (13.56%) p.Q354*
# Missense substitution	3963 (78.90%) p.D391N
# Inframe insertion	10 (0.20%) p.A307_L308insASFLS, p.D281_R282insP
# Inframe deletion	76 (1.51%) p.N345_G356del, 	p.E287del, p.G465_N466del
# Frameshift insertion	107 (2.13%) p.D393Rfs*78, p.E336*, p.F341Efs*7
# Frameshift deletion	298 (5.93%) p.K382Nfs*40, p.E343Gfs*2
cosmic0 <- combine_files("D:/OneDrive - KU Leuven/phd/NMTseq/metadata/ref_variant/COSMIC")
dim(cosmic0)
colnames(cosmic0)

# Create a new column based on patterns
cosmic1 <- cosmic0 %>%
  mutate(Mutation_Type = case_when(
    grepl("^p\\.[A-Z]\\d+\\*$", AA.Mutation) ~ "Nonsense substitution",
    grepl("^p\\.[A-Z]\\d+[A-Z]$", AA.Mutation) ~ "Missense substitution",
    grepl("^p\\.[A-Z]\\d+_[A-Z]\\d+ins[A-Z]+$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.[A-Z]\\d+_[A-Z]\\d+ins\\?$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.[A-Z]\\d+_[A-Z]\\d+dup$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.[A-Z]\\d+dup$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.\\d+_\\d+ins\\?$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.\\(\\d+_\\d+\\)ins\\?$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.\\([A-Z]\\d+\\)ins\\?$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.\\(\\d+_\\d+\\)ins[A-Z]+$", AA.Mutation) ~ "Inframe insertion",
    grepl("^p\\.[A-Z]\\d+_[A-Z]\\d+del$", AA.Mutation) ~ "Inframe deletion",
    grepl("^p\\.[A-Z]\\d+del$", AA.Mutation) ~ "Inframe deletion",
    grepl("^p\\.\\(\\d+_\\d+\\)del\\?$", AA.Mutation) ~ "Inframe deletion",
    grepl("^p\\.[A-Z]\\d+_[A-Z]\\d+delins[A-Z*]+$", AA.Mutation) ~ "Inframe deletion",
    grepl("^p\\.[A-Z]\\d+[A-Z*]fs\\*\\d+\\b", AA.Mutation) ~ "Frameshift", 
    grepl("*fs*", AA.Mutation) ~ "Frameshift", 
    TRUE ~ "Other"
  ))
cosmic1[intersect(which(cosmic1$Mutation_Type=="Frameshift"), which(grepl("ins|dup$", cosmic1$CDS.Mutation))),'Mutation_Type'] = "Frameshift insertion"
cosmic1[intersect(which(cosmic1$Mutation_Type=="Frameshift"), which(grepl("del$", cosmic1$CDS.Mutation))),'Mutation_Type'] = "Frameshift deletion"
cosmic <- cosmic1[which(cosmic1$Histology != "Other"),c(1,19,6,7,20,11,12,16)]
coor <- str_split_fixed(cosmic$Genomic.Co.ordinates, '\\..', 2)[,1]
coor2 <- str_split_fixed(cosmic$Genomic.Co.ordinates, '\\..', 2)[,2]
cosmic[,c('chr','site')] <- str_split_fixed(coor, ':', 2)
cosmic$chr <- paste('chr', cosmic$chr, sep = '')
cosmic$mutsites <- paste(cosmic$chr, cosmic$site, sep='_')
#cosmic$mutsites <- paste(cosmic$chr, coor2, sep='_')
cosmic <- cosmic[,-c(9,10)]
rm(cosmic0, cosmic1, coor)
tail(cosmic)

# Take a look at the final dataframe
table(cosmic$Mutation_Type)
cosmic[cosmic$Mutation_Type=="Other",c("AA.Mutation","CDS.Mutation")]
unique(cosmic$Histology)
unique(cosmic$Histology.Subtype.1)
head(cosmic)

intersect(var_sites,cosmic$mutsites)
cosmic_df0 <- cosmic[which(cosmic$mutsites %in% var_sites),]
cosmic_df1 <- data.frame(mutsites = cosmic_df0$mutsites, snv = paste("COSMIC",cosmic_df0$Genomic.Co.ordinates,cosmic_df0$CDS.Mutation,cosmic_df0$Mutation_Type,cosmic_df0$Histology.Subtype.1,cosmic_df0$Somatic.Status,sep = ':'))
cosmic_df1 <- cosmic_df1%>%
  group_by(mutsites) %>%
  mutate(occurrence = row_number()) %>%
  ungroup()
cosmic_df1 <- as.data.frame(cosmic_df1)
cosmic_df1 <- reshape(cosmic_df1, idvar="mutsites", timevar = "occurrence", direction="wide")
cosmic_df <- data.frame(mutsites = cosmic_df1$mutsites ,snv = str_split_fixed(paste(cosmic_df1$snv.1,cosmic_df1$snv.2,cosmic_df1$snv.3,cosmic_df1$snv.4,sep=";"), ';NA',2)[,1])
cosmic_df

# Merge with the query data
m.cosmic <- merge(all.vc[,c('cell','mutation',"mutsites")], cosmic_df, by="mutsites")

### Read in the variants downloaded from TCGA
read_TCGA <- function(jsonfile){
  TCGA_gene <- as.data.frame(fromJSON(jsonfile))[,c(4,2,1)]
  chrs <- str_split_fixed(TCGA_gene$genomic_dna_change, ":",2)[,1]
  part2 <- str_split_fixed(TCGA_gene$genomic_dna_change, ":",2)[,2]
  poss <- str_split_fixed(part2, "g.", 2)[,2]
  positions <- sapply(str_extract_all(poss, "\\d+"), function(x) ifelse(length(x) > 0, x[1], NA))
  mutsites.gene <- paste(chrs, positions, sep = '_')
  TCGA_gene$mutsites <- mutsites.gene
  return(TCGA_gene)
}

TCGA <- read_TCGA("D:/OneDrive - KU Leuven/phd/NMTseq/metadata/ref_variant/TCGA/TCGA_5genes.json")
all.vc[which(var_sites %in% TCGA$mutsites),]
TCGA[which(TCGA$mutsites %in% var_sites),]

list_trans <- lapply(TCGA[,'consequence'], function(x){return(x$transcript)})

get_most_frequent_element <- function(df) {
  # df = list_trans[[2]]
  anno = df$annotation
  polyphen_impact = paste('polyphen_impact',tail(names(sort(table(anno$polyphen_impact))), 1),sep='=')
  sift_impact = paste('sift_impact',tail(names(sort(table(anno$sift_impact))), 1),sep='=')
  vep_impact = paste('vep_impact',tail(names(sort(table(anno$vep_impact))), 1),sep='=')
  consequence_type = paste('consequence_type',tail(names(sort(table(df$consequence_type))), 1),sep='=')
  string = paste(polyphen_impact,sift_impact,vep_impact,consequence_type,sep=';')
  return(string)
}
TCGA$consequence = unlist(lapply(list_trans, get_most_frequent_element))

TCGA_df0 <- TCGA[which(TCGA$mutsites %in% var_sites),-2] %>%
  mutate(consequence = sapply(consequence, paste, collapse = ", "))
TCGA_df <- data.frame(mutsites = TCGA_df0$mutsites, snv = paste('TCGA',TCGA_df0$genomic_dna_change, TCGA_df0$consequence,sep=':'))
rm(TCGA_df0)
head(TCGA_df)
m.tcga <- merge(all.vc[,c('cell','mutation',"mutsites")],TCGA_df, by="mutsites")
head(TCGA)
### Germline mutations from clinivar
combine_files <- function(dir){
  read_clinivar <- function(file){
    gene_cl <- read.csv(file, sep = '\t')
    mutsites.gene <- paste(gene_cl$GRCh38Chromosome, gene_cl$GRCh38Location, sep = '_')
    mutsites.gene <- str_split_fixed(mutsites.gene, ' - ', 2)[,1]
    mutsites.gene <- paste('chr', mutsites.gene, sep = '')
    gene_cl$mutsites <- mutsites.gene
    # remove entries with 'benign' from 'Clinical significance'
    gene_cl <- gene_cl[-grep("enign", gene_cl$Clinical.significance..Last.reviewed.),]
    return(gene_cl)
  }
  
  list_files <- list.files(path = dir, full.names = TRUE, include.dirs=FALSE)
  nrfiles <- length(list_files)
  print(paste('Number of files:',nrfiles))
  df_list = list()
  i=1
  for (file in 1:length(list_files)) {
    file = list_files[i]
    df <- read_clinivar(file)
    df_list[[i]] <- as.data.frame(df)
    i=i+1
  }
  df_full <- bind_rows(df_list)
  return(df_full)
}

cl <- combine_files('D:/OneDrive - KU Leuven/phd/NMTseq/metadata/ref_variant/clinvar')
cl_df0 <- cl[which(cl$mutsites %in% var_sites),c('mutsites','Name','dbSNP.ID','Condition.s.','Clinical.significance..Last.reviewed.')]
cl_df1 <- data.frame(mutsites = cl_df0$mutsites, snv = paste("Clinivar",cl_df0$Name, cl_df0$dbSNP.ID, cl_df0$Condition.s.,cl_df0$Clinical.significance..Last.reviewed.,sep = ':'))
cl_df1 <- cl_df1%>%
  group_by(mutsites) %>%
  mutate(occurrence = row_number()) %>%
  ungroup()
cl_df1 <- as.data.frame(cl_df1)
cl_df1 <- reshape(cl_df1, idvar="mutsites", timevar = "occurrence", direction="wide")
cl_df <- data.frame(mutsites = cl_df1$mutsites ,snv = str_split_fixed(paste(cl_df1$snv.1,cl_df1$snv.2,cl_df1$snv.3,cl_df1$snv.4,sep=";"), ';NA',2)[,1])
rm(cl, cl_df0, cl_df1)
m.cl <- merge(all.vc[,c('cell','mutation',"mutsites")],cl_df, by="mutsites")

### Referece mutation sites for the genes of interest in GWAS Catalog
combine_files <- function(dir){
  read_GWAS <- function(file){
    gene <- read.csv(file, sep = '\t')
    gene[,c('chr','pos')] <- str_split_fixed(gene$locations, ':', 2)
    gene[,c('rs','risk_allele')] <- str_split_fixed(gene$riskAllele, '-', 2)
    mut.gene <- paste(gene$chr, gene$pos, sep = '_')
    mut.gene <- paste('chr', mut.gene, sep = '')
    gene$mutsites <- mut.gene
    return(gene)
  }
  
  list_files <- list.files(path = dir, full.names = TRUE, include.dirs=FALSE)
  nrfiles <- length(list_files)
  print(paste('Number of files:',nrfiles))
  df_list = list()
  i=1
  for (file in 1:length(list_files)) {
    file = list_files[i]
    df <- read_GWAS(file)
    df_list[[i]] <- as.data.frame(df)
    i=i+1
  }
  df_full <- bind_rows(df_list)
  return(df_full)
}

gwas <- combine_files("D:/OneDrive - KU Leuven/phd/NMTseq/metadata/ref_variant/dbSNP")
gwas_df0 <- gwas[which(gwas$mutsites %in% var_sites),c('mutsites','riskAllele','traitName')]
gwas_df1 <- data.frame(mutsites = gwas_df0$mutsites, snv = paste("GWAS",gwas_df0$riskAllele,gwas_df0$traitName,sep = ':'))
gwas_df1 <- gwas_df1%>%
  group_by(mutsites) %>%
  mutate(occurrence = row_number()) %>%
  ungroup()
gwas_df1 <- as.data.frame(gwas_df1)
gwas_df1 <- reshape(gwas_df1, idvar="mutsites", timevar = "occurrence", direction="wide")
gwas_df1_snv <- gwas_df1[,-1]
snv = gwas_df1_snv[,1]
for(i in 2:ncol(gwas_df1_snv)){
  snv = paste(snv, gwas_df1[,i],sep=";")
}
gwas_df <- data.frame(mutsites = gwas_df1$mutsites ,snv = str_split_fixed(snv,';NA',2)[,1])
rm(gwas, gwas_df0, gwas_df1)
m.gwas <- merge(all.vc[,c('cell','mutation',"mutsites")],gwas_df, by="mutsites")

#### Combine all variant annotations from all the databases
m.all <- rbind(m.cosmic, m.tcga, m.cl, m.gwas)
m.all <- m.all[order(m.all$cell),]
m.all <- m.all[order(m.all$mutsites),]
m.all$SNV <- paste(m.all$mutation, m.all$snv, sep='. DB={')
m.all$SNV <- paste(m.all$SNV, '}')
head(m.all)
dim(m.all)
unique(m.all$cell)
unique(m.all$mutation)
write.csv(m.all[,c(2,5)], "D:/OneDrive - KU Leuven/phd/NMTseq/metadata/mutations_all_5genes.csv",row.names = FALSE)

m.all$mutation = paste('[',m.all$mutation, sep = '')
m.all$mutation = paste(m.all$mutation, ']:',sep = '')
m.all_collapsed <- m.all %>%
  group_by(mutation, cell) %>%
  summarize(snv = paste(snv, collapse = ";"))
head(m.all_collapsed)

m.all.new <- m.all_collapsed %>%
  group_by(cell)  %>%
  summarize(SNV_DB = paste0(mutation, snv, sep='', collapse = "\n")) 
m.all.new$SNV_DB[3]
m.all.new <- as.data.frame(m.all.new)
head(m.all.new)

### Integrate with depth of the loci obtained by samtools depth
# For all the sites with at least one mutation detected in one cell, samtools depth was used to compute all the depth of single sites
depth.all <- read.csv('D:/OneDrive - KU Leuven/phd/NMTseq/variant/VariantCalling/merged_depth.tsv', sep = '\t', header = FALSE)
depth.all[,c('cell','chr')]=str_split_fixed(depth.all$V1, ' ', 2)
depth.all <- depth.all[,c(4,5,2,3)]
colnames(depth.all) <- c('cell','chr','pos','depth')
depth.all$chrpos <- paste(depth.all$chr, depth.all$pos,sep='_')

# Wrapper for dataframe processing
get_depth_df <- function(mutations){
  if (length(mutations)==1){
    common_mut_df = all.vc
  } else {
    common_mut_df = all.vc[which(all.vc$mut %in% mutations),]
  }
  common_mut_df <- common_mut_df %>%
    group_by(mut) %>%
    summarize(cells = paste(cell, collapse = ","))%>%
    ungroup()
  common_mut_df <- as.data.frame(common_mut_df)
  common_mut_df$chr <- str_split_fixed(common_mut_df$mut,'_',3)[,1]
  common_mut_df$pos <- str_split_fixed(common_mut_df$mut,'_',3)[,2]
  common_mut_df$chrpos <- paste(common_mut_df$chr,common_mut_df$pos,sep='_')
  
  # Depth data
  depth.all.sub <- depth.all[which(paste(depth.all$chr, depth.all$pos) %in% paste(common_mut_df$chr,common_mut_df$pos)),]
  
  # Merge depth.all.sub to common_mut_df
  # Merge the two data frames
  merged_df <- common_mut_df %>%
    left_join(depth.all.sub, by = "chrpos")
  
  # Pivot the data frame to the desired format
  pivoted_df <- merged_df %>%
    select(chrpos, mut, cells, cell, depth) %>%
    pivot_wider(names_from = cell, values_from = depth, values_fill = 0)
  
  # Replace NA values with 0
  pivoted_df[is.na(pivoted_df)] <- 0
  pivoted_df <- as.data.frame(pivoted_df)
  
  # View the resulting data frame
  pivoted_df_cells <- sort(grep('GC',colnames(pivoted_df),value = TRUE))
  sorted_col_names <- c('chrpos','mut','cells',pivoted_df_cells)
  pivoted_df <- pivoted_df[,sorted_col_names]
  
  # Include metadata
  meta = 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/RNA_Meta.csv'
  make_df <- function(file){
    df = read.csv(file)
    colnames(df)[1] = 'cell'
    return(df)
  }
  meta = make_df(meta)[,c(1,9,10,11)]
  #meta$ident = sub("^(GC\\d+).*", "\\1", meta$cell)
  
  # Add summaries like average values
  pivoted_df$avg_all <- rowMeans(pivoted_df[,pivoted_df_cells])
  cells_tumor_CD45neg <- meta[intersect(which(meta$Tumor=='Y'),which(meta$CD45=='Negative')),'cell']
  cells_tumor_CD45pos <- meta[intersect(which(meta$Tumor=='Y'),which(meta$CD45=='Positive')),'cell']
  cells_normal_CD45neg <- meta[intersect(which(meta$Tumor=='N'),which(meta$CD45=='Negative')),'cell']
  cells_normal_CD45pos <- meta[intersect(which(meta$Tumor=='N'),which(meta$CD45=='Positive')),'cell']
  depth_df_tumor_CD45neg = pivoted_df[,which(sorted_col_names %in% cells_tumor_CD45neg)]
  depth_df_tumor_CD45pos = pivoted_df[,which(sorted_col_names %in% cells_tumor_CD45pos)]
  depth_df_normal_CD45neg = pivoted_df[,which(sorted_col_names %in% cells_normal_CD45neg)]
  depth_df_normal_CD45pos = pivoted_df[,which(sorted_col_names %in% cells_normal_CD45pos)]
  pivoted_df$avg_depth_tumor_CD45neg = rowMeans(depth_df_tumor_CD45neg)
  pivoted_df$avg_depth_tumor_CD45pos = rowMeans(depth_df_tumor_CD45pos)
  pivoted_df$avg_depth_normal_CD45neg = rowMeans(depth_df_normal_CD45neg)
  pivoted_df$avg_depth_normal_CD45pos = rowMeans(depth_df_normal_CD45pos)
  
  # Seperate the cell names in the "cells" column
  pivoted_df$var_cells <- strsplit(pivoted_df$cells, ",")
  
  # Add notes on normal/tumor, CD45+/-, having coverage in normal tissue or not
  pivoted_df$nr_cells_carrying_this_variant = NA
  pivoted_df$Note = ""
  for (i in 1:nrow(pivoted_df)){
    pivoted_df$nr_cells_carrying_this_variant[i] <- length(pivoted_df$var_cells[i][[1]])
    
    nr_CD45pos <- length(intersect(pivoted_df$var_cells[i][[1]], c(cells_tumor_CD45pos, cells_normal_CD45pos)))
    if (nr_CD45pos>0){
      pivoted_df$Note[i] = paste('variant in CD45+',pivoted_df$Note[i], sep=';')
    }
    
    nr_normal <- length(intersect(pivoted_df$var_cells[i][[1]], c(cells_normal_CD45neg, cells_normal_CD45pos)))
    if (nr_normal>0){
      pivoted_df$Note[i] = paste('variant in normal tissue',pivoted_df$Note[i], sep=';')
    }
    
    sum_depth_normal <- pivoted_df$avg_depth_normal_CD45neg[i]+pivoted_df$avg_depth_normal_CD45pos[i]
    if (sum_depth_normal==0){
      pivoted_df$Note[i] = paste('no_coverage_in_normal_tissue',pivoted_df$Note[i], sep=';')
    } else if (sum_depth_normal<0.05){
      pivoted_df$Note[i] = paste('low_coverage_in_normal_tissue',pivoted_df$Note[i], sep=';')
    }
    # sum_depth_CD45pos <- pivoted_df$avg_depth_normal_CD45pos[i]+pivoted_df$avg_depth_tumor_CD45pos[i]
    # if (sum_depth_CD45pos==0){
    #   pivoted_df$Note[i] = paste('no_coverage_in_CD45pos',pivoted_df$Note[i], sep=';')
    # } else if (sum_depth_CD45pos<0.05){
    #   pivoted_df$Note[i] = paste('low_coverage_in_CD45pos',pivoted_df$Note[i], sep=';')
    # }
  }
  pivoted_df <- pivoted_df[,-which(colnames(pivoted_df)=='var_cells')]
  pivoted_df <- pivoted_df %>% relocate(chrpos,mut,cells,nr_cells_carrying_this_variant,avg_all,avg_depth_tumor_CD45pos,avg_depth_tumor_CD45neg,avg_depth_normal_CD45neg,avg_depth_normal_CD45pos,Note)
  colnames(pivoted_df)[1:3] <- c('location','variant','cells_carrying_this_variant')
  print("Nr of variants detected:")
  print(nrow(pivoted_df))
  print("Nr of variants of interest:")
  print(length(which(pivoted_df$Note=="")))
  # Add the gene names to the df
  pivoted_df$chr = str_split_fixed(pivoted_df$location,'_',2)[,1]
  pivoted_df$pos = as.numeric(str_split_fixed(pivoted_df$location,'_',2)[,2])
  pivoted_df$gene <- NA
  pivoted_df$gene <- ifelse(pivoted_df$chr == TP53[1] & pivoted_df$pos >= TP53[2] & pivoted_df$pos <= TP53[3], "TP53", pivoted_df$gene)
  pivoted_df$gene <- ifelse(pivoted_df$chr == KRAS[1] & pivoted_df$pos >= KRAS[2] & pivoted_df$pos <= KRAS[3], "KRAS", pivoted_df$gene)
  pivoted_df$gene <- ifelse(pivoted_df$chr == EGFR[1] & pivoted_df$pos >= EGFR[2] & pivoted_df$pos <= EGFR[3], "EGFR", pivoted_df$gene)
  pivoted_df$gene <- ifelse(pivoted_df$chr == ERBB2[1] & pivoted_df$pos >= ERBB2[2] & pivoted_df$pos <= ERBB2[3], "ERBB2", pivoted_df$gene)
  pivoted_df$gene <- ifelse(pivoted_df$chr == BRAF[1] & pivoted_df$pos >= BRAF[2] & pivoted_df$pos <= BRAF[3], "BRAF", pivoted_df$gene)
  pivoted_df <- pivoted_df[,!(names(pivoted_df) %in% c('chr','pos'))]
  pivoted_df <- pivoted_df %>% select(location, variant, gene, everything())
  
  return(pivoted_df)
}

# Same mutations in different cells in one patients (frequency>1)
all.vc.tb <- as.data.frame(table(all.vc[,'mut']))
dim(all.vc.tb)
head(all.vc.tb)
common_mut <- all.vc.tb[which(all.vc.tb$Freq>1),'Var1']
length(common_mut)

# Mutations that appeared at least twice (frequency>=2)
depth_df <- get_depth_df(common_mut)

# Get the mutations that appear >=2 times in one patient 
# (remove the ones that appear across patients but only once in each patient)
head(depth_df)
colnames(depth_df)[12:ncol(depth_df)]
# Get the cell GC codes
#ident <- sub("^(GC\\d+).*", "\\1", colnames(depth_df)[12:ncol(depth_df)])

# Get the patients corresponding to the cells
patient = 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/data sheet _GC code & patient number.csv'
patient = read.csv(patient)
colnames(patient) <- c('ident','patient')
head(patient)
p <- character(length(ident))
for (i in seq_along(ident)) {
  p[i] <- patient[which(patient$ident == ident[i]), 'patient']
}
p

# Get the indices for which there is more than 1 read in a cell for a position
#depth_df_c <- depth_df[,(12:ncol(depth_df))]
#depth_df_c_boo <- depth_df_c>0
# Apply the function to each row using apply
#indices_list <- apply(depth_df_c_boo, 1, function(row) {which(row)})
# Convert the result to a list of numeric vectors
#indices_list <- lapply(indices_list, as.numeric)

cells <- strsplit(depth_df$cells, ",")

# Function to extract prefixes
extract_prefixes <- function(lst) {
  sapply(lst, function(x) {
    elements <- unlist(strsplit(x, "_"))
    prefixes <- elements[1]
    return(prefixes)
  })
}

# Apply the function to each element in the list
cells_gc <- lapply(cells, extract_prefixes)
ls_patient = list()
with_multiple_occ_same_patient = c()

for (i in 1:length(cells_gc)){
  ident_ls = cells_gc[i][[1]]
  n_cells <- length(ident_ls)
  p <- c()
  for (j in seq_along(ident_ls)) {
    p = c(p, patient[which(patient$ident == ident_ls[j]), 'patient'])
  }
  # Get if there is one patient that appeared twice (more than one cell that show variant in the same patient)
  
  more_than_once_patient <- p[duplicated(p)]
  if (length(more_than_once_patient)>0){
    with_multiple_occ_same_patient = c(with_multiple_occ_same_patient, TRUE)
  } else {
    with_multiple_occ_same_patient = c(with_multiple_occ_same_patient, FALSE)
  }
  
  ls_patient[[i]] = p
  i=i+1
}

depth_df$with_multiple_occ_same_patient = with_multiple_occ_same_patient
depth_df_reoccuring <- depth_df[which(depth_df$with_multiple_occ_same_patient),]
write.csv(depth_df_reoccuring, 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/common_mut_depth.csv', row.names = FALSE)

# Subset to the varaints of interest (exclude the shared variant that are shared with non-tumor/CD45-)
depth_df_reoccuring_sub <- depth_df_reoccuring[which(depth_df_reoccuring$Note==''),]
# Get a compiled list of cells that passed filtering
ls_all_cells_shared_variant <- unlist(strsplit(depth_df_reoccuring_sub$cells_carrying_this_variant, ","))
df_all_cells_shared_variant <- as.data.frame(table(ls_all_cells_shared_variant))
colnames(df_all_cells_shared_variant) <- c('cell', 'nr_shared_variants_in_cell')
df_all_cells_shared_variant
# Combine with QC result (so that we can filter out the cells failing QC)
qc = 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/qc_202310.csv'
qc = read.csv(qc)
colnames(qc)[1] <- 'cell'
df_all_cells_shared_variant <- merge(df_all_cells_shared_variant, qc, by='cell')
df_all_cells_shared_variant <- df_all_cells_shared_variant[order(-df_all_cells_shared_variant$nr_shared_variants_in_cell), ]
write.csv(df_all_cells_shared_variant, 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/cells_with_sv.csv', row.names = FALSE)

# Mutations that appeared at least once (frequency>=1)
depth_df0 <- get_depth_df('all')
print(paste("Before filtering: Nr of variants detected in at least one cell:",nrow(depth_df0)))
write.csv(depth_df0, 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/all_mut_depth.csv', row.names = FALSE)

# Select the sub data with the variants of interest (=variants that are not filtered out by depth or prescence in CD45+/non-tumor)
depth_df_voe <- depth_df0[which(depth_df0$location %in% depth_df0[which(depth_df0$Note==""),'location']),]
print(paste("After filtering: Nr of variants of interest:",nrow(depth_df_voe)))

# Check the number of occurrence of the cells in the variants of interest
all_cells_with_voe <- unlist(strsplit(depth_df_voe$cells_carrying_this_variant, ","))
all_cells_with_voe <- as.data.frame(table(all_cells_with_voe)) %>% arrange(desc(Freq))
colnames(all_cells_with_voe) <- c("cell","nr_VOE")
head(all_cells_with_voe)
print(paste("Nr of cells with at least one variant of interest:",nrow(all_cells_with_voe)))

# Merge with QC results
qc = 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/qc_202310.csv'
qc = read.csv(qc)
colnames(qc)[1] <- 'cell'

all_cells_with_voe <- merge(all_cells_with_voe, qc, by='cell')
write.csv(all_cells_with_voe, 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/all_cells_with_voe.csv', row.names = FALSE)

## Now get the VOEs again, but allow "variant in normal tissue" to occur in notes
depth_df_voe <- depth_df0[-grep('variant in CD45+', depth_df0$Note),]
depth_df_voe <- depth_df_voe[-grep('no_coverage_in_CD45pos', depth_df_voe$Note),]
depth_df_voe <- depth_df_voe[-grep('low_coverage_in_CD45pos', depth_df_voe$Note),]
print(paste("After filtering: Nr of variants of interest:",nrow(depth_df_voe)))
voe_loc <- depth_df_voe$location

# Find out which variants of interest are annotated in the databases
m.all$annotated_loc <- paste(str_split_fixed(m.all$mutation,':',6)[,2], str_split_fixed(m.all$mutation,':',6)[,3], sep='_')
voe_loc[which(voe_loc %in% m.all$annotated_loc)]
voe_df <- m.all[which(m.all$annotated_loc %in% voe_loc[which(voe_loc %in% m.all$annotated_loc)]),]
voe_df_collapsed <- voe_df %>%
  group_by(mutsites, cell) %>%
  summarize(SNV_DB = paste(snv, collapse = ";"))
voe_df_final <- voe_df_collapsed %>%
  group_by(mutsites, SNV_DB)  %>%
  summarize(cells_carrying_var = paste(cell, collapse = ",")) 
head(voe_df_final)
unique(voe_df_final$cells_carrying_var)
print(paste("After filtering: Nr of database-annotated variants:",length(unique(voe_df_final$mutsites))))
voe_df_final$SNV_DB[3]
write.csv(voe_df_final, 'D:/OneDrive - KU Leuven/phd/NMTseq/metadata/voe_DB_annotated2.csv', row.names = FALSE)


# --------------------------------------------------------------------------------
#### Coverage
#### Number of mutations detected / total "site coverage"  within the 3 genes of interest
### Process the coverage obtained by samtools coverage
cov_df <- read.csv('D:/OneDrive - KU Leuven/phd/NMTseq/rna_20231005/merged_coverage.tsv', sep='\t', row.names = 1, header = FALSE)
cov_df <- cov_df[,-c(1,2,3,10,11,12,19,20,21)]
cov_cols <- c('numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq')
colnames(cov_df)[c(1:6)] = paste(cov_cols, 'TP53', sep = '_')
colnames(cov_df)[c(7:12)] = paste(cov_cols, 'KRAS', sep = '_')
colnames(cov_df)[c(13:18)] = paste(cov_cols, 'EGFR', sep = '_')
head(cov_df)
cov_df$cell <- rownames(cov_df)

all.vc.sum <- all.vc %>%
  group_by(cell, gene) %>%
  summarize(sum_coverage_var = sum(coverage_var))
all.vc.sum <- as.data.frame(all.vc.sum)
all.vc.sum
all.vc.sum <- all.vc.sum %>%
  pivot_wider(names_from = gene, values_from = sum_coverage_var, 
              names_prefix = "cov_var_", values_fill = 0)
all.vc.sum <- as.data.frame(all.vc.sum)

# Generate a df with coverage data and the varaints calling data
merged_coverage <- merge(cov_df, all.vc.sum, by='cell', all=TRUE)
merged_coverage <- merged_coverage[,c(1,2,3,21,4,5,6,7,8,9,22,10,11,12,13,14,15,20,16,17,18,19)]
head(merged_coverage)
# The NAs in perc_var simply means no variants detected in this gene; thus we replace them with zeros
merged_coverage[is.na(merged_coverage)] <- 0
dim(merged_coverage)
# High coverage but no variant detected -> high confidence that it's not a cancer cell
genelenTP53 = 7687538-7661779+1
genelenKRAS = 25250936-25205246+1
genelenEGFR = 55211628-55019017+1
merged_coverage$perc_var_TP53 <- as.numeric(merged_coverage$cov_var_TP53/(genelenTP53*(merged_coverage$meandepth_TP53)))
merged_coverage$perc_var_KRAS <- as.numeric(merged_coverage$cov_var_KRAS/(genelenKRAS*(merged_coverage$meandepth_KRAS)))
merged_coverage$perc_var_EGFR <- as.numeric(merged_coverage$cov_var_EGFR/(genelenEGFR*(merged_coverage$meandepth_EGFR)))
# Take a look at the results of the percentage of variants discovered
boxplot(merged_coverage$perc_var_TP53)
summary(merged_coverage$perc_var_TP53)
summary(merged_coverage$perc_var_KRAS)
summary(merged_coverage$perc_var_EGFR)
summary(merged_coverage$coverage_TP53)
summary(merged_coverage$coverage_KRAS)
summary(merged_coverage$coverage_EGFR)


merged_coverage <- merged_coverage[,c(1,2,3,21,23,4,5,6,7,8,9,22,24,10,11,12,13,14,15,20,25,16,17,18,19)]
head(merged_coverage)
# Combine with the variant annotation result
merged.all <- merge(merged_coverage, m.all.new, by='cell', all=TRUE)

# Identify the cells with adequate coverage for all 3 genes but no variant detected
merged.all$Note_var <- ""
sel = which(merged.all$perc_var_TP53==0)
merged.all$Note_var[sel] = paste(merged.all$Note_var[sel],"No_variants_TP53",sep=";")
sel = which(merged.all$perc_var_KRAS==0)
merged.all$Note_var[sel] = paste(merged.all$Note_var[sel],"No_variants_KRAS",sep=";")
sel = which(merged.all$perc_var_EGFR==0)
merged.all$Note_var[sel] = paste(merged.all$Note_var[sel],"No_variants_EGFR",sep=";")
merged.all$Note_var <- sub("^;|;$", "", merged.all$Note_var)
#my_intersect <- intersect(which(merged.all$perc_var_TP53==0),which(merged.all$perc_var_KRAS==0))
#my_intersect <- intersect(my_intersect,which(merged.all$perc_var_EGFR==0))
#merged.all$Note[my_intersect] <- "No_variants_in_all_3_genes"

# Identify the cells with No_coverage_in_all_3_genes
merged.all$Note_cov <- ""
sel = which(merged.all$covbases_TP53==0)
merged.all$Note_cov[sel] = paste(merged.all$Note_cov[sel],"No_cov_TP53",sep=";")
sel = which(merged.all$covbases_KRAS==0)
merged.all$Note_cov[sel] = paste(merged.all$Note_cov[sel],"No_cov_KRAS",sep=";")
sel = which(merged.all$covbases_EGFR==0)
merged.all$Note_cov[sel] = paste(merged.all$Note_cov[sel],"No_cov_EGFR",sep=";")
merged.all$Note_cov <- sub("^;|;$", "", merged.all$Note_cov)
#my_intersect <- intersect(which(is.na(merged.all$perc_var_TP53)),which(is.na(merged.all$perc_var_KRAS)))
#my_intersect <- intersect(my_intersect,which(is.na(merged.all$perc_var_EGFR)))
#merged.all$Note[my_intersect] <- "No_coverage_in_all_3_genes"

# Identify the cells with high percentage of variants
merged.all$Note_mut <- ""
sel = which(merged.all$perc_var_TP53>0.1)
merged.all$Note_mut[sel] = paste(merged.all$Note_mut[sel],"High_var_TP53",sep=";")
sel = which(merged.all$perc_var_KRAS>0.1)
merged.all$Note_mut[sel] = paste(merged.all$Note_mut[sel],"High_var_KRAS",sep=";")
sel = which(merged.all$perc_var_EGFR>0.1)
merged.all$Note_mut[sel] = paste(merged.all$Note_mut[sel],"High_var_EGFR",sep=";")
merged.all$Note_mut <- sub("^;|;$", "", merged.all$Note_mut)

merged.all <- merged.all %>% relocate(Note)
merged.all <- merged.all %>% relocate(SNV_DB)
merged.all <- merged.all %>% relocate(cell,SNV_DB,Note_mut,Note_var,Note_cov)
write.csv(merged.all, "D:/OneDrive - KU Leuven/phd/NMTseq/metadata/merged_coverage_variants.csv", row.names = FALSE)

length(which(merged.all$covbases_TP53==0))/nrow(merged.all)


### Plotting the coverage
colnames(merged.all)
library(raincloudplots)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)

plot_df <- merged.all[,c('coverage_TP53','coverage_KRAS','coverage_EGFR')]
rownames(plot_df) <- merged.all$cell


plot_df <- gather(plot_df, key = gene, value = coverage)
plot_df$gene <- sub("coverage_", "", plot_df$gene)
head(plot_df)

ggplot(plot_df, aes(x = gene, y = coverage, fill = gene))+ 
  ylim(0,1.5)+ # Limits of y axis
  geom_flat_violin(aes(fill = gene), position = position_nudge(x = 0, y = 0), 
                   adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_boxplot(alpha = .5, colour = "black", position = position_nudge(x = -0.15, y = 0), 
               outlier.size=0.03, lwd=0.2, width=0.25)+ #, outlier.shape = NA
  theme(axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border = element_blank(),panel.background = element_blank())

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "area",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}


GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
          },
  
  draw_group = function(data, panel_scales, coord) {
    # Find the points for the line to go all the way around
    data <- transform(data, xminv = x,
                      xmaxv = x + violinwidth * (xmax - x))
    
    # Make sure it's sorted properly to draw the outline
    newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                     plyr::arrange(transform(data, x = xmaxv), -y))
    
    # Close the polygon: set first and last point the same
    # Needed for coord_polar and such
    newdata <- rbind(newdata, newdata[1,])
    
    ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
  },
  
  draw_key = draw_key_polygon,
  
  default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                    alpha = NA, linetype = "solid"),
  
  required_aes = c("x", "y")
  )

### Alternatively, use TCGAbiolinks to access TCGA data
# library(TCGAbiolinks)
# query <- GDCquery(
#   project = c("TCGA-LUAD","APOLLO-LUAD", "CGCI-HTMCP-LC"), 
#   data.category = "Simple Nucleotide Variation", 
#   access = "open",
#   data.type = "Masked Somatic Mutation", 
#   workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
# )
# GDCdownload(query)
# data <- GDCprepare(query)
# head(data)
# data$mutsite <- paste(data$Chromosome, data$Start_Position, sep = '_')
# 
# GDC_goe <- as.data.frame(data[data$Hugo_Symbol %in% c("TP53","KRAS","EGFR"),])
# rm(query)
# rm(data)
# 
# all.vc[which(var_sites %in% GDC_goe$mutsite),]
# GDC_goe[which(GDC_goe$mutsite %in% var_sites),]

# Save the R image
save.image("vcf_20231017.RData")
