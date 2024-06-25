library(dplyr)
library(anndata)
library(ggpubr)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]

# Activity scores (output from miReact)
cell <- read_h5ad(celltypes) 
# xm <- readRDS("C:/Users/paula/Documents/Aarhus master/msc_thesis/workflow/data/hs.seqXmot.counts.utr3_mrs_7mer.rds")
subs_dfs_paths <- gsub("annot1_", "annot2_", celltypes)

#correlation value
cell$obs$corr_viralload <- NA_character_
cell$obs$corr_viralload <- apply(cell$X, 1, function(row) {
  cor(row[cell$var$Is_infected == TRUE], cell$var[cell$var$Is_infected == TRUE, "viral_load"], method = "spearman")
})

#correlation p-value
cell$obs$corr_viralload_pval <- apply(cell$X, 1, function(row) {
  hu <- cor.test(row[cell$var$Is_infected == TRUE], cell$var[cell$var$Is_infected == TRUE, "viral_load"], method = "spearman")
  return(hu$p.value)
})

write_h5ad(cell, file = subs_dfs_paths)

print(table(cell$var$Is_infected))
print(summary(cell$obs$corr_viralload))
print(summary(cell$obs$corr_viralload_pval))