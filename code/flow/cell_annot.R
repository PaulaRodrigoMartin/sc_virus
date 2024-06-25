library(dplyr)
library(anndata)
library(gprofiler2)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]
print(paste0("printargs2: ",args[2]))
# print(paste0("printargs12: ",args[12]))
target <- args[2]

# Activity scores (output from miReact)
cell <- read_h5ad(celltypes) 
xm <- readRDS(target)
# xm <- readRDS("C:/Users/paula/Documents/Aarhus master/msc_thesis/workflow/data/hs.seqXmot.counts.utr3_mrs_7mer.rds")
subs_dfs_paths <- gsub("df_", "annot1_", celltypes)

# calculate median activity healthy and infected 
cell$obs$healthy_med_act <- apply(cell$X[, !cell$var$Is_infected], 1, median)
cell$obs$infected_med_act <- apply(cell$X[, cell$var$Is_infected], 1, median)

# calculate mean activity healthy and infected
cell$obs$healthy_mean_act <- apply(cell$X[, !cell$var$Is_infected], 1, mean)
cell$obs$infected_mean_act <- apply(cell$X[, cell$var$Is_infected], 1, mean)
cell$obs$mean_act <- apply(cell$X, 1, mean)
cell$obs$log2fc <- log2((cell$obs$infected_mean_act + 1e-10)/(cell$obs$healthy_mean_act + 1e-10))

# calculate sd
cell$obs$healthy_sd_act <- apply(cell$X[, !cell$var$Is_infected], 1, sd)
cell$obs$infected_sd_act <- apply(cell$X[, cell$var$Is_infected], 1, sd)
cell$obs$sd_act <- apply(cell$X, 1, sd)

# calculate inf - hea values
cell$obs$median_diff <- cell$obs$infected_med_act - cell$obs$healthy_med_act
cell$obs$mean_diff <- cell$obs$infected_mean_act - cell$obs$healthy_mean_act

# WRS + BH corrected pvals
cell$obs$WRS_pval <- rep(NA, 16384)
# Perform Wilcoxon rank sum test between infected/healthy and store the pvalues
p_values <- apply(cell$X, 1, function(row){
  wil <-  wilcox.test(row[cell$var$Is_infected], row[!cell$var$Is_infected], conf.level = 0.99, alternative = "two.sided")
  return(wil$p.value)
})
cell$obs$WRS_pval <- p_values
cell$obs$WRS_pval_corrected <- p.adjust(cell$obs$WRS_pval, method = "BH") ##########

# Perform t test between infected/healthy (because the means are normally distributed) and store the pvalues
p_val <- apply(cell$X, 1, function(row){
  wil <-  t.test(row[cell$var$Is_infected], row[!cell$var$Is_infected], conf.level = 0.99, alternative = "two.sided")
  return(wil$p.value)
})
cell$obs$ttest_pval <- p_val
# Multiple testing correction for pvalues
cell$obs$ttest_pval_corrected <- p.adjust(cell$obs$ttest_pval, method = "BH") ##########

# Target genes
cell$obs$target_genes <- NA_character_
cell$obs$target_genes <- apply(xm > 0, 1, function(row) {
  hu <- colnames(xm)[row]
  hu <- paste(na.omit(hu), collapse = ", ")
  return(hu)
})
print("Target genes annotation --> done")

write_h5ad(cell, file = subs_dfs_paths)
# To then access the target names
# rbp$kmers <- lapply(strsplit(rbp$kmers, "|"), trimws)

