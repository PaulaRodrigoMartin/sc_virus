library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(anndata)
library(gprofiler2)
library(Biostrings)
library(gridExtra)
library(grid)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]

print(celltypes)


subs_dfs_paths_1 <- gsub("ranked_7mers_", "filtered_7mers_", celltypes)
subs_dfs_paths_2 <- gsub("ranked_7mers_", "filtered_table_", celltypes)
png2 <- gsub(".rds", ".png", subs_dfs_paths_2)

b <- readRDS(celltypes) 


###############################
############### Hypothesis 1
###############################
a <- b
a <- a[order(a$corr_viralload, decreasing = TRUE), ]

x <- length(which(a$is_mirnatarget[1:100] == T))/ nrow(a) *100
print(paste("miRNA enrichment in top 100 ranked by viral load: ", x))

#highly correlated with viral load
if (dim(a)[1]>300){
    a_100 <- a[1:300,]
} else {
    a_100 <- a
}


# Subset by most extreme log2fc
cand_plasm <- a_100[a_100$log2fc < summary(a_100$log2fc)[2] | a_100$log2fc > summary(a_100$log2fc)[5],]
print(paste("dimensions after filtering by logFC: ", dim(cand_plasm))) 

cand_plasm$logttest <- -log10(cand_plasm$ttest_pval_corrected)
cand_plasm <- cand_plasm[cand_plasm$logttest > summary(cand_plasm$logttest)[5],]

# x <- length(which(cand_plasm$is_mirnatarget[1:100] == T))/ nrow(cand_plasm) *100
# print(paste("miRNA enrichment in top 100 filtered by log2fc: ", x))

## doesnt hit in the sarscov2 genome
# cand_plasm_novir <- cand_plasm[cand_plasm$sarshitnumber == 0 | cand_plasm$sarshitnumber_rev == 0,]
print(paste("final dim after checking sars no hit: ", dim(cand_plasm))) 


d <- cand_plasm[,c(1, 8, 14, 23, 25, 26, 31, 33)]
j <- c("mir target", "3'UTR hits", "Effect size\n(log2FC)", "pval\n(infected vs uninfected)", "correlation coefficient \nwith viral load", "correlation pval", "n'hits in SARS-CoV-2 \ngenome (5'-3')", "n'hits in SARS-CoV-2 \ngenome (3'-5')")
colnames(d) <- j
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))

df <- data.frame(Column1 = c(1, 1),Column2 = c(1, 1))

if (nrow(d)>0){
    saveRDS(cand_plasm, file = subs_dfs_paths_1)
} else {
    saveRDS(df, file = subs_dfs_paths_1)
}

hei <- 50*nrow(d)
wid <- 200*ncol(d)


if (hei > 0 & wid > 0){
    png(png2, height = hei, width = wid)
    grid.table(d, theme=tt)
    dev.off()
} else {
    png(png2, height = 500, width = 500)
    grid.table(df, theme=tt)
    dev.off()
}

print("Hypothesis 1 built") 