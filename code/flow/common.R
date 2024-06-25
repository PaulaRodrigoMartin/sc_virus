library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(anndata)

args<-commandArgs(trailingOnly = TRUE) 

a <- args[1]
b <- args[2]
c <- args[3]
d <- args[4]
e <- args[5]
f <- args[6]
g <- args[7]
h <- args[8]
i <- args[9]
j <- args[10]
# k <- args[11]

subs_dfs_paths <- gsub("filtered_7mers_", "common_7mers_", a)
subs <- gsub("_Secretory", "", subs_dfs_paths)
print(subs)


a1 <- readRDS(a) 
b1 <- readRDS(b) 
c1 <- readRDS(c) 
d1 <- readRDS(d) 
e1 <- readRDS(e) 
f1 <- readRDS(f) 
g1 <- readRDS(g) 
h1 <- readRDS(h) 
i1 <- readRDS(i) 
j1 <- readRDS(j) 
# k1 <- readRDS(k)

file_names <- list(a1, b1, c1, d1, e1, f1, g1, h1, i1, j1)
all_rownames <- character()
# Loop through each file
for (file in file_names) {
    all_rownames <- c(all_rownames, rownames(file))
}

# Remove duplicates
unique_rownames <- unique(all_rownames)

celltypes <- c("Secretory","Mono","Macro","T","DC","Ciliated","Neu", "NK","Plasma","Squamous")

# Create an empty data frame with rows for each string and columns for each celltype
df <- data.frame(matrix(0, ncol = length(celltypes), nrow = length(unique_rownames)))
colnames(df) <- celltypes
rownames(df) <- unique_rownames

for (col in seq_len(ncol(df))) {
  for (row in rownames(df)) {
    if (row %in% rownames(file_names[[col]])) {
      df[row, col] <- df[row, col] + 1}}}
df$RowSum <- rowSums(df)

saveRDS(df, file = subs)