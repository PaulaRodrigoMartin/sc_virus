library(RColorBrewer)
library(dplyr)
library(ggvenn)

args<-commandArgs(trailingOnly = TRUE) 

sec <- args[1]
cil <- args[2]
sq <- args[3]
annot <- args[4]
mir <- args[5]


sec <- readRDS(sec)
cil <- readRDS(cil)
sq <- readRDS(sq)
annot <- read.csv(annot, header = TRUE, row.names = 1)
mir <- readRDS(mir)

file_names <- list(sec,cil,sq)
all_rownames <- character()
# Loop through each file
for (file in file_names) {
    all_rownames <- c(all_rownames, rownames(file))
}

# Remove duplicates
unique_rownames <- unique(all_rownames)

celltypes <- c("Secretory","Ciliated","Squamous")

# Create an empty data frame with rows for each string and columns for each celltype
df <- data.frame(matrix(0, ncol = length(celltypes), nrow = length(unique_rownames)))
colnames(df) <- celltypes
rownames(df) <- unique_rownames

for (col in seq_len(ncol(df))) {
  for (row in rownames(df)) {
    if (row %in% rownames(file_names[[col]])) {
      df[row, col] <- TRUE}}}

print(head(df))

conqmir_venn <- df
myCol <- brewer.pal(11, "Set3")
ggvenn(
  conqmir_venn[,cytok], 
  fill_color = myCol,
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = F
  )