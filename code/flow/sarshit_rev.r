library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(anndata)
library(gprofiler2)
library(Biostrings)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]
sars <- args[2]

# Activity scores (output from miReact)
b <- readRDS(celltypes) 
sequences <- readDNAStringSet(sars)
rseq <- reverseComplement(sequences$`hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124`)

subs_dfs_paths <- gsub("annot4_", "annot5_", celltypes)


###### Does a kmer hit the sarscov2 genome?
matches <- logical(nrow(b))

# Loop through each row of selected_rows
for (i in 1:nrow(b)) {
  # Extract the sequence (row name)
  sequence <- rownames(b)[i]
  
  # Check if the sequence exists within the bigger sequence
  if (grepl(sequence, rseq, fixed = TRUE)) {
    matches[i] <- TRUE  # Set match to TRUE if found
  } else {
    matches[i] <- FALSE  # Set match to FALSE if not found
  }
}

# Add the matches vector as a new column named 'sarshit' in selected_rows
b$sarshit_rev <- matches

######## How many times does a sequence hit?
counts <- numeric(nrow(b))

# Loop through each row of selected_rows
for (i in 1:nrow(b)) {
  # Extract the sequence (row name)
  sequence <- rownames(b)[i]
  
  # Count the number of occurrences of the sequence in the bigger sequence
  match_positions <- gregexpr(sequence, rseq, fixed = TRUE)
  count <- sum(unlist(match_positions) >= 0)  # Count non-negative match positions
  
  # Store the count in the counts vector
  counts[i] <- count
}

# Add the counts vector as a new column named 'sarshitnumber' in selected_rows
b$sarshitnumber_rev <- counts

saveRDS(b, file = subs_dfs_paths)
