library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(anndata)
library(gprofiler2)
library(Biostrings)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]
sars <- args[2]

# Activity scores (output from miReact)
b <- readRDS(celltypes) 

# direction 5'--> 3'by consensus
# SARS-CoV-2 genome is +ssRNA virus, meaning that it's genome is encoded 3' --> 5'
sequences <- readDNAStringSet(sars) 
subs_dfs_paths <- gsub("annot3_", "annot4_", celltypes)


###### Does a kmer hit the sarscov2 genome?
matches <- logical(nrow(b))

# Loop through each row of selected_rows
for (i in 1:nrow(b)) {
  # Extract the sequence (row name)
  sequence <- rownames(b)[i]
  
  # Check if the sequence exists within the bigger sequence
  if (grepl(sequence, sequences$'hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124', fixed = TRUE)) {
    matches[i] <- TRUE  # Set match to TRUE if found
  } else {
    matches[i] <- FALSE  # Set match to FALSE if not found
  }
}

# Add the matches vector as a new column named 'sarshit' in selected_rows
b$sarshit <- matches

######## How many times does a sequence hit?
counts <- numeric(nrow(b))

# Loop through each row of selected_rows
for (i in 1:nrow(b)) {
  # Extract the sequence (row name)
  sequence <- rownames(b)[i]
  
  # Count the number of occurrences of the sequence in the bigger sequence
  match_positions <- gregexpr(sequence, sequences$'hCoV-19/Wuhan/WIV04/2019|EPI_ISL_402124', fixed = TRUE)
  count <- sum(unlist(match_positions) >= 0)  # Count non-negative match positions
  
  # Store the count in the counts vector
  counts[i] <- count
}

# Add the counts vector as a new column named 'sarshitnumber' in selected_rows
b$sarshitnumber <- counts

saveRDS(b, file = subs_dfs_paths)
