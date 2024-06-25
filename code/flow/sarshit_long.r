library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(anndata)
library(gprofiler2)
library(Biostrings)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]
sars <- args[2]

# Activity scores (output from miReact)
b <- readRDS(celltypes) 
sequences <- readDNAStringSet(sars) # 5'--> 3'
subs_dfs_paths <- gsub("longermer_downreg_filtered_", "annot_longmer_", celltypes)

#################33
## Functions
###################

#Create PFM, regular expression and alignment to 3'UTRs
# Align longer kmers as regular expressions accounting for substitutions to the 3'UTR sequences to further filtering.

list_names <- colnames(long_and_7mers)
b <- character(length = length(list_names))
a <- character(length = length(list_names))#init empty vector

row_name <- list()
for (i in seq_along(list_names)) {
  row_name <- append(row_name,sub(" .*", "", list_names[i]))
}
row_name <- unlist(row_name)

# work with matrix because in a dataframe we cannot have duplicated rownames
matrix_data <- matrix(
  data = NA_character_,
  ncol = 4,  
  nrow = length(row_name),  
  dimnames = list(row_name, c("regexp", "celltype", "align_UTR", "n_UTR_aligned"))
)

## BUILD REGEXP, ALIGN IT TO 3'UTRs, 
for (i in seq_along(list_names)) {
  b <- sub(" .*", "", list_names[i]) #extract the longmer
  cellt <- sub(".*\\s", "", list_names[i]) #extract the celltype
  a <- long_and_7mers[list_names[i]] #extract all the 7mers that have passed the filters and conform the longermer
  msa_result <- msa(as.character(a[[1]]), type="dna") ## kmers that conform the longermer
  t <- Biostrings::consensusMatrix(DNAStringSet(msa_result), baseOnly=TRUE)
  t <- t[-nrow(t), ]

  final_list <- list()
  for (col in 1:ncol(t)) {
    col_list <- create_regexp(t[, col]) # also create_regexp_GUwpbble()
    final_list <- c(final_list, col_list)
  }

  # Convert the final list to a single string --> regular expression
  regexpression <- paste(final_list, collapse = "")
  matrix_data[i, "regexp"] <- regexpression
  matrix_data[i, "celltype"] <- cellt
  # Align the regexpression to the 3'UTRs
  x <- align.motif(regexpression, seqlist, lf = 3, rf = 3)
  if (!is.null(x)){
    matrix_data[i, "align_UTR"] <- T
    matrix_data[i, "n_UTR_aligned"] <- as.numeric(length(x))
  } else {
    matrix_data[i, "align_UTR"] <- F
    matrix_data[i, "n_UTR_aligned"] <- 0
  }
}

# Save matrix_data --> shows the initial longmer, the updated longmers' regexpression and how many times it has aligned to some 3'UTR
write.csv(matrix_data, file = "res/final_longmers.csv", row.names = T, col.names = T)


#Build logos for visualization
pattern_nt <- "[A-Za-z]"
pattern_sub <- "\\|"
plot_list <- list()
for (i in 1:nrow(matrix_data)) {
  if (matrix_data[i, "align_UTR"] == T & as.numeric(matrix_data[i, "n_UTR_aligned"]) < 500) {
    x <- align.motif_updated(matrix_data[i,"regexp"], seqlist, lf = 10, rf = 10)
    if (length(x)>1){
      letter_count <- str_count(matrix_data[i,"regexp"], pattern_nt)
      pipe_count <- str_count(matrix_data[i,"regexp"], pattern_sub)
      gg <- ggplot() + 
      geom_logo(x, method = 'prob') + 
      annotate('segment', x = 10, xend= ((letter_count-pipe_count)+9), y=0, yend=0, size=1) + 
      theme_logo()+
      ggtitle(paste0(matrix_data[i,"regexp"], "  in  " , matrix_data[i,"celltype"], "  found  ", matrix_data[i,"n_UTR_aligned"], " times")) + theme(plot.title = element_text(hjust = 0.5))
      plot_list[[i]] <- gg
    }
  }
}

# Eliminate the plots that are empty 
plot_list <- Filter(Negate(is.null), plot_list)


# Save the logos

pdf("/results/plots/logos_good.pdf", height = 50, width =20) # adjust path to file, height and width as needed
grid.arrange(grobs = plot_list, ncol = 1)  # adjust the number of columns and rows as needed
dev.off()


###### Does a longmer hit the sarscov2 genome?
matches <- logical(nrow(b))


######################################
# transpose longermer df 
######################################

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
