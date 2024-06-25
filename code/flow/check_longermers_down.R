library(dplyr)
library(Biostrings)
library(stringr)
library(remotes)
library(ggplot2)

# library(gridExtra)
library(msa)
library(tools)
# require(ggseqlogo)

args<-commandArgs(trailingOnly = TRUE) 

downreg <- args[1]
down <- args[2]
mat <- args[3]
annot <- args[4]

# Activity scores (output from miReact)
longmers <- readRDS(downreg)
cell <- readRDS(down)
# xm <- readRDS("C:/Users/paula/Documents/Aarhus master/msc_thesis/workflow/data/hs.seqXmot.counts.utr3_mrs_7mer.rds")
subs_dfs_paths_down <- gsub("longermer_downreg_", "longermer_downreg_filtered_", downreg)

# Get 7mers from where the longermers have been constructed
cell_down <- cell[cell$diffexpressed_corr == "DOWN_corrected",]
kmers <- rownames(cell_down)

# Read miRNA activity matrix output from miReact
all <- readRDS(mat)

# Read annotation table
ann <- read.csv(annot, header = TRUE, row.names = 1)


## Need to read downreg and check for any of the celltypes. The one that is present is going to be j
CELLTYPES <- c("Secretory","Mono","Macro","T","DC","Ciliated","Neu", "NK","Plasma","Squamous")

j <- NULL
# Loop through each element in CELLTYPES
for (celltype in CELLTYPES) {
  # Check if the celltype is a substring of downreg
  if (grepl(celltype, downreg)) {
    j <- celltype
  }
}

######################### 
# FUNCTIONS
#########################

# Function that removes kmers that have lower correlation than a certain threshold
filter_and_remove_low_occurrences <- function(cor_df, threshold = 0) {
  filtered_df <- subset(cor_df, Freq >= threshold)
  
  # Count occurrences of each variable in Var1
  var1_counts <- table(filtered_df$Var1)
  var2_counts <- table(filtered_df$Var2)
  # Identify variables that appear less than half of the unique values
  low_occurrence_vars <- names(var1_counts[var1_counts <= (length(unique(filtered_df$Var1)) / 2)])
  low_occurrence_vars_2 <- names(var2_counts[var2_counts <= (length(unique(filtered_df$Var2)) / 2)])
  # Remove rows where Var1 is a low-occurrence variable
  filtered_df <- subset(filtered_df, !(Var1 %in% low_occurrence_vars))
  filtered_df <- subset(filtered_df, !(Var2 %in% low_occurrence_vars_2))
  
  return(filtered_df)
}

align.motif <-
  function(motif,seqlist,lf,rf){ # include
    vec.res<-c()
    for(i in seqlist){
      vec <- gregexpr(motif,i$sequence)[[1]]
      if(vec[1]==-1)next
      for(k in vec){
        #	cat(substr(i$sequence,k-lf,k+nchar(motif)+rf-1),"\n")
        vec.res <- c(vec.res,substr(i$sequence,k-lf,k+nchar(motif)+rf-1))
      }}
    return(vec.res)
  }

align.motif_updated <- function(motif, seqlist, lf, rf) {
  vec.res <- c()
  for (i in seqlist) {
    vec <- gregexpr(motif, i$sequence)[[1]]
    if (vec[1] == -1) next
    for (k in vec) {

      # Calculate left and right flanking lengths
      left_flank <- max(lf - k, 0)
      right_flank <- max(k + nchar(motif) + rf - 1 - nchar(i$sequence), 0)
      
      # Extract the substring with padded flanks
      aligned_motif <- paste0(
        strrep(" ", left_flank),
        substr(i$sequence, max(k - lf + 1, 1), min(k + nchar(motif) + rf - 1, nchar(i$sequence))),
        strrep(" ", right_flank)
      )
      vec.res <- c(vec.res, aligned_motif)
    }
  }
  return(vec.res)
}

# Function that overlaps substrings that align partly
find_overlapping_substrings <- function(str1, str2) {
  max_overlap <- 0
  overlap <- min(nchar(str1), nchar(str2))
  
  while (overlap > 0) {
    if (substr(str1, nchar(str1) - overlap + 1, nchar(str1)) == substr(str2, 1, overlap)) { #check if end substr from string1 matches beginning substr of str2
      max_overlap <- overlap
      break
    }
    overlap <- overlap - 1 #if no match with the big overla, -1 and repeat
  }
  
  return(max_overlap)
}


# Function that runs "find_overlapping_substrings" iteratively until certain conditions are met
repeat_until_condition_met <- function(km, threshold) {
  while (TRUE) {
    #############################1st step -- build overlp matrix, rows are str1 and cols str 2
    overlap_matrix <- matrix(0, nrow = length(km), ncol = length(km))
    colnames(overlap_matrix) <- km
    rownames(overlap_matrix) <- km
    # calculate overlaps and store in matrix
    for (i in 1:length(km)) {
      for (j in 1:length(km)) {
        if (i != j) {
          overlap_matrix[i, j] <- find_overlapping_substrings(km[i], km[j])
        }
      }
    }
    # check if max val in overlap_matrix  < threshold
    if (max(overlap_matrix) < threshold) {
      break  
    }
    
    ########################### 2nd step
    over6 <- list()
    
    for (i in 1:nrow(overlap_matrix)) {
      for (j in 1:ncol(overlap_matrix)) {
        if (overlap_matrix[i, j] > (max(overlap_matrix) - 1)) {
          # append row (str1) as key and column as value to dict
          over6[[rownames(overlap_matrix)[i]]] <- colnames(overlap_matrix)[j]
        }
      }
    }
    
    ########################### 3rd step
    mer8 <- list()
    
    for (key in names(over6)) {
      value <- over6[[key]]
      ###################################
      if (length(key) == length(value)){
        last_letter <- substr(value, (length(value) - (length(value)-max(overlap_matrix)))+1,nchar(value)+1)
      } else if (length(key) > length(value)){
        last_letter <- substr(value, (length(value) - (length(value)-max(overlap_matrix)))+1,nchar(value)+1)
      } else {
        last_letter <- substr(value, (length(value) - (length(value)-max(overlap_matrix)))+1,nchar(value)+1)
      } 
      ################################### neeed to write the condition where str2 is contained in str1 but str1 is not elongated
      modified_key <- paste0(key, last_letter)  # append last part of str2 to str1
      mer8[[modified_key]] <- value
    }
    
    # 4th step
    keys_vector <- names(over6)
    values_vector <- unlist(over6)
    
    combined_vector <- c(keys_vector, values_vector)
    
    uni <- unique(combined_vector)
    km <- km[!km %in% uni]
    
    longer <- names(mer8)
    km <- c(km, longer)
  }
  return(km)
}


######################### 
# ANALYSIS
#########################
print("Starting analysis")

plot_list <- list()
long_and_7mers <- list()

print(paste("Number of 7mers: ", length(longmers)))

for (i in longmers){
  seven <- kmers 
  ## distance matrix to know which kmers have formed the longermer
  print("1")
  dist_matrix <- adist(i, seven, costs = c(insertion = 7, deletion = 7, substitution =4))
  similarity_matrix <- 1 / (1 + dist_matrix)# Convert the distance matrix into a similarity matrix
  rownames(similarity_matrix)<-i
  colnames(similarity_matrix)<-seven
  print("2")
  b <- colnames(similarity_matrix)[which(similarity_matrix == max(similarity_matrix))]
  ## Calculate correlation matrix, to discard kmers that dont have correlation with activities or add kmers that have a mismatch and correlate in activity
  cells<- rownames(ann)[which(ann$typecell == j)] 
  print("3")
  acts <- all[which(rownames(all) %in% b),cells]
  ac <- t(acts)
  print("4")
  ac <- as.data.frame(ac)
  print("5")
  row_means <- rowMeans(acts)

  ## Normalize the dataset by subtracting the mean of each row
  corr <- cor(ac, method = "spearman")
  cor_df <- as.data.frame(as.table(corr))
  random_indices <- sample(1:nrow(all), 200)
  random_hu <- all[random_indices, cells]
  random_hul <- t(random_hu)
  random_hul <- as.data.frame(random_hul)
  print("6")

  # Calculate correlation for random 100 to set the threshold
  ju_random <- cor(random_hul, method = "spearman")
  corr_values_random <- as.data.frame(as.table(ju_random))
  quantile_value <- quantile(corr_values_random$Freq, 0.8) #threshold

  ## set threshold to eliminate the kmers that dont correlate
  print("7")
  filtered_df <- filter_and_remove_low_occurrences(cor_df, threshold = quantile_value)
  if (!is.null(filtered_df) && nrow(filtered_df) >=9){ ## at least its composed with 3 kmers that have more correlation in their activities than 0.5
    # plot <- ggplot(filtered_df, aes(Var1, Var2, fill = Freq, label = sprintf("%.2f", Freq))) +
    # geom_tile(color = "white") +
    # geom_text(color = "black", size = 5) +
    ##scale_fill_gradient(low = "blue", high = "red", name = "Correlation") +
    # scale_fill_distiller(palette = "RdBu", direction = -1, name = "Correlation",limits = c(0, 1)) +
    # theme_minimal() +
    # labs(title = paste("Correlation 7mers of longmer", i ,"in", j), x = " ", y = " ")
    # plot_list[paste(i,j)] <- plot
    # print(plot)

    long_and_7mers[[i]] <- unique(c(filtered_df$Var1)) ## 7mers that are part of the longmer
    }}

# For each long and 7mers now build the longermers again, its going to be super fast
threshold <- 5


# df <- as.data.frame(long_and_7mers)
# Convert each element of the list to a character vector
long_and_7mers_char <- lapply(long_and_7mers, as.character)

# Determine the maximum length of the vectors
max_length <- max(sapply(long_and_7mers_char, length))

# Extend each vector to the maximum length by adding NAs
long_and_7mers_char <- lapply(long_and_7mers_char, function(x) {
  length(x) <- max_length
  return(x)
})

# Create a data frame from the list of character vectors
df <- do.call(data.frame, long_and_7mers_char)

# Set the column names of the data frame to the names of the list elements
colnames(df) <- names(long_and_7mers)
# print(df)
h <- long_and_7mers

# print(paste("Number of 7mers after filtering: ", length(h)))
# print("9")

# Print the number of 7mers after filtering
print(paste("Number of 7mers after filtering: ", length(long_and_7mers)))

for (i in 1:ncol(df)){
  df["old_longmer",i] <- colnames(df)[i]
  yeh <- c(as.character(h[[i]]))
  final_km <- repeat_until_condition_met(yeh, threshold)
  lengths <- sapply(final_km, nchar)
  # Find the index of the longest string (first occurrence in case of ties)
  longest_index <- which(lengths == max(lengths))[1]
  # Retrieve the longest string
  longest_string <- final_km[longest_index]
  colnames(df)[i] <- longest_string
}


saveRDS(df, subs_dfs_paths_down)
