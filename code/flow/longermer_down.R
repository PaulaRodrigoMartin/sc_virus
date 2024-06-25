library(dplyr)
library(Biostrings)
library(stringr)
library(remotes)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]

# Activity scores (output from miReact)
cell <- readRDS(celltypes) 
# xm <- readRDS("C:/Users/paula/Documents/Aarhus master/msc_thesis/workflow/data/hs.seqXmot.counts.utr3_mrs_7mer.rds")
subs_dfs_paths_down <- gsub("annot4_", "longermer_downreg_", celltypes)

######################### 
# FUNCTIONS
#########################

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


# Function that aligns a motif to a sequence list and retrieves the matches with the specified amount of flanking bases
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

# Same function as before but accounting for motifs that align in the beginning or end of the sequence, creating shorter alignments
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


# Function to create lists for each column and finally create a regular expression
create_regexp <- function(column) {
  values_above_zero <- names(column[column > 0])
  
  if (length(values_above_zero) == 1) {
    return(list(values_above_zero))
  } else if (length(values_above_zero) > 1) {
    list_string <- paste("(", paste(values_above_zero, collapse = "|"), ")", sep = "")
    return(list(list_string))
  } else {
    return(list())
  }
}


#Function to create a regular expression but only allowing for G|U wobble
create_regexp_GUwobble <- function(column) {
  values_above_zero <- names(column[column > 0])
  
  if (length(values_above_zero) == 1) {
    return(list(values_above_zero))
  } else if (length(values_above_zero) > 1) {
    if ('G' %in% values_above_zero && 'T' %in% values_above_zero) {
      relevant_values <- values_above_zero[values_above_zero %in% c('G', 'T')]
      list_string <- paste("(", paste(relevant_values, collapse = "|"), ")", sep = "")
      return(list(list_string))
    } else {
      max_value_name <- names(column)[which.max(column)]
      return(list(max_value_name))
    }
  } else {
    return(list())
  }
}

# Find local maxima
find_peaks <- function(x) {
  peak_indices <- which(diff(sign(diff(x$y))) == -2) + 1
  data.frame(x = x$x[peak_indices], y = x$y[peak_indices])
}


######################### 
# ANALYSIS
#########################

## Construct longermers
threshold <- 5 ## min overlap length

## Divide dataset in UP regulated and DOWN regulated for analysis
cell_down <- cell[cell$diffexpressed_corr == "DOWN_corrected",]

cell_down_log2FC <- cell_down[order(cell_down$log2fc, decreasing = FALSE), ]
g <- rownames(cell_down_log2FC[1:100,])
cell_down_corrviralload <- cell_down[order(cell_down$corr_viralload, decreasing = TRUE), ]
h <- rownames(cell_down_corrviralload[1:100,])
cell_down$logttest <- -log10(cell_down$ttest_pval_corrected)
cell_down_ttest <- cell_down[order(cell_down$logttest, decreasing = TRUE), ]
y <- rownames(cell_down_ttest[1:100,])

unique_rownames <- unique(c(g,h,y))

cell_down_top <- cell_down[unique_rownames,]


print("downregulated longermers doing...")
## DOWN-REGULATED kmers
dict_down <- list()

final_km <- repeat_until_condition_met(rownames(cell_down_top), threshold) # cell[100] if want to check first 100 kmers, remember to first order them by relevance
longmers <- list()
for (i in final_km){
  if (nchar(i)>12){
    longmers <-append(longmers,i)}}
dict_down <- longmers
print("downregulated longermers done!")

saveRDS(dict_down, file = subs_dfs_paths_down)