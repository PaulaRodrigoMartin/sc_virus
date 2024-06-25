library(dplyr)
library(anndata)
library(Biostrings)
library(stringr)
args<-commandArgs(trailingOnly = TRUE) 
###############This is what i need to modify
activity_mireact <- args[1]
# print(activity_mireact)
annotations  <- args[2]
mir_annotations <- args[3]
rbp_annotations <- args[4]
xmot_annot <- args[5]
fasta_file_path <- args[6]


# Activity scores (output from miReact)
a <- readRDS(activity_mireact) 

# Annotation table fed to mireact
annot <- read.csv(annotations, header = TRUE, row.names = 1)

# Annotated miRNA data
mir<- readRDS(mir_annotations)

# Annotated rbp data
rbp <- readRDS(rbp_annotations)

# Motif occurrence across genes
xmot <- readRDS(xmot_annot)

subs_dfs_paths <- c()
for (i in unique(annot$typecell)){
  subs_dfs_paths <- append(subs_dfs_paths, paste0( "log/df_", i, ".h5ad"))
}

##Transform miReact output in an anndata object
# Check if cellIDs are in the same order
if (identical(colnames(a), rownames(annot))){
  # Create anndata object from annotations
  ad <- AnnData(
    X = a,
    obs = data.frame(mirna_target = rep(NA, 16384), mirna_seed =rep(NA, 16384), rbp =rep(NA, 16384), genes = rep(NA, 16384), GO_terms = rep(NA, 16384),  row.names = rownames(a)),
    var = data.frame(annot)
  )
  
  ## Add mirna annotations (seed and target sites) to obs (kmers)
  # Target site
  print("mirna target or seed site annotation...\n")
  replace_value <- match(rownames(a),mir$seed7target)
  names <- mir[replace_value,]$name
  ad$obs$mirna_target <- names
  ad$obs$is_mirnatarget <- FALSE
  ad$obs$is_mirnatarget[replace_value] <- TRUE
  
  # Seed site
  replace_value <- match(rownames(a),mir$seed7)
  names <- mir[replace_value,]$name
  ad$obs$mirna_seed <- names
  ad$obs$is_mirnaseed <- FALSE
  ad$obs$is_mirnaseed[replace_value] <- TRUE
  print("mirna target or seed site annotation done!\n")

  print("RBP annotation...\n")
  # RBP annotations
  ad$obs$rbp <- FALSE
  for (i in 1:nrow(rbp)){
    align <- c()
    rbp_name <- rbp[i,"RBP"]
    for (j in seq_along(rbp$kmers[[i]])){
      result <- grep(rbp$kmers[[i]][[j]], rownames(a))
      align <- append(align, result)
    }
    ad$obs[align, "rbp"] <- rbp_name
  }
  print("RBP annotation done!\n")
  # n' hits in 3'UTRs
  print("n' hits in 3'UTRs ...\n")
  ad$obs$n_hits_3UTR <- rowSums(xmot)
  print("n' hits in 3'UTRs done!\n")

  # n' hits in SARS-CoV-2 genome
  # calculate hits in sars-cov-2 virus
  # print("n' hits in SARS-CoV-2...\n")

  # sequences <- readDNAStringSet(fasta_file_path)
  # print("SARS-CoV-2 read...\n")
  # ad$obs$sarshit <- NA_character_
  # for (i in 1:nrow(ad$obs)){
  #   ad$obs$sarshit[i] <- vcountPattern(rownames(ad$obs)[i], sequences)
  # }
  # print("n' hits in SARS-CoV-2 done!\n")

  # Subset anndata per cell type and store each subsetted dataframe
  subs_dfs <- c()
  for (i in unique(annot$typecell)){
    dfcell <- ad[,ad$var$typecell == i] # subset dataframe
    assign(paste0("df_", i), dfcell) # Assign the subsetted dataframe to a variable with a dynamic name
    subs_dfs <- append(subs_dfs, paste0("df_", i))
  }
  
  # Save each dataframe
  for (i in 1:11) {
    write_h5ad(get(subs_dfs[i]), file = subs_dfs_paths[i])
  }
  # Also create some figures -->
    # Distribution of hits across 3'UTRs
    # Distribution of hits across SARS-CoV-2 virus
} else {
  print("Order of cellIDs in activity matrix and annotation is not the same")
}
