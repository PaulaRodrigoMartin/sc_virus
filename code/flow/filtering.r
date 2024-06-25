library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(anndata)
library(gprofiler2)

args<-commandArgs(trailingOnly = TRUE) 

celltypes <- args[1]

# Activity scores (output from miReact)
b <- read_h5ad(celltypes) 
subs_dfs_paths <- gsub("annot2_", "annot3_", celltypes)
subs_dfs_pathss <- gsub("h5ad", "rds", subs_dfs_paths)

h <- b$obs
# Add a column to specify if UP- or DOWN- regulated (log2fc respectively positive or negative)
h$diffexpressed <- "NO"
h$diffexpressed_corr <- "NO"
h$diffexpressed_corr[h$log2fc > 0.6 & h$ttest_pval_corrected < 0.05] <- "UP_corrected"
h$diffexpressed_corr[h$log2fc < -0.6 & h$ttest_pval_corrected < 0.05] <- "DOWN_corrected"
h$delabel <- ifelse(rownames(h) %in% rownames(h)[head(order(h$ttest_pval_corrected), 10)], rownames(h), NA)

selected_rows <- h[h$diffexpressed_corr %in% c("DOWN_corrected", "UP_corrected"), ]
print(paste0("Significant down and upregulated motifs: ",dim(selected_rows)[1]))
print(table(selected_rows$diffexpressed_corr))

cat("Significant down and upregulated motifs: ")
cat(dim(selected_rows)[1])

selected_rows$GO_terms <- NA_character_
print(" GO terms annotation --> doing")
for (i in 1:nrow(selected_rows)){
  genes <- strsplit(selected_rows$target_genes[i], ", ")#[[1]]
  go_terms <- gost(query = genes, organism = "hsapiens")
# gostres <- paste(na.omit(go_terms$result$term_name), collapse = "|")
# Extract GO terms
  go_terms_list <- go_terms$result$term_name
  go_terms_list <- paste(na.omit(go_terms_list), collapse = "|")
  selected_rows$GO_terms[i] <- go_terms_list
}
# h$GO_terms <- NA_character_
# h$GO_terms <- apply(h$target_genes, 1, function(row) {
#   if (h$WRS_pval_corrected[row] < 0.01){
#     genes <- strsplit(h$target_genes[row], ", ")[[1]]
#     gostres <- gost(query = genes, organism = "hsapiens")
#     gostres <- list(gostres$result$term_name)
#   } else {
#     gostres <- NA_character_
#   }
#   return(gostres)
# })
print(" GO terms annotation --> finished")
print(head(unique(selected_rows$GO_terms)))

saveRDS(selected_rows, file = subs_dfs_pathss)

############################## filtering by more variables ######################################

# mirna_target <- selected_rows[!is.na(selected_rows$mirna_target), ]

# cat(subs_dfs_pathss)
# cat("How many of them are mirnas?")
# cat(dim(mirna_target))

# #mirna_target
# q_05 <- quantile(mirna_target$corr_viralload, 0.05)
# q_95 <- quantile(mirna_target$corr_viralload, 0.95)

# # Subset mirna_seed$cluster based on the quantiles
# signif_corr_viralload_target <- mirna_target[mirna_target$corr_viralload <= q_05 | mirna_target$corr_viralload >= q_95,]

# if (nrow(signif_corr_viralload_target) == 0) {
#   q_05 <- quantile(mirna_target$corr_viralload, 0.1)
#   q_95 <- quantile(mirna_target$corr_viralload, 0.9)
#   signif_corr_viralload_target <- mirna_target[mirna_target$corr_viralload <= q_05 | mirna_target$corr_viralload >= q_95,]
#     }

# print("Motifs after mirna + corr viralload filtering")
# print(table(signif_corr_viralload_target$diffexpressed_corr))

# saveRDS(mirna_target, file = subs_dfs_pathss)
# # saveRDS(signif_corr_viralload_target, file = subs_dfs_pathss)

# subs_dfs_pathsss <- gsub("rds", "txt", subs_dfs_pathss)
# subs_dfs_pathssss <- gsub("filtered", "numbers", subs_dfs_pathsss)

# file_conn <- file(subs_dfs_pathssss, "w")

# # Write the information to the file
# cat("Significant down and upregulated motifs: \n", file = file_conn)
# cat(dim(selected_rows)[1], "\n" ,file = file_conn)
# cat("out of 16384 motifs \n", file = file_conn)
# cat("\n")
# cat("Of which ___ are mirna targets \n", file = file_conn)
# cat(dim(mirna_target)[1], "\n" ,file = file_conn)
# cat("deregulation percentage \n", file =file_conn)
# cat(dim(selected_rows)[1]/16384*100, "\n", file =file_conn)
# cat("mirna percentage \n", file =file_conn)
# cat(dim(mirna_target)[1]/dim(selected_rows)[1]*100, file =file_conn)

# # Close the file connection
# close(file_conn)


# ## save selected_rows to run longermer pipeline
# ## with selected_rows, do GO terms

# ## save mirna_target + mirna_seed for further filtering