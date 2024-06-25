library(tidyverse)
# library(ggplot2)
# library(Rtsne)
# library(RColorBrewer)
library(anndata)

annot2_cil <- read_h5ad("/home/paulilokiestudia/jsp_student_projects/sc_covid_PiB2023/msc/log/annot2_Ciliated.h5ad")
annot2_sec <- read_h5ad("/home/paulilokiestudia/jsp_student_projects/sc_covid_PiB2023/msc/log/annot2_Secretory.h5ad")
annot2_sq <- read_h5ad("/home/paulilokiestudia/jsp_student_projects/sc_covid_PiB2023/msc/log/annot2_Squamous.h5ad")
mir <- readRDS("/home/paulilokiestudia/jsp_student_projects/sc_covid_PiB2023/msc/data/humir.rds")

files <- c(annot2_cil, annot2_sec, annot2_sq)

for (i in seq_along(files)){
  h <- files[i]
  h <- h[[1]]$obs
  h$diffexpressed_corr <- "NO"
  h$diffexpressed_corr[h$log2fc > 0.6 & h$ttest_pval_corrected < 0.05] <- "UP_corrected"
  h$diffexpressed_corr[h$log2fc < -0.6 & h$ttest_pval_corrected < 0.05] <- "DOWN_corrected"
  c <- h[h$diffexpressed_corr %in% c ("UP_corrected"),]
  d <- length(rownames(c))
  e <- mir[which(mir$seed7target %in% rownames(c)),]
  f <- dim(e)[1]
  cat("------------------\n")
  cat(i)
  cat("\nUPREG\n")
  cat(d)
  cat("\n-mir--\n")
  cat(f)
  c <- h[h$diffexpressed_corr %in% c ("DOWN_corrected"),]
  d <- length(rownames(c))
  e <- mir[which(mir$seed7target %in% rownames(c)),]
  f <- dim(e)[1]
  cat("------------------\n")
  cat(i)
  cat("\nDOWNREG\n")
  cat(d)
  cat("\n-mir--\n")
  cat(f)
  c <- h[h$diffexpressed_corr %in% c ("UP_corrected", "DOWN_corrected"),]
  d <- length(rownames(c))
  e <- mir[which(mir$seed7target %in% rownames(c)),]
  f <- dim(e)[1]
  cat("------------------\n")
  cat(i)
  cat("\nTOTAL\n")
  cat(d)
  cat("\n-mir--\n")
  cat(f)
}