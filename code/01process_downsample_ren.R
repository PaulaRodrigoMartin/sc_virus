library(dplyr)
library(anndata)
library(matrixStats)

nik <- read_h5ad("/faststorage/project/sc_covid_PiB2023/data/sc_covid/COVID19_ALL.h5ad") 
infe <- read_h5ad("/faststorage/project/sc_covid_PiB2023/data/sc_covid/covid_19_virus.h5a")


## chekc for cells that are duplicated (remember underscore and upperscore, so change them)
# do it for the healthy cells
new_rownames <- sub("-", "_", rownames(nik$obs))
rownames(nik$obs) <- new_rownames ## if thereś duplicates it should retrieve an error

# now check if the cells in the infected dataset re also in the common one
#head(rownames(infe))
common_values <- intersect(rownames(nik$obs), rownames(infe$obs))
length(common_values) #1078 of the cells are duplicated, meaning they appear in the infected AND whole datasets, but the rest of the infected cells don appear in the whole dataset (????)
colnames(infe$obs)[colnames(infe$obs) == "samples"] <- "sampleID"


## First, eliminate the cells that have no corresponding infected cells. Rename the cell names to be the same in both datasets
##subset nik
niti <- subset(nik, nik$obs$majorType == c("Macro", "CD4", "CD8", "NK", "Neu", "Plasma", "Epi", "Mono", "B", "DC") | rownames(nik$obs) %in% common_values) ## add the second argument because there is wrong annotation either in infected or in the whole dataset
niti <- subset(niti, niti$obs$celltype != "Epi-AT2") # get rid of epithelials AT2 that are not in infected
## make new column "typecell" and collapse CD4 and CD8 to T
niti$obs$typecell <- niti$obs$majorType
niti$obs$typecell <- gsub("CD4", "T", niti$obs$majorType)
niti$obs$typecell <- gsub("CD8", "T", niti$obs$typecell)
## find indx that contain the epi-... cells
Idx_squ <- which(niti$obs$celltype == 'Epi-Squamous')  # majortype Epi
Idx_cil <- which(niti$obs$celltype == 'Epi-Ciliated')  # majortype Epi
Idx_sec <- which(niti$obs$celltype == 'Epi-Secretory')

## add in the typecell column the subsequent names for the cells
niti$obs$typecell[c(Idx_squ)] <- "Squamous"
niti$obs$typecell[c(Idx_cil)] <- "Ciliated"
niti$obs$typecell[c(Idx_sec)] <- "Secretory"

niti$obs$typecell <- as.factor(niti$obs$typecell)

## Equal the names of the cells, do it in the infected because yes
infe$obs$cellType <- gsub("Macrophage", "Macro", infe$obs$cellType)
infe$obs$cellType <- gsub("Neutrophil", "Neu", infe$obs$cellType)


## Check
common_values2 <- intersect(rownames(niti$obs), rownames(infe$obs)) ## dim(common_values2) --> 1078 supernice sanity check
idx_common <- which(rownames(infe$obs) %in% common_values) ## cellnames in infe$obs that are the same in the whole dataset

setequal(common_values, common_values2) ## TRUE so they contain exactly the same 
#### Properly check if the sampleIDs are the same

## Add information to whole dataset to the duplicated cells and eliminate them from the infe (the information has been translated to the full dataset)
# from now on i work with niti (whole dataset minus cells that are not represented in the infected dataset) and infe (infected dataset)
colnames(infe$obs) <- c("sampleid" ,"patients" ,"n_counts" ,"typecell")
names(infe$obsm) <- c("X_PCA" ,"X_TSNE") ## names() because its a list

##check
y <- length(unique(rownames(nik))) == nrow(nik)
unique(y) ##TRUE so all rows are unique in nik

z <- length(unique(rownames(infe))) == nrow(infe)
unique(z) ##TRUE

##eliminate the duplicated genes that they have anotated (colnames(infe))
x <- which(duplicated(colnames(infe))==T)
infe <- infe[,-x]

## eliminate the cells that are common in both datasets from infe and then concatenate
infe_notcommon <- infe[!(rownames(infe) %in% common_values),]

## add the cells that are not common
full <- concat(
  list(niti, infe_notcommon),
  axis = 0L,
  join = "outer"
)

## add the cells that are common with mapping
common_rownames <- intersect(rownames(infe$obs), rownames(niti$obs))
indices <- which(rownames(full) %in% common_rownames)
idx_n <- which(rownames(infe) %in% common_rownames)
n_cou <- infe$obs$n_counts[idx_n]


indices_full<- data.frame(indices, n_cou)
#full$obs$indices_full <- 1:nrow(full$obs)

# This one takes a lot of time
for (i in 1:nrow(indices_full)){
  zeta <- indices_full[i,"indices"]
  full$obs[zeta,"n_counts"] <- indices_full[i,"n_cou"]
  print(i)
}

full$obs$Is_infected <- F

full$obs <- full$obs %>%
  mutate(Is_infected = ifelse(is.na(n_counts), FALSE, TRUE))

####################################################################################################################
full <- read_h5ad("/faststorage/project/sc_covid_PiB2023/data/sc_covid/COVID19_VIRUS_HEA.h5ad")
fu  <- t(full)
seqs <- readRDS("/faststorage/project/jsp_student_projects/sc_covid_PiB2023/miReact/seqs/hs.utr3.seqs.rds")
idx <- match(rownames(fu),seqs$gsym)
# remove genes without 3'UTRs from expression matrix
fu <- fu[!is.na(idx),]
# rename genes to ensemble geneIDs
rownames(fu) <- seqs$gid[idx[!is.na(idx)]]

cou <- sum(grepl("ENS", rownames(fu)))
# Percentage of genes with ensemble ID
cou/length(colnames(full))

tags <- c("Mono"   ,   "B"    ,     "T"   ,      "NK"   ,     "Neu"  ,     "DC"   ,     "Macro"  ,   "Plasma"  ,  "Ciliated"  ,"Secretory", "Squamous" )
fu <- t(fu)
summary(fu$obs$typecell) ## Itś typecell, not celltype

## Subset the dataset for the patients for which we have infected cells
unique(fu$obs$patients) %in% unique(fu$obs$PatientID)

# Replace NA values in fu$obs$PatientID with corresponding values from fu$obs$patients
fu$obs$PatientID[is.na(fu$obs$PatientID)] <- fu$obs$patients[is.na(fu$obs$PatientID)]

subset_fu <- fu[fu$obs$PatientID %in% unique(fu$obs$patients), ]
fu<-subset_fu


ratios <- numeric()
for (celltype in tags) {
  
  # Calculate the ratio for the current celltype
  ratio <- (sum(fu$obs$Is_infected[fu$obs$typecell == celltype] == TRUE)) / 
    (sum(fu$obs$typecell == celltype))
  
  # Store the ratio in the 'ratios' vector
  ratios[celltype] <- ratio
}
ratios

# Add viral load
viralload <- infe$X[,23559] 
cellIDs <- names(viralload)

# Initialize a viral load column with NA
fu$obs$viral_load <- NA

# Match cell IDs and assign viral load values
matching_indices <- match(rownames(fu$obs), cellIDs)
matched_values <- viralload[na.omit(matching_indices)]

# Assign matched values to corresponding rows in the viral load column
fu$obs$viral_load[!is.na(matching_indices)] <- matched_values
fu$obs$viral_load[!fu$obs$Is_infected] <- NA

vi <- t(as.matrix(fu$X))

## DO ANNOTATION FILES, it would be nice to add more annotation

mi <- data.frame(
  row.names = rownames(fu),
  Is_infected = fu$obs$Is_infected,
  typecell = fu$obs$typecell,
  viral_load = fu$obs$viral_load,
  patientID = fu$obs$PatientID,
  sub_celltype = fu$obs$celltype
)

## SAVE DOWNSAMPLED FILES
write.csv(mi, file = "/faststorage/project/jsp_student_projects/sc_covid_PiB2023/msc/scr/annotation_patients_more.csv") #annotation table
write_h5ad(anndata = fu, filename = "/faststorage/project/jsp_student_projects/sc_covid_PiB2023/msc/scr/th_mireact_downsampled_patients.h5ad") #subset dataset (11 celltypes, 3085 infected cells)
# X <- t(as.matrix(fu_downsampled$X))
saveRDS(vi, file = "/faststorage/project/sc_covid_PiB2023/2/data2/th_mireact_exp_downsampled.rds") #expression matrix of subset dataset
