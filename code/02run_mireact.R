#!/usr/bin/env Renv
## USAGE: mireact_downsampled.R a output_file

library(dplyr)
source("/faststorage/project/jsp_student_projects/sc_covid_PiB2023/miReact/code/mireact.R")
efile <- "/faststorage/project/jsp_student_projects/sc_covid_PiB2023/msc/datax/th_mireact_downsampled_patients.rds"
annot <- read.csv("/faststorage/project/jsp_student_projects/sc_covid_PiB2023/msc/data/annotation_patients_more.csv"
, header = TRUE, row.names = 1)
# This following command starts miReact.
mireact(exp=efile, species="hs", motifs=7, seq.type="utr3",
        out.file="/faststorage/project/jsp_student_projects/sc_covid_PiB2023/msc/res/mirnaActivity_th.rds",out.meonly=T,
        mail="paulilokiestudia@gmail.com", install.dir = "/faststorage/project/jsp_student_projects/sc_covid_PiB2023/miReact", 
        annotation = annot)