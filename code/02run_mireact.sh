#!/bin/bash
#SBATCH --partition normal
#SBATCH --account sc_covid_PiB2023
#SBATCH --mem 60g
#SBATCH -t 6:00:00
Rscript /faststorage/project/jsp_student_projects/sc_covid_PiB2023/msc/scr/02run_mireact.R
