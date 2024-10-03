# National-Workshop-Microbiome-Informatics

## Introduction

This workshop will help the newbies to learn how to analyse the 16s rRNA amplicon data using R.
![IMG-20240911-WA0009](https://github.com/user-attachments/assets/237bb872-dcb2-4ec8-8304-7c81798bb017)

## Pre-requisities
R should be installed in the user system/PC. R installation steps are given in https://cran.r-project.org/bin/windows/base/R-4.4.1-win.exe. 

RStudio should also be installed, https://download1.rstudio.org/electron/windows/RStudio-2024.04.2-764.exe

The following R packages are required for the analysis the 16s rRNA amplicon data:
## R packages

  if (!requireNamespace("BiocManager", quietly = TRUE))
                                install.packages("BiocManager")
  
  BiocManager::install("dada2")
  
  BiocManager::install("phyloseq")
  
  install.packages("microeco")
  
  install.packages("microeco")
  
  install.packages("ggplot2")
  
  install.packages("ggpubr")
  
  install.packages("ggdendro")

## Workflow
![workflow_v2](https://github.com/user-attachments/assets/751d7781-eddb-4539-869f-cd26a90d7c56)



## Raw FASTQ files and reference taxonomy file (Silva_138)
Download amplicon sequencing raw data (.zip) files for analyses, following this link:
https://github.com/RoyDibakar/National-Workshop-Microbiome-Informatics/archive/refs/heads/main.zip

Download the reference taxonomy file using the following link :
https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz

### Unzip the folder; keep the reference taxonomy file and the FASTQ files in the same folder

The R codes can be downloded from the following link https://github.com/RoyDibakar/National-Workshop-Microbiome-Informatics/blob/main/Codes_Microbiome_Informatics.R

## Team
Dibakar Roy, Paramita Roy, and Sudipto Saha
