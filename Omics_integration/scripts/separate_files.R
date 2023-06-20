### This script separate the smaple files based on subtype information
library(tidyverse)
library(stringr)
library(readr)
library(data.table)
library(TCGAbiolinks)

# Get subtype information
subtypes <- PanCancerAtlas_subtypes() %>% 
  filter(cancer.type == "BRCA") %>% 
  select(pan.samplesID, Subtype_mRNA)


# Create folder
Dir <- "../results/Single_omics_files"
subDir <- c("Normal", "Basal", "Her2", "LumA", "LumB")
setwd(file.path(Dir))
for (sub in subDir) {
  if (!file.exists(sub)){
    dir.create(file.path(Dir, sub))
  } 
}


# Move files into corresponding folder
files <- list.files(pattern = "csv")
for (name in files){
  subtype <- subtypes[which(str_sub(subtypes$pan.samplesID, 1, 15) == str_sub(name, 1, 15)),]$Subtype_mRNA
  file.rename(from = file.path(name),
              to = file.path(subtype, name))
}
