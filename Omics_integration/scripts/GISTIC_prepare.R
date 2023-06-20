### This script prepares the CNV segment file and marker file for running web-based GISTIC 2.0
library(tidyverse)
library(readr)
library(data.table)

# Load data --------------------------------------------------------------------
BRCA_CNV_seg <- get(load("../data/BRCA_cnv.exp.rda"))
rm(data)

# Prepare data for GISTIC analysis ---------------------------------------------
# Segment file
Seg_file <- BRCA_CNV_seg[,c(7,2,3,4,5,6)]
write.table(Seg_file, file = "../data/CNV_segment_file.txt", 
            sep = "\t", row.names = F,quote = F)

# Marker file
Marker <- fread("../data/snp6.na35.remap.hg38.subset.txt",
                data.table = F)

# Keep probesets with freqcnv = FALSE
Marker <- Marker[Marker$freqcnv == "FALSE",1:3]
colnames(Marker) <- c("Marker_Name", "Chromosome", "Marker_Position")
write.table(Marker, file = "../data/CNV_marker_file.txt", 
            sep = "\t", row.names = F, quote = F)
