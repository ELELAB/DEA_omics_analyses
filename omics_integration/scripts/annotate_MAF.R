### This script annotates MAF file with CScape-somatic scores
library(tidyverse)
library(purrr)
library(stringr)
library(readr)
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(seqminer)
source("cscape_functions.R")


# Load original MAF file -------------------------------------------------------
BRCA_mu <- read.csv("../data/mutations.csv")



# Data wrangling ---------------------------------------------------------------
# Lift MAF to an old genome old needed by CScape
mu_hg19 <- LiftMAF(Infile = BRCA_mu, Current_Build = 'GRCh38')

# Prepare input data for CScape
cscape_in <- MAFtoCscape(mu_hg19)

# Run CScape
cscape_somatic_output <- RunCscape_somatic(input = cscape_in,
                                           coding_file = "../data/css_coding.vcf.gz",
                                           noncoding_file = "../data/css_noncoding.vcf.gz")


# Add cscape_columns in case one type is not found
cscape_cols <- c(Coding_score = NA, Noncoding_score = NA, Remark = NA)
cscape_somatic_output <- cscape_somatic_output %>%
  add_column(., !!!cscape_cols[setdiff(names(cscape_cols), names(.))]) %>%
  rename_with(.cols = c("Coding_score","Noncoding_score","Remark"),
              .fn = ~ paste("CScape_", ., sep = "")) %>%
  mutate(Variant_Type = "SNP")

# merge cscape results
mu_annotated_19 <- mu_hg19 %>%
  separate(Chromosome, into = c(NA, "Chr"), sep = 3, remove = FALSE, convert = TRUE)%>%
  mutate(Mutant = case_when(Reference_Allele == Tumor_Seq_Allele1 ~ Tumor_Seq_Allele2,
                            Reference_Allele == Tumor_Seq_Allele2 ~ Tumor_Seq_Allele1)) %>%
  left_join(cscape_somatic_output,
            by = c("Start_Position" = "Position",
                   "Variant_Type",
                   "Chr",
                   "Mutant")) %>%
  mutate(CScape_Mut_Class = case_when((CScape_Coding_score > 0.5 | CScape_Noncoding_score > 0.5) ~ "Driver",    #Driver
                                      (CScape_Coding_score <= 0.5 | CScape_Noncoding_score <= 0.5 ~ "Passenger"), #Passenger
                                      TRUE ~ "Unclassified")) %>%  #When no score is found
  unique()

#Lift back to 38
mu_annotated <- LiftMAF(Infile = mu_annotated_19,
                        Current_Build = "GRCh37")

# Save annotated mutation data
save(data = mu_annotated, file = "../data/BRCA_mu_annotated.rda")
