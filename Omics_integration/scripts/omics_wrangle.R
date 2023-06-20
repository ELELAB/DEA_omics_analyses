### This script wrangles different omics data to a gene-level base
library(tidyverse)
library(purrr)
library(stringr)
library(readr)
library(data.table)
source("data_wrangle_functions.R")

# Load preprocessed data -------------------------------------------------------
# Expression
dataDEGs <- get(load("../data/BRCA_dataDEGs.rda"))
dataFilt_hugo <- get(load("../data/BRCA_dataFilt_HUGO.rda"))

# CNV
BRCA_CNV_gistic <- read.table(file = "../data/GISTIC_all_data_by_genes.txt",
                              sep="\t",
                              header=TRUE)

# Mutation
BRCA_mu <- read.csv("../data/mutations.csv")
mu_annotated <- get(load("../data/BRCA_mu_annotated.rda"))

# Methylation
probes <- read_tsv("../data/HM450.hg38.manifest.gencode.v36.tsv")
BRCA_met <- get(load("../data/BRCA_DNAMet.rda"))
rm(met_df)

# Common samples
common_samples <- read.csv("../data/common_samples.csv")$x




# Wrangle data -----------------------------------------------------------------
# Expression
dataFilt_hugo <- find_replicates(as.data.frame(dataFilt_hugo))
colnames(dataFilt_hugo) <- str_sub(colnames(dataFilt_hugo), 1, 15)

# CNV
CNV_common_samples <- CNV_wrangle(BRCA_CNV_gistic, common_samples)

# Methylation
BRCA_met <- find_replicates(BRCA_met)
Methy_common_samples <- Methylation_wrangle(methy_data = BRCA_met,
                                            probe_data = probes,
                                            samples = common_samples)

# Mutation
# Wrangle to mutation types based on consequence table 
mutation_list <- Mu_wrangle_type(mu_df = BRCA_mu,
                                 samples = common_samples,
                                 cscape = FALSE)

# Combine with CScape-somatic
mutation_list_cscape <- Mu_wrangle_type(mu_df = mu_annotated,
                                        samples = common_samples,
                                        cscape = TRUE)


# Save image -------------------------------------------------------------------
rm(list = c("BRCA_CNV_gistic", "BRCA_met", "probes", "mu_annotated",
            "CNV_wrangle","Methylation_wrangle", "Mu_wrangle_number", "Mu_wrangle_type", "find_replicates"))
save.image("../data/wrangled_data.RData")
