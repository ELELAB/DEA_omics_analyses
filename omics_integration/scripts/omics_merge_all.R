### This script merges multi-omics data types into input including all samples
library(stringr)
library(readr)
library(data.table)
library(tidyverse)
library(purrr)

# Load data --------------------------------------------------------------------
load("../data/wrangled_data.RData")



# Merge omics data of all samples ----------------------------------------------
# Find different samples
diff_samples <- setdiff(str_sub(colnames(Methy_common_samples)[-1], 1, 15), 
                        colnames(CNV_common_samples)[-1])
colnames(CNV_common_samples)[-1] <- paste0("Cn_", colnames(CNV_common_samples[,-1]))
colnames(Methy_common_samples)[-1] <- paste0("Me_", colnames(Methy_common_samples[,-1]))

# Delete different samples
Methy_common_samples <- Methy_common_samples[,-which(str_sub(colnames(Methy_common_samples), 1, 18) == paste0("Me_", diff_samples))]
mutation_list_cscape <- mutation_list_cscape[-which(names(mutation_list_cscape) == diff_samples)]

# Change column names in mutation list of data frames
for (i in 1:length(names(mutation_list_cscape))){
  colnames(mutation_list_cscape[[i]])[-1] <- paste0(colnames(mutation_list_cscape[[i]][,-1]),
                                                   "_",
                                                   names(mutation_list_cscape)[i])
}

# Pivot mutation data wider
mutation <- mutation_list_cscape %>% 
  reduce(full_join, by = "Hugo_Symbol") %>% 
  replace(is.na(.), 0)  

# Extract target from dataDEG
DEGs <- dataDEGs %>% 
  rownames_to_column(., var = "Gene.Symbol") %>% 
  select(Gene.Symbol, logFC)

# Extract target from expression data
no_genes <- data.frame(setdiff(rownames(dataFilt_hugo), rownames(dataDEGs)), 
                       rep('Non_sign', length(setdiff(rownames(dataFilt_hugo), rownames(dataDEGs)))))
colnames(no_genes) <- colnames(DEGs)
Genes <- rbind(DEGs, no_genes)  


# Merge 3 data sets with non-DE genes
omics_merged_nonDE <- inner_join(CNV_common_samples, Methy_common_samples, by = c("Gene.Symbol" = "genesUniq")) %>% 
  inner_join(mutation, by = c("Gene.Symbol" = "Hugo_Symbol")) %>% 
  inner_join(Genes, by = "Gene.Symbol") %>% 
  inner_join((BRCA_mu %>% 
                distinct(Hugo_Symbol, Chromosome)),
             by = c("Gene.Symbol" = "Hugo_Symbol"))


write.csv(omics_merged_nonDE,
          file = "../results/multiomics_all_samples.csv",
          row.names = FALSE)
