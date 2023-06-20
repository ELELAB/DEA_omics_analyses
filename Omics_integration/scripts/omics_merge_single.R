### This script merges multi-omics data types into input for individual samples
library(stringr)
library(readr)
library(data.table)
library(tidyverse)
library(purrr)

# Load data --------------------------------------------------------------------
load("../data/wrangled_data.RData")



# Merge omics data for each sample ---------------------------------------------
omics_list <- list()
# Extract common genes in CNV and Methylation
common_genes <- intersect(CNV_common_samples$Gene.Symbol, Methy_common_samples$genesUniq)

# Assign genes with target variable
common_genes_Filt <- intersect(rownames(dataFilt_hugo), common_genes)

# Merge data
for (i in 1:length(mutation_list_cscape)){
  sample <- names(mutation_list_cscape)[i]
  if (sample %in% colnames(CNV_common_samples) & sample %in% str_sub(colnames(Methy_common_samples)[-1], 1, 15)){
    # Extract omics data of the sample
    mu_df <- mutation_list_cscape[[i]]
    CNV_df <- CNV_common_samples %>%
      select(Gene.Symbol, sample) %>%
      rename(cnv = sample)
    methy_df <- Methy_common_samples %>%
      select(genesUniq, contains(sample))
    str_sub(colnames(methy_df)[-1], 1, 15) = "Me"
    
    # Merge CNV and Methylation
    temp <- inner_join(x = CNV_df,
                       y = methy_df,
                       by = c("Gene.Symbol" = "genesUniq")) %>%
      filter(Gene.Symbol %in% common_genes_Filt)
    
    # Merge with mutation data
    merged_df <- left_join(x = temp,
                           y = mu_df,
                           by = c("Gene.Symbol" = "Hugo_Symbol")) # Only keep the overlapping genes among omics data types
    
    # Replace NA with 0
    merged_df[is.na(merged_df)] <- 0
    
    # Extract target from dataFilt
    Expr_df <- as.data.frame(dataFilt_hugo) %>%
      rownames_to_column(var = "Gene.Symbol") %>%
      select(Gene.Symbol, sample) %>%
      filter(Gene.Symbol %in% merged_df$Gene.Symbol)
    
    # Add target and chromosome information
    merged_df <- inner_join(merged_df, Expr_df, by = "Gene.Symbol") %>%
      rename(Counts = sample) %>%
      inner_join((BRCA_mu %>% 
                    distinct(Hugo_Symbol, Chromosome)),
                 by = c("Gene.Symbol" = "Hugo_Symbol")) %>% 
      column_to_rownames(var = "Gene.Symbol")
    
    
    # Fill in the list
    omics_list[[sample]] <- merged_df
  }
}

# Save data frames to separate files
mapply(
  write.table,
  x = omics_list,
  file = paste0("/data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/Omics_integration/results/Single_omics_files/", names(omics_list), ".csv"),
  MoreArgs=list(row.names=TRUE, sep=",")
)
