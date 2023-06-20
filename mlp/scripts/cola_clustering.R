### This script implements clustering for truly predicted DEGs from 
### MLP binary classification by using cola package.

library(cola)
library(tidyverse)

# Load data --------------------------------------------------------------------
DEGs <- read.csv("../results/biclass_trueDEGs_data.csv")
BRCA_DEGs <- get(load("../../Omics_integration/data/BRCA_dataDEGs.rda"))
rm(dataDEGs)


# Preprocessing ----------------------------------------------------------------
# Add labels to DEGs
BRCA_DEGs <- BRCA_DEGs %>% 
  mutate(down_up = case_when(logFC < 0 ~ "down",
                              logFC > 0 ~ "up")) %>% 
  rownames_to_column("Gene.Symbol") %>% 
  right_join(DEGs, by = "Gene.Symbol") %>% 
  select(Gene.Symbol, down_up)

# Transpose the matix
DEG_mat <- t(as.matrix(DEGs %>% 
                         column_to_rownames(var = "Gene.Symbol") %>% 
                         select(-c(gene_type))))

# Adjust matrix???
DEG_mat <- adjust_matrix(DEG_mat)

# Clustering with single partition method --------------------------------------
kmeans_res = consensus_partition(DEG_mat,
                          top_value_method = "ATC",
                          top_n = c(1, 2, 5, 8, 9),
                          partition_method = "skmeans",
                          max_k = 6,
                          p_sampling = 1,
                          partition_repeat = 50,
                          anno = data.frame(direction = BRCA_DEGs$down_up))


# Visualization ----------------------------------------------------------------
# Statistics along with different k
select_partition_number(kmeans_res)

# The heatmap for the consensus matrix with a certain k
consensus_heatmap(kmeans_res, k = 2)

# Membership heatmap
membership_heatmap(kmeans_res, k = 2)

# Dimension reduction plot
dimension_reduction(kmeans_res, k = 2)

# Get signature
get_signatures(kmeans_res, k = 2)


# Retrieve the genes in clusters -----------------------------------------------
cluster1 <- get_classes(kmeans_res, k = 2) %>% 
  rownames_to_column("Gene.Symbol") %>% 
  filter(class == 1) %>% 
  pull(Gene.Symbol)

cluster2 <- get_classes(kmeans_res, k = 2) %>% 
  rownames_to_column("Gene.Symbol") %>% 
  filter(class == 2) %>% 
  pull(Gene.Symbol)

# Save the genes
write.csv(cluster1,
          file="../results/cola_cluster1.csv",
          row.names = F)

write.csv(cluster2,
          file="../results/cola_cluster2.csv",
          row.names = F)
