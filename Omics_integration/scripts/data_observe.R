### This script visualizes the raw data and finds common samples and 
### samples having replicates in omics data types

library(tidyverse)
library(stringr)
library(UpSetR)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(TCGAbiolinks)


# Download omics data ----------------------------------------------------------
# Clinical data
dataClin <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

# Expression
BRCA_DEGs <- get(load("../data/BRCA_dataDEGs.rda"))
rm(dataDEGs)
BRCA_dataFilt_hugo <- get(load("../data/BRCA_dataFilt_HUGO.rda"))
rm(dataFilt_hugo)

# CNV
BRCA_CNV <- get(load("../data/BRCA_cnv.exp.rda"))
rm(data)

# Mutation
BRCA_mu <- read.csv("../data/mutations.csv")

# Methylation
BRCA_methy <- get(load("../data/BRCA_DNAMet.rda"))
rm(met_df)
# Change column names into right format
colnames(BRCA_methy) <- gsub("\\.", "-", colnames(BRCA_methy))



# Observe data -----------------------------------------------------------------
# 1) Target distribution
temp_DEG <- BRCA_DEGs %>% 
  mutate(FC = 2^(logFC),
         Type = case_when(logFC >= 1 ~ "Up-regulated",
                          logFC <= -1 ~ "Down-regulated"))
# logFC distribution
p1 <- ggplot(temp_DEG, aes(x=logFC)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666") 

# logFC distribution between groups
p2 <- ggplot(temp_DEG, aes(x=logFC, fill=Type)) +
  geom_histogram(binwidth=.5, alpha=.5, position="identity")

p1 + p2

# 2) Find overlapping samples in all data types
# Build up lists
sample_list <- list("Expression" = unique(str_sub(colnames(BRCA_dataFilt_hugo), 1, 15)),
                     "CNV" = unique(str_sub(BRCA_CNV$Sample, 1, 15)),
                     "Mutation" = unique(str_sub(BRCA_mu$Tumor_Sample_Barcode, 1, 15)),
                     "Methylation" = unique(str_sub(colnames(BRCA_methy), 1, 15)))



# UpSet Plot
pdf(file = "/data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/Omics_integration/figures/sample_upset_plot.pdf",
    width = 12, height = 8)

# Input matrix of UpSet function
sample_matrix <- make_comb_mat(sample_list)

sample_upset <- UpSet(sample_matrix,
                   set_order = c("Expression", "CNV", "Mutation", "Methylation"),
                   pt_size = unit(5, "mm"),
                   lwd = 3,
                   height = unit(4, "cm"),
                   comb_col = c("deepskyblue3", "darkorange", "firebrick", "darkolivegreen")[comb_degree(sample_matrix)],
                   top_annotation = upset_top_annotation(sample_matrix,
                                                         height = unit(12, "cm"),
                                                         bar_width = 0.7,
                                                         axis_param = list(side = "left", at = seq(0,max(comb_size(sample_matrix))+100,100)),
                                                         annotation_name_side = "left",
                                                         annotation_name_gp = gpar(cex = 1),
                                                         annotation_name_offset = unit(1.5,"cm")),
                   right_annotation = upset_right_annotation(sample_matrix,
                                                             width = unit(3, "cm"),
                                                             gp = gpar(fill = "darkseagreen"),
                                                             axis_param = list(at = seq(0,max(set_size(sample_matrix))+100,200)),
                                                             annotation_name_offset = unit(1.5, "cm")),
                   
                   row_names_gp = gpar(fontsize = 18))

# Draw UpSet plot
sample_upset <- draw(sample_upset)

# Add intersection set sizes to each bar
decorate_annotation("intersection_size", {
  grid.text(comb_size(sample_matrix)[column_order(sample_upset)],
            x = 1:(ncol(sample_matrix)),
            y = unit(comb_size(sample_matrix)[column_order(sample_upset)], "native") + unit(1.5, "mm"),
            gp = gpar(fontsize = 12, fontface = "bold"),
            just = "bottom",
            default.units = "native")
})

# Add title
grid.text(label = "Upset plot showing samples of four omics data types",
          x = unit(17, "cm"), y = unit(18, "cm"),
          gp = gpar(fontsize = 20))

#End with dev.off()
dev.off()


# Extract overlapping patients
common_samples <- extract_comb(sample_matrix, "1111")
# (All of them are tumor samples!!)

write.csv(common_samples,
          file="../data/common_samples.csv",
          row.names = F)


# 3) Mutation data visualzation
sample_mu_df <- BRCA_mu %>% 
  # Add sample code information
  mutate(Sample = str_sub(Tumor_Sample_Barcode, 1, 15)) %>% 
  # Keep input samples
  filter(Sample %in% common_samples)

transcript_effect <- c("5'Flank", "5'UTR", "Nonsense_Mutation",
                       "Splice_Site", "Splice_Region", "Intron",
                       "3'UTR", "3'Flank", "RNA")
mu_list <- list()
for (effect in transcript_effect){
  mu_list[[effect]] <- nrow(sample_mu_df[which(sample_mu_df$Variant_Classification == effect),])
}
