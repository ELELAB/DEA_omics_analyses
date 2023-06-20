#-------------------- TCGA biolinks differential expression analysis -------------------------
library(TCGAbiolinks)
library(limma)

# Source revised DEA function
source("TCGAanalyze_DEA_revised.R")


# Import the filtered data file in "TCGA_data/Expression"
dataFilt <- get(load("../data/BRCA_dataFilt_HUGO.rda"))
  
# Divide dataFilt into tumor samples and normal samples
bar.codes <- colnames(dataFilt)
tum.samples <- TCGAquery_SampleTypes(bar.codes, "TP")
norm.samples <- TCGAquery_SampleTypes(bar.codes, "NT")
  
TP <- dataFilt[,which(colnames(dataFilt) %in% tum.samples)]
NT <- dataFilt[,which(colnames(dataFilt) %in% norm.samples)]
  
# Perform DEA
dataDEGs <- TCGAanalyze_DEA_revised(mat1 = NT,
                            mat2 = TP,
                            pipeline = "limma", 
                            batch.factors = "TSS",
                            Cond1type = "Normal",
                            Cond2type = "Tumor",
                            fdr.cut = 0.05,
                            logFC.cut = 1,
                            voom = TRUE)
save(dataDEGs, file= paste0("../data/BRCA_dataDEGs.rda"))
