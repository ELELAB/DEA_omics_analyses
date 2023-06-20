# DEA_omics_analyses

## Project Title
Interpretation of gene expression in breast cancer through AI-guided analysis of multi-omics data

## Description
Our project aims to use genomics and epigenomics data from The Cancer Genome Atlas (TCGA) to predict and understand gene dysregulation in breast cancer patients. Integrating multi-omics data is necessary to fully comprehend the complexity of cancer biology. We processed and integrated each omics layer into a gene-level dataset and used machine learning algorithms such as symbolic regression, LASSO regression, random forest, and multilayer perceptron (MLP) to identify important features and interactions contributing to the prediction. By combining multi-omics data and using various machine learning approaches, we aim to gain insights into the molecular mechanisms behind gene dysregulation and breast cancer progression, paving the way for development of potential therapeutic targets.

## Contents
**Omics_integration:**
Include folders related to implement the multi-omics integraion for generating the input data set. The data and resulting files are not included here due to the size.

**SR_RF_LASSO:**
Include folders related to machine learning modeling on the multi-omics data by using symbolic regression, random forest and LASSO.

**MLP:**
Include folders related to MLP modeling, consensus clustering and functional enrichment analysis.

## Requirements
The user needs to install and activiate an environment containing R packages and Python modules that need to be used. Here the environment we used on the server is under: 
`/data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/QLattice_env/`

## Main packages 
- R/Bioconductor: limma (3.40.6), TCGAbiolinks (2.12.5), cola (2.4.0), MoonlightR (1.24.0)
- Python: feyn (3.0.2), PyTorch (1.13.1+cu117), Scikit-learn (1.0.2)