# DEA_omics_analyses

## Title
Interpretation of gene expression in breast cancer through AI-guided analysis of multi-omics data

## Description
In this project, we aim to use genomics and epigenomics data from The Cancer Genome Atlas (TCGA) to predict and understand the dysregulation of genes’ transcriptional levels in breast cancer patients. Single omics data can only reveal a portion of the biological complexity of cancer, thus integration of multi-omics data is required to provide a comprehensive view of the molecular mechanism of carcinogenesis. To achieve this goal, we will first process and integrate each omics layer into a gene-level base, creating an input dataset with the target related to differential expression. We will then use QLattice, a supervised machine learning tool using symbolic regression, to find out explanations about the most important omics features and specific interactions between them that contribute to the prediction. Additionally, we will use other basic learning-based algorithms, such as LASSO regression and random forest for comparison of the performance and validation of the feature importance. To further enhance our understanding, we will employ multilayer perceptron (MLP), a deep learning method that can create more complex models with high performance for a complementary study, which allows us to benchmark and evaluate both the predictive power and interpretability of different approaches. Overall, by combining multi-omics data and using a range of machine learning algorithms, we aim to gain a deeper understanding of the molecular mechanisms driving dysregulation of gene expression, thus obtaining insights into cancer progression and potential targets for therapeutic intervention.

## Contents
**Omics_integration:**
Include folders related to implement the multi-omics integraion for generating the input data set. The data and resulting files are not included here due to the size.\

**SR_RF_LASSO:**
Include folders related to machine learning modeling on the multi-omics data by using symbolic regression, random forest and LASSO.\

**MLP:**
Include folders related to MLP modeling, consensus clustering and functional enrichment analysis.

## Requirements
The user needs to install and activiate an environment containing R packages and Python modules that need to be used. Here the environment we used on the server is under: 
`/data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/QLattice_env/`

## Main packages 
- R/Bioconductor: limma (3.40.6), TCGAbiolinks (2.12.5), cola (2.4.0), MoonlightR (1.24.0)
- Python: feyn (3.0.2), PyTorch (1.13.1+cu117), Scikit-learn (1.0.2)