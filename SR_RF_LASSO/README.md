## Introduction 
This project is to investigate the effect of omics data types on differential gene expression of breast cancer samples compared with normal samples. This folder is generated with the purpose to apply `QLattice` which uses symbolic regression to model the multi-omics data, and benchmark it with the other two machine learning methods (Random Forest and LASSO) in terms of performance and interpretability.

## Requirements
The user needs to activiate an environment containing R packages and Python modules that need to be used. The environment is activated using:
`conda activate /data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/QLattice_env/`

## Structure
**scripts:**
Include all the scripts/jupyter notebooks needed to model the multi-omics data\

**models:**
Include the resulting models of QLattice regarding different binary classifications\

**figures:**
* `omics_samples_merged.csv`: the input data with all the common samples. The genes include NonDEGs, up-regulated and down-regulated DEGs.
* `Merged_omics_files`: the per-sample-based data within each subtype folder\

**figures:**
Include important figures before and after data modeling

## Running the Scripts
The user can first change the directory to `/data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/` and use a romote connection to the jupyter notebook from your local PC following the process in the [link](https://amber-md.github.io/pytraj/latest/tutorials/remote_jupyter_notebook).

* `auxiliary.py`: Manually-defined functions used for preprocessing the multi-omics data and modeling.

* `aggregated_nonDE_QLattice.ipynb`: QLattice modeling for the aggregated input gene data that're classified into DEGs and NonDEGs.

* `aggregated_nonDE_benchmark.ipynb`: Bechmarking with the other ML algorithms on the aggregated dataset where genes are classified into NonDEGs(0) and DEGs(1).

* `aggregated_onlyDE_QLattice.ipynb`: QLattice modeling for the aggregated input DEG data that're classified into down-regulated DEGs(0) and up-regulated-DEGs(1).
