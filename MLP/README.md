## Introduction 
This folder is generated with the purpose to apply Multilayer Perceptron (MLP) to model the multi-omics data, and implement biological investigations by using consensus partitioning and functional enrichment analysis.

## Requirements
The user needs to activiate an environment containing R packages and Python modules that need to be used. The environment is activated using:
`conda activate /data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/QLattice_env/`

## Structure
**scripts:**
Include all the scripts/jupyter notebooks needed to model the multi-omics data.\

**models:**
Include the resulting selected models of MLP in terms of binary and multiclass classifications.\

**figures:**
Include important figures after modeling.\

**results:**
Include the truly predicted DEGs and Non-DEGs from MLP, gene clusters from cola, and the resulting tables after FEA on the clusters.


## Running the Scripts
The user can first change the directory to `/data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/` and use a romote connection to the jupyter notebook from your local PC following the process in the [link](https://amber-md.github.io/pytraj/latest/tutorials/remote_jupyter_notebook).

* `functions.py`: Manually-defined functions used for preprocessing the multi-omics data and modeling.

* `MLP_agg_biclass_NonDEG.ipynb`: MLP modeling for the aggregated input gene data that have binary classes, namely DEGs and NonDEGs.

* `MLP_agg_biclass_onlyDEG.ipynb`: MLP modeling for the aggregated input DEG data that have binary classes, namely down-regulated DEGs(0) and up-regulated-DEGs(1).
s
* `MLP_agg_multiclass.ipynb`: MLP modeling for the aggregated input gene data that have 3 classes, namely NonDEGs, down-regulated DEGs and up-regulated-DEGs.

* `cola_clustering.R`: Implementation of clustering on truly predicted DEGs from MLP binary classification by using cola.

* `moonlight_FEA.R`: Functional enrichment analysis on the clusters from cola by using FEA() function in Moonlight.
