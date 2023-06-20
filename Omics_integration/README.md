## Introduction 
This project is to investigate the effect of omics data types on differential gene expression of breast cancer samples compared with normal samples. This folder is generated with the purpose to integrate the multi-omics data types, and create the input data sets ready for the following modeling process.

## Requirements
The user needs to activiate an environment containing R packages and Python modules that need to be used. The environment is activated using:
`conda activate /data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/QLattice_env/`

## Structure
**scripts:**
Include all the scripts needed to produce the integrated omics data\

**data:**
Include the data needed for preprocessing and building the individual gene-level omics data\

**results:**
* `omics_all_samples.csv`: the input data with all the common samples. The genes include NonDEGs, up-regulated and down-regulated DEGs.
* `Single_omics_files`: the per-sample-based data within each subtype folder\s

## Running the Scripts
The user should be in the `/data/user/shared_projects/moonlight_bc_paper/moonlight_QLattice/Omics_integration/scripts` directory and run the scripts from there. 

* `TCGAanalyze_DEA_revised.R`: function used to do the differential expression analysis

* `TCGA_DEA.r`: changes the threshold of the DEG selection that aligns with Moonlihgt project. The output is under `../data/BRCA_dataDEGs.rda`

* `data_observe.R`: checks for the omics data before preprocessing and finds common samples. The output is under `../data/common_samples.csv`

* `GISTIC_prepare.R`: prepares the data needed to run GISTICS web-based tool for the CNV data. The output data are `../data/CNV_segment_file.txt` and `../data/CNV_marker_file.txt`. The reference file downloaded from CGD webpage for running GISTIC is `../data/snp6.na35.remap.hg38.subset.txt`. The web-based tool for running GISTIC is `https://cloud.genepattern.org/gp/pages/index.jsf`.

* `cscape_functions.R`: functions used to annotate the MAF file with CScape-somatic scores for each mutation

* `data_wrangle_functions.R`: functions to be used in omics data wrangling

* `annotate_MAF.R`: applies the functions to annotate the MAF fie. The original MAF file is `../data/mutations.csv`. The CScape-somatic score data used in the script are `../data/css_coding.vcf.gz` and `../data/css_noncoding.vcf.gz`.

* `omics_wrangle.R`: applies the wrangling function to each omics data. To run this script, one can run the `run.sh` file with command line by using:
`bash run.sh`.
The probe annotation file used here for Methylation wrangling is `../data/HM450.hg38.manifest.gencode.v36.tsv`, and the output is saved under `../data/wrangled_data.RData`

* `omics_merge_all.R`:merges the wrangled omics data together to produce the input with all the common samples. 

* `omics_merge_single.R`:merges the wrangled omics data together to produce the input for each single sample.

* `separate_files.R`: separate omics files for each sample based on subtype information, and move the files to the corresponding subtype folders. 

