# Pipelines

## Downstream analysis pipelines
- Run `run_pipeline.R` to run the general downstream analysis pipelines on all samples from the first and second datasets. This pipeline uses `Cellranger` outputs as inputs.

- Run the function `01_preprocess_01_datasets.R` to perform preprocessing steps and integration of all samples in the first dataset.

- Run the function `01_preprocess_02_datasets.R` to perform preprocessing steps and integration of all samples in the second dataset.

- Run `01_run_remove_and_integrate_clusters_1st_2nd_Datasets.R` to remove some clusters and re-integrate the all samples in the 1st and 2nd datasets.

## Clone data analysis