# Run the first downstream pipeline standard analysis
Rscript /media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/1st_2nd_BSimons_Datasets/run_pipeline.R

# Preprocess the 1st single cell VDJ dataset
Rscript /media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/1st_2nd_BSimons_Datasets/01_preprocess_01_datasets.R

# Preprocess the 2nd single cell VDJ dataset
Rscript /media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/1st_2nd_BSimons_Datasets/01_preprocess_02_datasets.R

# Run final preprocessing step
Rscript /media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/1st_2nd_BSimons_Datasets/01_run_remove_and_integrate_clusters_1st_2nd_Datasets.R
