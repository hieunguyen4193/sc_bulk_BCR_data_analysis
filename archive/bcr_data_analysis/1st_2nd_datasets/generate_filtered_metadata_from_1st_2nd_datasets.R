#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "1st_2nd_BSimons_Datasets"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
s.obj.1st <- readRDS(file.path(path.to.main.output, 
                               "01_output",
                               "1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1", 
                               "1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1.rds"))
metadata.1st <- s.obj.1st@meta.data %>% rownames_to_column("barcode")

s.obj.2nd <- readRDS(file.path(path.to.main.output, 
                               "01_output", 
                               "2nd_dataset_removed_5_6.without_reInt.res1", 
                               "2nd_dataset_removed_5_6.without_reInt.res1.rds"))
metadata.2nd <- s.obj.2nd@meta.data %>% rownames_to_column("barcode")

writexl::write_xlsx(metadata.1st, file.path(path.to.main.output, "01_output", "metadata_1st_dataset.xlsx"))
writexl::write_xlsx(metadata.2nd, file.path(path.to.main.output, "01_output", "metadata_2nd_dataset.xlsx"))