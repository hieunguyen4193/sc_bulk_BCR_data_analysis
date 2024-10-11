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
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis"
PROJECT <- "20240805_BSimons"
path.to.main.output <- file.path(outdir, PROJECT, "s8_output")
s.obj <- readRDS(file.path(path.to.main.output, 
                               "01_output", 
                               "2nd_dataset_removed_5_6.without_reInt.res1", 
                               "2nd_dataset_removed_5_6.without_reInt.res1.rds"))
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")

writexl::write_xlsx(meta.data, file.path(path.to.main.output, "metadata.xlsx"))