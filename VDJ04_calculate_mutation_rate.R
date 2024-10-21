gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"

outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
thres <- 0.85

all.clonedf.files <- Sys.glob(file.path( outdir, "VDJ_output", "*", sprintf("VDJ_output_%s", thres), "preprocessed_files", "*clonedf.xlsx"))
all.clone.files <- Sys.glob(file.path(outdir, "VDJ_output", "*", sprintf("VDJ_output_%s", thres), "preprocessed_files", "*.split_clones.xlsx" ))

tmpdf <- readxl::read_excel(all.clonedf.files[[1]])
