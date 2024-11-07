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

output <- list()
# for (PROJECT in c("240826_BSimons",
#                   "220701_etc_biopsies")){
for (PROJECT in c("240826_BSimons")){
  print("#####----------------------------------#####")
  print(sprintf("Working on project %s", PROJECT))
  print("#####----------------------------------#####")
  path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
  path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  output[[PROJECT]] <- run_preprocessing_all_bulk_VDJ_data(path.to.mid.output = path.to.mid.output, 
                                                           path.to.save.output = path.to.save.output,
                                                           thres = thres,
                                                           PROJECT = PROJECT,
                                                           rerun = TRUE)  
}


