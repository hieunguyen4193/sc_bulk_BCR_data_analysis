gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/00_helper_functions.R")
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
new.pkgs <- c("svglite")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }   
}

new.pkgs <- c("msa")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() ==  FALSE){
    BiocManager::install(pkg, update = FALSE)
  }
}
library("msa")

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis"
PROJECT <- "1st_2nd_BSimons_Datasets"
thres <- 0.85

path.to.VDJ.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres))
path.to.save.output <- file.path(path.to.VDJ.output, "preprocessed_files")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
