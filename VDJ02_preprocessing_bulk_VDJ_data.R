gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/00_helper_functions.R")
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis"

run_preprocessing_bulk_VDJ_data(outdir = outdir, PROJECT = "220701_etc_biopsies", thres = 0.85, ref.gene = "IMGT")
run_preprocessing_bulk_VDJ_data(outdir = outdir, PROJECT = "240826_BSimons", thres = 0.85, ref.gene = "IMGT")

#####----------------------------------------------------------------------#####
##### read REFERENCE GENES
#####----------------------------------------------------------------------#####
# s.V.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHV.fasta")) 
# names(s.V.genes) <- lapply(names(s.V.genes), function(x){
#   str_split(x, "[|]")[[1]][[2]]
# })
# 
# s.J.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHJ.fasta")) 
# names(s.J.genes) <- lapply(names(s.J.genes), function(x){
#   str_split(x, "[|]")[[1]][[2]]
# })