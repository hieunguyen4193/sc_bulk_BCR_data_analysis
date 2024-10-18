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
PROJECT <- "240826_BSimons"

path.to.save.fasta <- file.path(outdir, "FASTA_output", PROJECT,  sprintf("VDJ_output_%s", thres))
dir.create(path.to.save.fasta, showWarnings = FALSE, recursive = TRUE)

path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

clone.output <- run_preprocessing_all_bulk_VDJ_data(
  path.to.mid.output = path.to.mid.output, 
  path.to.save.output = path.to.save.output,
  thres = thres,
  PROJECT = PROJECT)  

clonesets <- clone.output$clonesets
input.clone <- "IGHV14-3*01_IGHJ2*01_36_1"
ref.gene <- "10x"
ref.gene.config <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/ref_gene_config.R"
