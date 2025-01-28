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

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.04.output <- file.path(outdir, "VDJ_output", "04_output")
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

thres <- 0.85

#####---------------------------------------------------------------------------#####
##### GENERATE A BIG DATAFRAME CONTAINING ALL CLONES FROM ALL DATASETS
#####---------------------------------------------------------------------------#####
all.clone.files <- Sys.glob(file.path(outdir, "VDJ_output", "*", sprintf("VDJ_output_%s", thres), "preprocessed_files", "clonesets*.split_clones.xlsx" ))

dataset.origin <- list(
  `1st_2nd_BSimons_Datasets` = "sc",
  `220701_etc_biopsies` = "bulk",
  `240805_BSimons` = "sc",
  `240826_BSimons` = "bulk",
  `241002_BSimons` = "sc",
  `241031_BSimons` = "bulk",
  `241104_BSimons` = "sc"
)

names(all.clone.files) <- to_vec(
  for (item in all.clone.files) str_replace(str_replace(basename(item), "clonesets_", ""), ".split_clones.xlsx", "")
)

selected.cols <- c(
  "id",
  "VJseq.combi",
  "V.gene",
  "J.gene",
  "D.gene",
  "nSeqFR1",
  "nSeqCDR1",
  "nSeqFR2",
  "nSeqCDR2",
  "nSeqFR3",
  "nSeqCDR3",
  "nSeqFR4",
  "aaSeqFR1",
  "aaSeqCDR1",
  "aaSeqFR2",
  "aaSeqCDR2",
  "aaSeqFR3",
  "aaSeqCDR3",
  "aaSeqFR4",
  "VJ.len.combi",
  "targetSequences",
  "uniqueMoleculeCount")

file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv")