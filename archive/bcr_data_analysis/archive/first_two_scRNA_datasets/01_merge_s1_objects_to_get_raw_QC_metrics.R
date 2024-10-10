#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.main.src <-  "/home/hieunguyen/CRC1382/src_2023/BSimons/CRC1382_BSimons_project"
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
PROJECT <- "BCR_dataset_1st_2nd"
outdir <- "/home/hieunguyen/CRC1382/outdir"

path.to.main.output <- file.path(outdir, "BSimons", "OUTPUT", "THESIS_OUTPUT_20231026")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.input <- file.path(outdir, "BSimons/OUTPUT/1st_round") # inputs taken from the downstream pipeline

path.to.save.output <- file.path(path.to.main.output, "01_output", "raw_seurat_objects")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

all.samples <- Sys.glob(file.path(path.to.main.input, "*"))
names(all.samples) <- basename(all.samples)

first_dataset <- to_vec(for (item in names(all.samples)) if (grepl("Sample", item) == FALSE) item)

first.data.list <- list()

for (sample in first_dataset){
  print(sprintf("Working on sample: %s", sample))
  tmp.s.obj <- readRDS(file.path(all.samples[[sample]], "s1_output", "BSimons_1st_and_2nd_dataset_added_YFP.original.output.s1.rds"))
  first.data.list[[str_replace(sample, "_1st_round", "")]] <- tmp.s.obj
}

s.obj.1st <- merge(x = first.data.list[[1]], 
                   y = unlist(first.data.list[2:length(first.data.list)]),
                   merge.data = FALSE, 
                   add.cell.ids = names(first.data.list), 
                   project = PROJECT)
saveRDS(s.obj.1st, file.path(path.to.save.output, "merge_1st_datasets.raw.rds"))


second_dataset <- to_vec(for (item in names(all.samples)) if (grepl("Sample", item) == TRUE) item)

second.data.list <- list()

for (sample in second_dataset){
  print(sprintf("Working on sample: %s", sample))
  tmp.s.obj <- readRDS(file.path(all.samples[[sample]], "s1_output", "BSimons_1st_and_2nd_dataset_added_YFP.original.output.s1.rds"))
  second.data.list[[str_replace(sample, "_1st_round", "")]] <- tmp.s.obj 
}

s.obj.2nd <- merge(x = second.data.list[[1]], 
                   y = unlist(second.data.list[2:length(second.data.list)]),
                   merge.data = FALSE, 
                   add.cell.ids = names(second.data.list), 
                   project = PROJECT)
saveRDS(s.obj.2nd, file.path(path.to.save.output, "merge_2nd_datasets.raw.rds"))


