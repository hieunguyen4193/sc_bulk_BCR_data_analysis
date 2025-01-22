gc()
rm(list = ls())
#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

library(viridis)
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
if ("ggthemes" %in% installed.packages() == FALSE){
  install.packages("ggthemes")
}
library(ggthemes)
#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

if (input.dataset %in% c("241002_BSimons", "241104_BSimons", "BonnData")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.01.output, save.figure.type, "module_scores"), showWarnings = FALSE, recursive = TRUE)

print(sprintf("Working on dataset %s", input.dataset))
s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
