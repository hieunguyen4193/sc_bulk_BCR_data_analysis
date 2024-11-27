gc()
rm(list = ls())

new.pkgs <- c("APackOfTheClones", "svglite")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}

library(APackOfTheClones)

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
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

input.dataset <- "1st_dataset"

print(sprintf("Working on dataset %s", input.dataset))
s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
DefaultAssay(s.obj) <- "RNA"

path.to.02.output <- file.path(outdir, "GEX_output", "02_output", input.dataset)
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)
reduction.name <- "INTE_UMAP"

s.obj <- RunAPOTC(seurat_obj = s.obj, reduction_base = reduction.name, clonecall = "CTstr")
