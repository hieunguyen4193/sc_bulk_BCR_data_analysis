gc()
rm(list = ls())
my_random_seed <- 42
set.seed(my_random_seed)

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/FHager_datasets/NEW_20240508"
source(file.path(path.to.project.src, "helper_functions.R"))
source(file.path(path.to.project.src, "config.R"))

#####-----------------------------------------------------------------------#####
##### install monocle
#####-----------------------------------------------------------------------#####
# if ("monocle3" %in% installed.packages()){
# remove.packages("monocle3")
# devtools::install_github("cysouw/qlcMatrix")
# install.packages("DDRTree")
# install.packages("densityClust")
# BiocManager::install("monocle.objSingleCell", update = FALSE)
# install.packages("fastICA")
# BiocManager::install("biocViews", update = FALSE)
# # remove.packages("BiocGenerics")
# BiocManager::install("HSMMSingleCell", update = FALSE)
# # install.packages("/media/hieunguyen/HD0/storage/offline_pkgs/monocle_2.30.0.tar.gz", type = "source", repos = NULL)
# install.packages("/media/hieunguyen/HD0/storage/offline_pkgs/monocle_2.30.0/monocle", type = "source", repos = NULL)
# BiocManager::install("tradeSeq", update = FALSE)
# }
#####-----------------------------------------------------------------------#####

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
library(devtools)
library(monocle)

for (dataset.name in c("1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1",
                       "2nd_dataset_removed_5_6.without_reInt.res1")){
  outdir <- "/media/hieunguyen/HD0/outdir/CRC1382"
  
  path.to.output <- file.path(outdir, "BSimons", "OUTPUT", "THESIS_OUTPUT_20231026")
  path.to.01.output <- file.path(path.to.output, "01_output", dataset.name)
  path.to.04.output <- file.path(path.to.output, "04_output", dataset.name)
  dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- readRDS(file.path(path.to.01.output, sprintf("%s.rds", dataset.name)))
  
  if (file.exists(file.path(path.to.04.output, "monocle_obj.rds")) == FALSE){
    monocle.obj <- run_monocle2(s.obj, path.to.04.output)
  } else {
    print("monocle object exists!")
    monocle.obj <- readRDS(file.path(path.to.04.output, "monocle_obj.rds"))  
  }
  
  
}

