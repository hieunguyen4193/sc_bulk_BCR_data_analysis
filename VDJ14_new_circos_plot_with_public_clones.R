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
library(circlize)

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

sc.projects <- c("240805_BSimons_filterHT_cluster_renamed")
bulk.projects <- c("240826_BSimons")
list.of.PROJECT <- c(sc.projects, bulk.projects)

thres <- 0.85 

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
path.to.14.output <- file.path(outdir, "VDJ_output", "14_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.14.output, showWarnings = FALSE, recursive = TRUE)

path.to.tree.05.output <- file.path(outdir, "GEX_output", "05_output", sc.projects[[1]])
path.to.tree.01.output <- file.path(outdir, "tree_analysis", bulk.projects[[1]], "01_output")

match.bulkdf <- read.csv(file.path(path.to.tree.01.output, "match_public_clonedf.csv"))
match.scdf <- readxl::read_excel(file.path(path.to.tree.05.output, "240805_BSimons_filterHT_cluster_renamed_public_clone_full.xlsx")) %>%
  subset(select = -c(`...1`)) %>% 
  subset(select = c(barcode, VJcombi_CDR3_0.85, dist_to_public_clone))

