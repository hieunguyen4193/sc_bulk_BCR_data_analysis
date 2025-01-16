#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

# if (packageVersion("Matrix") != "1.5.4.1"){
#   install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", type = "source", repos = NULL)
#   .rs.restartR()
# }
# 
# if (packageVersion("ggplot2") != "3.4.4"){
#   install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", type = "source", repos = NULL)
#   .rs.restartR()
# }

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "240805_BSimons_filterHT_cluster_renamed"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

s.obj <- readRDS(path.to.all.s.obj[[PROJECT]])

# DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
writexl::write_xlsx(meta.data, file.path(path.to.main.output, "metadata.xlsx"))