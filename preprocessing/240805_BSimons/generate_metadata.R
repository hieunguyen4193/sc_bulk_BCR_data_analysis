#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

if (packageVersion("Matrix") != "1.5.4.1"){
  install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", type = "source", repos = NULL)
  .rs.restartR()
}

if (packageVersion("ggplot2") != "3.4.4"){
  install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", type = "source", repos = NULL)
  .rs.restartR()
}

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "240805_BSimons"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

s.obj <- readRDS(file.path(outdir, PROJECT, "/data_analysis/04_output/quantile_0.85/all_samples/s8_output/240805_BSimons.output.s8.rds"))

DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
writexl::write_xlsx(meta.data, file.path(path.to.main.output, "metadata.xlsx"))