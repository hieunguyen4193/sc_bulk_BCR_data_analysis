#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis"
PROJECT <- "241002_BSimons"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

s.obj <- readRDS(file.path(outdir, PROJECT, "data_analysis/03_output/quantile_0.85/PP3/GEX_sample_PP3_seurat_object.rds"))

DimPlot(object = s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
writexl::write_xlsx(meta.data, file.path(path.to.main.output, "metadata.xlsx"))