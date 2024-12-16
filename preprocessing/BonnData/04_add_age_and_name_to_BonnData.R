gc()
rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240805_BSimons"
source(file.path(path.to.project.src, "config.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "BonnData"
sample.id <- "BonnData"
chosen.quantile <- 0.85

path.to.s.obj <- file.path(outdir, 
                           PROJECT, 
                           "data_analysis", 
                           "03_output", 
                           sprintf("quantile_%s", chosen.quantile), 
                           PROJECT, 
                           "GEX_sample_BonnData_seurat_object.rds")
s.obj <- readRDS(path.to.s.obj)

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  rowwise() %>%
  mutate(name = str_replace(HTO_classification, "CD45-", "")) %>%
  mutate(age = str_split(name, "-")[[1]][[1]]) %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, col.name = "age", metadata = meta.data$age)
s.obj <- AddMetaData(object = s.obj, col.name = "name", metadata = meta.data$name)

saveRDS(s.obj, file.path(outdir, 
                         PROJECT, 
                         "data_analysis", 
                         "03_output", 
                         sprintf("quantile_%s", chosen.quantile), 
                         PROJECT, 
                         "GEX_sample_BonnData_seurat_object.addedSampleID.rds"))