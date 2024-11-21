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
path.to.storage <- "/media/hieunguyen/GSHD_HN01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
PROJECT <- "241002_BSimons"

path.to.07.output <- file.path(outdir, "VDJ_output", "07_output", PROJECT)
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(path.to.all.s.obj[[PROJECT]])
DefaultAssay(s.obj) <- "RNA"
all.genes <- row.names(s.obj)

all.cells <- GetAssayData(object = s.obj, assay = "RNA", slot = "data")["YFP", ]

yfp.cells <- all.cells[all.cells != 0] %>% names()
non.yfp.cells <- all.cells[all.cells == 0] %>% names()

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  rowwise() %>%
  mutate(YFP = ifelse(barcode %in% yfp.cells, "yes", "no")) %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$YFP, col.name = "YFP")

yfp.de.markers <- FindMarkers(object = s.obj, group.by = "YFP", ident.1 = "yes", ident.2 = "no", assay = "RNA", test.use = "wilcox",
                              features = setdiff(all.genes, "YFP"))
saveRDS(yfp.de.markers, file.path(path.to.07.output, "yfp.de.markers.raw.rds"))

