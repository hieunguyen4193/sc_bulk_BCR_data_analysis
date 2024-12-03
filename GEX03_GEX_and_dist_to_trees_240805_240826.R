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

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
mouseid <- "m1"

path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
input.dataset <- "240805_BSimons"

if (input.dataset %in% c("241104_BSimons", "241002_BSimons")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

path.to.03.output <- file.path(outdir, "GEX_output", "03_output", input.dataset, mouseid)
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

print(sprintf("Working on dataset %s", input.dataset))
s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
DefaultAssay(s.obj) <- "RNA"
sample.list <- list(
  m1 = c("M1", "P1"),
  m2 = c("M2", "P2"),
  m3 = c("M3", "P3")
)

s.obj <- subset(s.obj, name %in% sample.list[[mouseid]])
path.to.distdf <- file.path(outdir, sprintf("tree_analysis/06_output/240826_BSimons_240805_BSimons/%s/scdistdf.csv", mouseid))

distdf <- read.csv(Sys.glob(path.to.distdf)) %>%
  subset(select = -c(X))

distdf.min <- distdf %>%
  group_by(barcode) %>%
  summarise(minDistTree = min(min_dist_to_a_tree))
  
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
umapdf <- s.obj@reductions[[reduction.name]]@cell.embeddings %>% as.data.frame() %>%
  rownames_to_column("barcode")

colnames(umapdf) <- c("barcode", "UMAP1", "UMAP2")
meta.data <- merge(meta.data, distdf.min, by.x = "barcode", by.y = "barcode", all.x = TRUE)
meta.data <- meta.data %>% column_to_rownames("barcode") %>%
  rowwise(with.dist = ifelse(is.na(gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
mouseid <- "m1"

path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
input.dataset <- "240805_BSimons"

if (input.dataset %in% c("241104_BSimons", "241002_BSimons")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

path.to.03.output <- file.path(outdir, "GEX_output", "03_output", input.dataset, mouseid)
dir.create(path.to.03.output, showWarnings = FALSE, recursive = TRUE)

print(sprintf("Working on dataset %s", input.dataset))
s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
DefaultAssay(s.obj) <- "RNA"
sample.list <- list(
  m1 = c("M1", "P1"),
  m2 = c("M2", "P2"),
  m3 = c("M3", "P3")
)

s.obj <- subset(s.obj, name %in% sample.list[[mouseid]])
path.to.distdf <- file.path(outdir, sprintf("tree_analysis/06_output/240826_BSimons_240805_BSimons/%s/scdistdf.csv", mouseid))

distdf <- read.csv(Sys.glob(path.to.distdf)) %>%
  subset(select = -c(X))

distdf.min <- distdf %>%
  group_by(barcode) %>%
  summarise(minDistTree = min(min_dist_to_a_tree))
  
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
umapdf <- s.obj@reductions[[reduction.name]]@cell.embeddings %>% as.data.frame() %>%
  rownames_to_column("barcode")

colnames(umapdf) <- c("barcode", "UMAP1", "UMAP2")
meta.data <- merge(meta.data, distdf.min, by.x = "barcode", by.y = "barcode", all.x = TRUE)

meta.data <- meta.data %>%
  rowwise() %>%
  mutate(with.dist = ifelse(is.na(minDistTree) == TRUE, "no", "yes"))
meta.data <- meta.data %>% 
  column_to_rownames("barcode") 
meta.data <- meta.data[row.names(s.obj@meta.data), ]

s.obj <- AddMetaData(object = s.obj, 
                     col.name = "minDistTree", 
                     metadata = meta.data$minDistTree)

feature.dist.plot <- FeaturePlot(object = s.obj, 
            reduction = reduction.name, 
            label = TRUE, 
            features = "minDistTree", 
            pt.size = 2, 
            order = TRUE) +
  scale_color_gradient(low = "gray28", high = "red")

violin.dist.plot <- meta.data %>% 
  ggplot(aes(x = seurat_clusters, y = dist, fill = seurat_clusters)) + 
  geom_boxplot()
) 
