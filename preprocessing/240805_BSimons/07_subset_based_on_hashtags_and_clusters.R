gc()
rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241002_BSimons"
source(file.path(path.to.project.src, "config.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"
chosen.seed <- 42
num.dim.integration <- 25 
num.PCA <- 25
num.dim.cluster <- 25
num.PC.used.in.Clustering <- 25
num.PC.used.in.UMAP <- 25
my_random_seed <- 42

all.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
all.quantiles <- c(0.99, 0.95, 0.90, 0.85)
all.integration.case <- list(
  all_samples = all.samples,
  mouse1 = c("M1", "P1"),
  mouse2 = c("M2", "P2"),
  mouse3 = c("M3", "P3")
)

#####----------------------------------------------------------------------#####
##### input arguments
#####----------------------------------------------------------------------#####

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "240805_BSimons"
output.version <- "20240820"
config.version <- "default"
chosen.quantile <- 0.85
integration.case <- "all_samples"

# outdir <- params$outdir
# PROJECT <- params$PROJECT
# chosen.quantile <- params$chosen.quantile
# integration.case <- params$integration.case

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.04.output <- file.path(path.to.main.output, "04_output", sprintf("quantile_%s", chosen.quantile))
path.to.05.output <- file.path(path.to.main.output, "05_output", sprintf("quantile_%s", chosen.quantile), integration.case)
path.to.07.output <- file.path(path.to.main.output, "07_output", sprintf("quantile_%s", chosen.quantile), integration.case)
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

path.to.input.s.obj <- file.path(path.to.04.output, integration.case, "s8_output", sprintf("%s.output.s8.rds", PROJECT))
s.obj <- readRDS(path.to.input.s.obj)

keep.hashtags <- list(
  M1 = c("HT1", "HT2", "HT3"),
  P1 = c("HT1", "HT2", "HT3", "HT5"),
  M2 = c("HT1", "HT2", "HT3"),
  P2 = c("HT1", "HT2", "HT3", "HT5"),          
  M3 = c("HT1", "HT2", "HT3"),
  P3 = c("HT2", "HT3", "HT4", "HT6")
)

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
keep.cells <- c()
for (sample.id in unique(meta.data$name)){
  tmp.metadata <- subset(meta.data, meta.data$name == sample.id)
  tmp.metadata <- subset(tmp.metadata, tmp.metadata$HTO_classification %in% keep.hashtags[[sample.id]])
  keep.cells <- c(keep.cells, tmp.metadata$barcode)
}

if (file.exists(file.path(path.to.07.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT))) == FALSE){
  s.obj <- subset(s.obj, cells = keep.cells)
  s.obj <- subset(s.obj, seurat_clusters %in% c(8, 11) == FALSE)
  
  for (sample.id in unique(s.obj$name)){
    print(sprintf("List of HT in this sample %s: %s", 
                  sample.id, 
                  paste0(unique(subset(s.obj@meta.data, s.obj@meta.data$name == sample.id)$HTO_classification), collapse = ", ")))
  }
  
  s.obj <- DietSeurat(s.obj)
  DefaultAssay(s.obj) <- "RNA"
  pca_reduction_name <- "RNA_PCA"
  umap_reduction_name <- "RNA_UMAP"
  
  s.obj <- NormalizeData(s.obj) # ---> use Log Normalized
  s.obj <- FindVariableFeatures(s.obj, selection.method = "vst")
  s.obj <- ScaleData(s.obj, features = rownames(s.obj))
  
  s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=pca_reduction_name)
  s.obj <- RunUMAP(s.obj, reduction = pca_reduction_name, 
                   dims = 1:num.PC.used.in.UMAP, reduction.name=umap_reduction_name,
                   seed.use = my_random_seed, umap.method = "uwot")
  # clustering 
  s.obj <- FindNeighbors(s.obj, reduction = pca_reduction_name, dims = 1:num.PC.used.in.Clustering)
  s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)
  
  s.obj.integrated <- s8.integration.and.clustering(s.obj = s.obj, 
                                                    path.to.output = file.path(path.to.07.output), 
                                                    save.RDS.s8 = TRUE,
                                                    PROJECT = PROJECT, 
                                                    num.dim.integration = num.dim.integration,
                                                    num.PCA = num.PCA,
                                                    num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                    num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                    cluster.resolution = cluster.resolution,
                                                    my_random_seed = 42,
                                                    umap.method = "uwot",
                                                    genes.to.not.run.PCA = NULL,
                                                    inte_pca_reduction_name = "INTE_PCA", 
                                                    inte_umap_reduction_name = "INTE_UMAP",
                                                    with.TSNE = FALSE,
                                                    k.filter = 200)
} else {
  print(sprintf("Reading existing data from %s", file.path(path.to.07.output, "s8_output")))
  s.obj <- readRDS(file.path(path.to.07.output, "s8_output", sprintf("%s.output.s8.rds", PROJECT)))
}

if (file.exists(file.path(path.to.07.output, 
                          "s8_output", 
                          sprintf("%s.renamedClusters.output.s8.rds", PROJECT))) == FALSE){
  before.p <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE)
  change.cluster.name <- list(
    `cluster_7` = 1, 
    `cluster_2` = 2, 
    `cluster_6` = 3,
    `cluster_5` = 4,
    `cluster_0` = 5, 
    `cluster_1` = 6,
    `cluster_8` = 7,
    `cluster_4` = 8,
    `cluster_3` = 9
  )
  meta.data.original <- s.obj@meta.data
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(seurat_clusters = change.cluster.name[[sprintf("cluster_%s", seurat_clusters)]]) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$seurat_clusters, col.name = "seurat_clusters")
  after.p <- DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, group.by = "seurat_clusters")
  saveRDS(s.obj, file.path(path.to.07.output, "s8_output", sprintf("%s.renamedClusters.output.s8.rds", PROJECT)))
} 

