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

library(viridis)
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

PROJECT <- "241002_241104_BSimons"

chosen.seed <- 42
my_random_seed <- 42
num.dim.integration <- 30
num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
cluster.resolution <- 0.5

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")

path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.s.obj1 <- path.to.all.s.obj[["241002_BSimons"]]
path.to.s.obj2 <- path.to.all.s.obj[["241104_BSimons"]]

pca_reduction_name <- "RNA_PCA"
umap_reduction_name <- "RNA_UMAP"

s.obj1 <- readRDS(path.to.s.obj1)
s.obj2 <- readRDS(path.to.s.obj2)

s.obj <- merge(s.obj1, s.obj2)

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
                                                  path.to.output = file.path(path.to.01.output), 
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

