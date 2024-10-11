gc()
rm(list = ls())

scrna_pipeline_src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/bcr_data_analysis/240805_data_analysis"
source(file.path(path.to.project.src, "config.R"))
#####----------------------------------------------------------------------#####
# CONFIGURATIONS 
#####----------------------------------------------------------------------#####
analysis.round <- "1st"
chosen.seed <- 42
my_random_seed <- 42

num.dim.integration <- 30
num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
cluster.resolution <- 0.5
#####----------------------------------------------------------------------#####
##### input arguments
#####----------------------------------------------------------------------#####
all.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
all.quantiles <- c(0.99, 0.95, 0.90, 0.85)
all.integration.case <- list(
  all_samples = all.samples,
  mouse1 = c("M1", "P1"),
  mouse2 = c("M2", "P2"),
  mouse3 = c("M3", "P3")
)

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "240805_BSimons"
output.version <- "20240820"
config.version <- "default"

for (chosen.quantile in all.quantiles){
  path.to.main.input <- file.path(outdir, PROJECT, output.version, config.version)
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis", output.version, config.version)
  
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  path.to.04.output <- file.path(path.to.main.output, "04_output", sprintf("quantile_%s", chosen.quantile))
  
  all.s.obj <- list()
  for (i in seq(1, length(all.samples))){
    sample.id <- all.samples[[i]]
    path.to.03.output <- file.path(path.to.main.output, "03_output", sprintf("quantile_%s", chosen.quantile), sample.id)
    print(sprintf("reading in data from the sample %s", sample.id))
    all.s.obj[[sample.id]] <- readRDS(file.path(path.to.03.output, sprintf("GEX_sample_%s_seurat_object.rds", sample.id)))
  }
  
  ##### integrate all samples
  for (integration.case in names(all.integration.case)){
    dir.create(file.path(path.to.04.output, integration.case), showWarnings = FALSE, recursive = TRUE)
    path.to.input.s.obj <- file.path(path.to.04.output, integration.case, "s8_output", sprintf("%s.output.s8.rds", PROJECT))
    if (file.exists(path.to.input.s.obj) == FALSE){
      pca_reduction_name <- "RNA_PCA"
      umap_reduction_name <- "RNA_UMAP"
      
      integrated.samples <- all.integration.case[[integration.case]]
      
      s.obj <- merge(all.s.obj[[ integrated.samples[[1]] ]], 
                     all.s.obj[integrated.samples[2:length(integrated.samples)]])
      
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
                                                        path.to.output = file.path(path.to.04.output, integration.case), 
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
      print(sprintf("File exsited at %s", path.to.input.s.obj))
    }
  }
}


