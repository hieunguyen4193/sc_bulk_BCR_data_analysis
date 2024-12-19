gc()
rm(list = ls())
# install.packages("dplyr")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.tar.gz", type = "source", repos = NULL)
new.pkgs <- c("APackOfTheClones", "svglite", "car", "ggpubr", "ggthemes", "dplyr")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

all.integration.cases <- list(
  Dataset1_2 = c("1st_dataset", "2nd_dataset"),
  Dataset1_2_Bonn = c("1st_dataset", "2nd_dataset", "BonnData"),
  Dataset1_Bonn = c("1st_dataset", "BonnData"),
  Dataset2_Bonn = c("2nd_dataset", "BonnData")
)

num.PCA <- 25
num.PC.used.in.UMAP <- 25
num.PC.used.in.Clustering <- 25
num.dim.cluster <- 25
num.dim.integration <- 25
cluster.resolution <- 0.5
my_random_seed <- 42

# integration.case <- "Dataset1_2"
for (integration.case in names(all.integration.cases)){
  print(sprintf("Working on integration case %s", integration.case))
  path.to.04.output <- file.path(outdir, "GEX_output", "04_output", integration.case)
  dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)
  
  data.list <- list()
  for (dataset.name in all.integration.cases[[integration.case]]){
    print(sprintf("reading in %s", dataset.name))
    data.list[[dataset.name]] <- readRDS(path.to.all.s.obj[[dataset.name]])
  }
  
  s.obj.merged <- merge(x = data.list[[1]],
                        y = data.list[2 : length(data.list)])
  
  s.obj.merged <- NormalizeData(s.obj.merged) # ---> use Log Normalized
  s.obj.merged <- FindVariableFeatures(s.obj.merged, selection.method = "vst")
  s.obj.merged <- ScaleData(s.obj.merged, features = rownames(s.obj.merged))
  
  s.obj.merged <- RunPCA(s.obj.merged, npcs = num.PCA, verbose = FALSE, reduction.name="RNA_PCA")
  s.obj.merged <- RunUMAP(s.obj.merged, reduction = "RNA_PCA", 
                          dims = 1:num.PC.used.in.UMAP, reduction.name="RNA_UMAP",
                          seed.use = my_random_seed, umap.method = "uwot")
  
  s.obj.merged.1st <- s8.integration.and.clustering(s.obj = s.obj.merged, 
                                                    path.to.output = path.to.04.output, 
                                                    save.RDS.s8 = TRUE,
                                                    PROJECT = integration.case, 
                                                    num.dim.integration = num.dim.integration,
                                                    num.PCA = num.PCA,
                                                    num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                                                    num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                                                    cluster.resolution = cluster.resolution,
                                                    my_random_seed = my_random_seed,
                                                    umap.method = "uwot",
                                                    inte_pca_reduction_name = "INTE_PCA",
                                                    inte_umap_reduction_name = "INTE_UMAP",
                                                    genes.to.not.run.PCA = NULL)
}



