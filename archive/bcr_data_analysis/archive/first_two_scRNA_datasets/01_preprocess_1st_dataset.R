#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

scrna_pipeline_src <- "/home/uk104163/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.main.src <-  "/home/uk104163/src_2023/BSimons/CRC1382_BSimons_project"
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
PROJECT <- "BCR_dataset_1st_2nd"
outdir <- "/media/outdir"

path.to.main.output <- file.path(outdir, "BSimons", "THESIS_OUTPUT_20231026")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.input <- file.path(outdir, "BSimons/OUTPUT/1st_round") # inputs taken from the downstream pipeline

path.to.save.output <- file.path(path.to.main.output, "01_output", "raw_seurat_objects")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

all.samples <- Sys.glob(file.path(path.to.main.input, "*"))
names(all.samples) <- basename(all.samples)

first_dataset <- to_vec(for (item in names(all.samples)) if (grepl("Sample", item) == FALSE) item)

first.data.list <- list()

for (sample in first_dataset){
  print(sprintf("Working on sample: %s", sample))
  tmp.s.obj <- readRDS(file.path(all.samples[[sample]], "s8a_output", "BSimons_1st_and_2nd_dataset_added_YFP.output.s8a.rds"))
  first.data.list[[str_replace(sample, "_1st_round", "")]] <- tmp.s.obj
}

s.obj.1st <- merge(x = first.data.list[[1]], 
               y = unlist(first.data.list[2:length(first.data.list)]),
               merge.data = FALSE, 
               add.cell.ids = names(first.data.list), 
               project = PROJECT)

#####----------------------------------------------------------------------#####
# PERFORM INTEGRATION FOR THE FIRST DATASET
#####----------------------------------------------------------------------#####

# configurations
chosen.seed <- 42
num.PC.used.in.UMAP <- 25 
num.dim.integration <- 25
num.PCA <- 25
num.dim.cluster <- 25
cluster.resolution <- 0.5
num.PC.used.in.Clustering <- 25

if (file.exists(file.path(path.to.save.output, "merged_all_first_dataset_BCR.rds")) == FALSE){
  chosen.assay <- "RNA"
  DefaultAssay(s.obj.1st) <- chosen.assay
  
  s.obj.1st <- NormalizeData(s.obj.1st) # ---> use Log Normalized
  s.obj.1st <- FindVariableFeatures(s.obj.1st, selection.method = "vst")
  s.obj.1st <- ScaleData(s.obj.1st, features = rownames(s.obj.1st))
  
  s.obj.1st <- RunPCA(s.obj.1st, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay))
  s.obj.1st <- RunUMAP(s.obj.1st, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                             reduction.name=sprintf("%s_UMAP", chosen.assay), seed.use = chosen.seed)
  s.obj.1st <- RunTSNE(s.obj.1st, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                             reduction.name=sprintf("%s_TSNE", chosen.assay), seed.use = chosen.seed)
  
  s.obj.1st <- FindNeighbors(s.obj.1st, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
  
  s.obj.1st <- FindClusters(s.obj.1st, resolution = cluster.resolution, random.seed = chosen.seed)
  
  saveRDS(object = s.obj.1st, file = file.path(path.to.save.output, "merged_all_first_dataset_BCR.rds")) 
} else {
  s.obj.1st <- readRDS(file.path(path.to.save.output, "merged_all_first_dataset_BCR.rds"))
}

if (file.exists(file.path(path.to.save.output, "merged_all_first_dataset_BCR.integrated.rds")) == FALSE){
  s.obj.1st.integrated <- s8.integration.and.clustering(s.obj.1st, 
                                                        path.to.save.output, 
                                                        FALSE,
                                                        PROJECT, 
                                                        num.dim.integration,
                                                        num.PCA,
                                                        num.PC.used.in.UMAP,
                                                        num.PC.used.in.Clustering,
                                                        cluster.resolution = cluster.resolution,
                                                        my_random_seed = 42,
                                                        umap.method = "uwot",
                                                        genes.to.not.run.PCA = NULL,
                                                        inte_pca_reduction_name = "INTE_PCA", 
                                                        inte_umap_reduction_name = "INTE_UMAP")
  saveRDS(s.obj.1st.integrated, file.path(path.to.save.output, "merged_all_first_dataset_BCR.integrated.rds"))
} else {
  s.obj.1st.integrated <- readRDS(file.path(path.to.save.output, "merged_all_first_dataset_BCR.integrated.rds"))
}

#####----------------------------------------------------------------------#####
# FIND MARKER GENES
#####----------------------------------------------------------------------#####
DefaultAssay(s.obj.1st.integrated) <- "RNA"
if (file.exists(file.path(path.to.save.output, "clusterMarkers_1st_dataset.rds")) == FALSE){
  cluster.markers.1st.dataset <- FindAllMarkers(object = s.obj.1st.integrated, test.use = "wilcox", assay = "RNA")
  saveRDS(cluster.markers.1st.dataset, file.path(path.to.save.output, "clusterMarkers_1st_dataset.rds"))
  cluster.markers.1st.dataset <- subset(cluster.markers.1st.dataset, (cluster.markers.1st.dataset$p_val_adj < 0.05) 
                                        & (cluster.markers.1st.dataset$avg_log2FC > 0))
  saveRDS(cluster.markers.1st.dataset, file.path(path.to.save.output, "clusterMarkers_1st_dataset.filtered.rds"))
} else {
  cluster.markers.1st.dataset <- readRDS(file.path(path.to.save.output, "clusterMarkers_1st_dataset.rds"))
}



