#####----------------------------------------------------------------------#####
#
# 01: PREPROCESSING THE FIRST AND SECOND DATASETS
#
# trnguyen@ukaachen.de
#####----------------------------------------------------------------------#####

##### clean up #####
# gc()
# rm(list = ls())

to_lower_gene_symbol <- function(s){
  list_string <- unlist(strsplit(s, ""))
  
  list_string[1] <- toupper(list_string[1])
  list_string[2:length(list_string)] <- tolower(list_string[2:length(list_string)])
  new.string <- paste(list_string, collapse = "")
  return(new.string)
}


scrna_pipeline_src <- "/home/uk104163/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.main.src <-  "/home/uk104163/src_2023/BSimons/BSimons_thesis_data"
#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
PROJECT <- "BCR_dataset_1st_2nd"

outdir <- "/media/outdir"

path.to.main.output <- file.path(outdir, "BSimons", "THESIS_OUTPUT_20231026")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.main.input <- file.path(outdir, "BSimons", "OUTPUT", "1st_round")

path.to.save.output <- file.path(path.to.main.output, "01_output", "raw_seurat_objects")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

all.samples <- Sys.glob(file.path(path.to.main.input, "*"))
names(all.samples) <- basename(all.samples)

second_dataset <- to_vec(for (item in names(all.samples)) if (grepl("Sample", item) == TRUE) item)

second.data.list <- list()

for (sample in second_dataset){
  print(sprintf("Working on sample: %s", sample))
  tmp.s.obj <- readRDS(file.path(all.samples[[sample]], "s8a_output", "BSimons_1st_and_2nd_dataset_added_YFP.output.s8a.rds"))
  second.data.list[[str_replace(sample, "_1st_round", "")]] <- tmp.s.obj 
}

s.obj.2nd <- merge(x = second.data.list[[1]], 
                   y = unlist(second.data.list[2:length(second.data.list)]),
                   merge.data = FALSE, 
                   add.cell.ids = names(second.data.list), 
                   project = PROJECT)

#####----------------------------------------------------------------------#####
# PERFORM INTEGRATION FOR THE SECOND DATASET
#####----------------------------------------------------------------------#####

# configurations
chosen.seed <- 42
num.PC.used.in.UMAP <- 25 
num.dim.integration <- 25
num.PCA <- 25
num.dim.cluster <- 25
cluster.resolution <- 0.5
num.PC.used.in.Clustering <- 25

if (file.exists(file.path(path.to.save.output, "merged_all_second_dataset_BCR.rds")) == FALSE){
  chosen.assay <- "RNA"
  DefaultAssay(s.obj.2nd) <- chosen.assay
  
  s.obj.2nd <- NormalizeData(s.obj.2nd) # ---> use Log Normalized
  s.obj.2nd <- FindVariableFeatures(s.obj.2nd, selection.method = "vst")
  
  remove.YFP.VariableFeatures <- to_vec(for (item in VariableFeatures(s.obj.2nd)) if (item != "YFP") item)
  
  s.obj.2nd <- ScaleData(s.obj.2nd, features = to_vec(for (item in rownames(s.obj.2nd)) if(item != "YFP") item))
  
  s.obj.2nd <- RunPCA(s.obj.2nd, npcs = num.PCA, verbose = FALSE, reduction.name=sprintf("%s_PCA", chosen.assay), 
                      features = remove.YFP.VariableFeatures)
  s.obj.2nd <- RunUMAP(s.obj.2nd, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                       reduction.name=sprintf("%s_UMAP", chosen.assay), seed.use = chosen.seed)
  s.obj.2nd <- RunTSNE(s.obj.2nd, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PCA, 
                       reduction.name=sprintf("%s_TSNE", chosen.assay), seed.use = chosen.seed)
  
  s.obj.2nd <- FindNeighbors(s.obj.2nd, reduction = sprintf("%s_PCA", chosen.assay), dims = 1:num.PC.used.in.Clustering)
  
  s.obj.2nd <- FindClusters(s.obj.2nd, resolution = cluster.resolution, random.seed = chosen.seed)
  
  saveRDS(object = s.obj.2nd, file = file.path(path.to.save.output, "merged_all_second_dataset_BCR.rds")) 
} else {
  s.obj.2nd <- readRDS(file.path(path.to.save.output, "merged_all_second_dataset_BCR.rds"))
}

if (file.exists(file.path(path.to.save.output, "merged_all_second_dataset_BCR.integrated.rds")) == FALSE){
  s.obj.2nd.integrated <- s8.integration.and.clustering(s.obj.2nd, 
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
  saveRDS(s.obj.2nd.integrated, file.path(path.to.save.output, "merged_all_second_dataset_BCR.integrated.rds"))
} else {
  s.obj.2nd.integrated <- readRDS(file.path(path.to.save.output, "merged_all_second_dataset_BCR.integrated.rds"))
}

count.YFP <- GetAssayData(s.obj.2nd, slot = "counts")["YFP", ]
count.YFP <- count.YFP[count.YFP != 0]
yfp.cells <- names(count.YFP)

#####----------------------------------------------------------------------#####
# FIND MARKER GENES
#####----------------------------------------------------------------------#####
DefaultAssay(s.obj.2nd.integrated) <- "RNA"
if (file.exists(file.path(path.to.save.output, "clusterMarkers_2nd_dataset.rds")) == FALSE){
  cluster.markers.2nd.dataset <- FindAllMarkers(object = s.obj.2nd.integrated, test.use = "wilcox", assay = "RNA")
  saveRDS(cluster.markers.2nd.dataset, file.path(path.to.save.output, "clusterMarkers_2nd_dataset.rds"))
  cluster.markers.2nd.dataset <- subset(cluster.markers.2nd.dataset, (cluster.markers.2nd.dataset$p_val_adj < 0.05) 
                                        & (cluster.markers.2nd.dataset$avg_log2FC > 0))
  saveRDS(cluster.markers.2nd.dataset, file.path(path.to.save.output, "clusterMarkers_2nd_dataset.filtered.rds"))
  
} else {
  cluster.markers.2nd.dataset <- readRDS(file.path(path.to.save.output, "clusterMarkers_2nd_dataset.rds"))
}



