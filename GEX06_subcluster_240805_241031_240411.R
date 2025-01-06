gc()
rm(list = ls())
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
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
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

for (PROJECT in c("240805_BSimons", "241002_BSimons", "241104_BSimons")){
  path.to.06.output <- file.path(outdir, "GEX_output", "06_output", PROJECT)
  dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
  s.obj <- readRDS(path.to.all.s.obj[[PROJECT]])
  if (PROJECT == "240805_BSimons"){
    reduction.name <- "INTE_UMAP"
  } else {
    reduction.name <- "RNA_UMAP"
  }
  
  p <- DimPlot(object = s.obj, reduction = reduction.name, label = TRUE, label.box = TRUE)
  ggsave(plot = p, filename = sprintf("UMAP_clusters_%s.svg", PROJECT),
         path = path.to.06.output, dpi = 300, width = 14, height = 10)
  sub.clusters <- list(
    `240805_BSimons` = list(
      group1 = c(2, 5, 9),
      group2 = c(0, 10, 3, 4, 6, 8),
      group3 = c(1)
    ),
    `241002_BSimons` = list(
      group1 = c(1, 2, 3, 6),
      group2 = c(0, 5),
      group3 = c(4)
    ),
    `241104_BSimons` = list(
      group1 = c(2, 10, 5, 0, 6, 8),
      group2 = c(1, 4, 7),
      group3 = c(3)
    )
  )
  for (g in names(sub.clusters[[PROJECT]])){
    print(sprintf("WORKING ON PROJECT %s, group %s", PROJECT, g))
    subset.s.obj <- subset(s.obj, seurat_clusters %in% sub.clusters[[PROJECT]][[g]])
    subset.metadata <- subset.s.obj@meta.data %>% rownames_to_column("barcode")
    print(head(subset.metadata$barcode))
    write.csv(subset.metadata, file.path(outdir, 
                                         "GEX_output", 
                                         "06_output", 
                                         PROJECT, 
                                         sprintf("cell_group_%s.csv", g)))    
  }
}

