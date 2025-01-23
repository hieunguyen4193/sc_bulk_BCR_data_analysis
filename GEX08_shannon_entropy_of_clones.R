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
if ("ggthemes" %in% installed.packages() == FALSE){
  install.packages("ggthemes")
}
library(ggthemes)
#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

input.dataset <- "Dataset1_2"

if (input.dataset %in% c("241002_BSimons", "241104_BSimons", "BonnData")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])

if (input.dataset == "Dataset1_2"){
  to.run.clusters <- c("seurat_clusters", "colonization")
  colonizations <- list( MM9_Ecoli = c("17_MM9_Ecoli", 
                                       "20_MM9_Ecoli"),
                         MM9_Ecoli_SPF = c("21_MM9_Ecoli_SPF"),
                         MM9 = c("MM9_S2",
                                 "MM9_S4"),
                         MM9_SPF = c( "MM9_SPF_S3",
                                      "MM9_SPF_S9"),
                         dataset2 = c("Sample_132",
                                      "Sample_133")
  )
  colonizationdf <- data.frame(SampleID = unlist(colonizations)) %>%
    rownames_to_column("colonization") %>%
    rowwise() %>%
    mutate(colonization = ifelse(colonization != "MM9_Ecoli_SPF", 
                                 paste(str_split(colonization, "")[[1]][1:nchar(colonization)- 1], collapse = ""), 
                                 "MM9_Ecoli_SPF"))
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(colonization = subset(colonizationdf, colonizationdf$SampleID == name)$colonization) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, col.name = "colonization", metadata = meta.data$colonization)
} else {
  to.run.clusters <- c("seurat_clusters")
}

for (clone.name in c("VJcombi_CDR3_0.85")){
  for (cluster.name in to.run.clusters){
    
    path.to.08.output <- file.path(outdir, "GEX_output", "08_output", input.dataset, clone.name, cluster.name)
    dir.create(path.to.08.output, showWarnings = FALSE, recursive = TRUE)
    
    
    N <- length(s.obj@meta.data[[cluster.name]])
    print(sprintf("Working on dataset %s", input.dataset))
    
    clonedf <- data.frame(clone = unique(s.obj@meta.data[[clone.name]])) %>%
      subset(is.na(clone) == FALSE) %>%
      rowwise() %>%
      mutate(total.count = subset(s.obj@meta.data, s.obj@meta.data$VJcombi_CDR3_0.85 == clone) %>% nrow())
    
    clonedf$shannon.entropy <- unlist(lapply(clonedf$clone, function(x){
      tmpdf <- subset(s.obj@meta.data, s.obj@meta.data[[clone.name]] == x)
      if (nrow(tmpdf) >= 10){
        countdf <- table(tmpdf[[cluster.name]]) %>% data.frame()
        colnames(countdf) <- c(cluster.name, "count")
        countdf <- countdf %>% rowwise() %>%
          mutate(p = count/sum(countdf$count))
        countdf <- subset(countdf, countdf$count != 0)
        shannon.entropy <- - sum(countdf$p * log2(countdf$p))/log2(N)  
      } else {
        shannon.entropy <- NA
      }
      return(shannon.entropy)
    }))
    
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
    meta.data <- merge(meta.data, clonedf, by.x = clone.name, by.y = "clone", all.x = TRUE) 
    meta.data <- meta.data %>%
      column_to_rownames("barcode")
    meta.data <- meta.data[row.names(s.obj@meta.data), ]
    s.obj <- AddMetaData(object = s.obj, col.name = sprintf("shannon.entropy.%s", cluster.name), metadata = meta.data$shannon.entropy)
    
    shannon.umap.plot <- FeaturePlot(object = s.obj, features = c(sprintf("shannon.entropy.%s", cluster.name)), reduction = reduction.name) &
      scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
    subset.s.obj <- subset(s.obj, cells = row.names(subset(s.obj@meta.data, is.na(s.obj@meta.data[[sprintf("shannon.entropy.%s", cluster.name)]]) == FALSE)))
    shannon.vln.plot <- VlnPlot(object = subset.s.obj, group.by = "seurat_clusters", features = c(sprintf("shannon.entropy.%s", cluster.name)))
    
    ggsave(plot = shannon.umap.plot, filename = sprintf("UMAP_Shannon_entropy.svg"), path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
    ggsave(plot = shannon.vln.plot, filename = sprintf("ViolinPlot_Shannon_entropy.svg"), path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
  }
}
