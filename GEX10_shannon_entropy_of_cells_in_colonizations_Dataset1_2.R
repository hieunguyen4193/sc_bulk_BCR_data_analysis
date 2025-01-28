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

path.to.10.output <- file.path(outdir, "GEX_output", "10_output", input.dataset)
dir.create(path.to.10.output, showWarnings = FALSE, recursive = TRUE)


meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")

countdf <- table(meta.data$seurat_clusters, meta.data$colonization) %>%
  as.data.frame()
colnames(countdf) <- c("cluster", "colonization", "count")

num.cluster <- length(unique(s.obj$seurat_clusters))
all.shannon.entropies <- list()
for (c in unique(countdf$colonization)){
  tmpdf <- subset(countdf, countdf$colonization == c)
  tmpdf <- tmpdf %>% rowwise() %>% 
    mutate(freq = count/sum(tmpdf$count))
  tmpdf <- subset(tmpdf, tmpdf$freq != 0)
  all.shannon.entropies[[c]] <- -sum(tmpdf$freq * log2(tmpdf$freq))/log2(num.cluster)
}
all.shannon.entropies <- data.frame(all.shannon.entropies) %>% t() %>% data.frame() %>% rownames_to_column("colonization")
colnames(all.shannon.entropies) <- c("colonization", "Shannon_entropy")
writexl::write_xlsx(all.shannon.entropies, file.path(path.to.10.output, "all_Shannon_entropy_colonization_over_clusters.xlsx"))
