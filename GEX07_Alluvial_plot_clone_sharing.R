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
library(ggpubr)
library(ggthemes)
library(APackOfTheClones)
library("gridExtra")

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
path.to.all.s.obj <- path.to.all.s.obj[setdiff(names(path.to.all.s.obj), c("BonnData"))]

sc.projects.with.ht <- c("240805_BSimons_filterHT_cluster_renamed")
# name_or_sampleHT <- "sample_ht"
name_or_sampleHT <- "name"

if (name_or_sampleHT == "sample_ht"){
  path.to.all.s.obj <- path.to.all.s.obj[sc.projects.with.ht]
}
clone.name <- "VJcombi_CDR3_0.85"
dataset.name <- "240805_BSimons_filterHT_cluster_renamed"
save.dev <- "png"

##### sample.list for each mouse:
if (grepl("240805_BSimons", dataset.name) == TRUE){
  sample.list <- list(
    M_samples = c("M1", "M2", "M3"),
    P_samples = c("P1", "P2", "P3")
  )
} else if (grepl("241002_BSimons", dataset.name) == TRUE){
  sample.list <- list(
    m3 = c("PP3")
  )
} else if (grpel("241104_BSimons", dataset.name) == TRUE){
  sample.list <- list(
    m7 = c("PP7")
  )
}

print(sprintf("Working on dataset %s", dataset.name))

s.obj <- readRDS(path.to.all.s.obj[[dataset.name]])
if (dataset.name %in% sc.projects.with.ht){
  meta.data <- s.obj@meta.data %>% 
    rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(sample_ht = sprintf("%s_%s", name, HTO_classification)) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$sample_ht, col.name = "sample_ht")      
}

DefaultAssay(s.obj) <- "RNA"

path.to.07.output <- file.path(outdir, 
                               "GEX_output", 
                               sprintf("07_output"), 
                               dataset.name, 
                               clone.name)
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)
if (dataset.name %in% c("241104_BSimons", "241002_BSimons")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

Idents(s.obj) <- "seurat_clusters"