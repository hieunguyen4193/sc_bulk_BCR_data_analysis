# gc()
# rm(list = ls())

scrna_pipeline_src <- "/home/uk104163/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))

path.to.main.src <-  "/home/uk104163/src_2023/BSimons/CRC1382_BSimons_project"

main.src.dir <- "/home/uk104163/src_2023/BSimons/BSimons_thesis_data"
source(file.path(main.src.dir, "00_helper_functions.R"))

outdir <- "/media/outdir"

path.to.output <- file.path(outdir, "BSimons", "THESIS_OUTPUT_20231026")
path.to.01.output <- file.path(path.to.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

s.obj.1st.dataset <- readRDS(file.path(path.to.01.output, "raw_seurat_objects", "merged_all_first_dataset_BCR.integrated.rds")) 
s.obj.2nd.dataset <- readRDS(file.path(path.to.01.output, "raw_seurat_objects", "merged_all_second_dataset_BCR.integrated.rds"))

#####----------------------------------------------------------------------#####
##### CONFIGURATIONS
#####----------------------------------------------------------------------#####
# always re integrate the dataset after removing a cluster
with.re.integration <- TRUE

if (with.re.integration == TRUE){
  re.integrate.status <- "with_reInt"
} else {
  re.integrate.status <- "without_reInt"
}

#####----------------------------------------------------------------------#####
##### DATASET 2
#####----------------------------------------------------------------------#####
# Remove cluster 5 and 6 from the 2nd dataset object merged_all_second_dataset_BCR.integrated.rds
# and re integration
removed.clusters <- c(5, 6)
cluster.resolution <- 1

if (file.exists(file.path(path.to.01.output, sprintf("2nd_dataset_removed_%s.%s.res%s", 
                                                     paste(removed.clusters, collapse =  "_"), re.integrate.status, cluster.resolution))) == FALSE){
  remove_clusters_from_an_object(s.obj.2nd.dataset, 
                                 clusters.to.be.removed = removed.clusters, 
                                 with.re.integration = FALSE, 
                                 save.dataset.name = "2nd_dataset",
                                 path.to.output = path.to.01.output, 
                                 cluster.resolution = cluster.resolution,
                                 generate.html = FALSE)
}

#####----------------------------------------------------------------------#####
##### DATASET 1
#####----------------------------------------------------------------------#####
# Remove cluster 7 and 9 from the 1st dataset object merged_all_first_dataset_BCR.integrated.rds
removed.clusters <- c(7, 9)
cluster.resolution <- 1

if (file.exists(file.path(path.to.01.output, sprintf("1st_dataset_removed_%s.%s.res%s", 
                                                     paste(removed.clusters, collapse =  "_"), re.integrate.status, cluster.resolution))) == FALSE){
  remove_clusters_from_an_object(s.obj.1st.dataset, 
                                 clusters.to.be.removed = removed.clusters, 
                                 with.re.integration = FALSE, 
                                 save.dataset.name = "1st_dataset",
                                 path.to.output = path.to.01.output, 
                                 cluster.resolution = cluster.resolution,
                                 generate.html = FALSE)
}


# Remove cluster 16 from the 1st dataset object 1st_dataset_removed_7_9.with_reInt.rds
removed.clusters <- c(16)
cluster.resolution <- 1

if (file.exists(file.path(path.to.01.output, sprintf("1st_dataset_removed_7_9_removed_%s.%s.res%s", 
                                                     paste(removed.clusters, collapse =  "_"), re.integrate.status, cluster.resolution))) == FALSE){
  s.obj.1st.dataset.tmp <- readRDS(file.path(path.to.01.output, 
                                             "1st_dataset_removed_7_9.without_reInt.res1",
                                             "1st_dataset_removed_7_9.without_reInt.res1.rds"))
  
  remove_clusters_from_an_object(s.obj.1st.dataset.tmp, 
                                 clusters.to.be.removed = removed.clusters, 
                                 with.re.integration = FALSE, 
                                 save.dataset.name = "1st_dataset_removed_7_9",
                                 path.to.output = path.to.01.output, 
                                 cluster.resolution = cluster.resolution,
                                 generate.html = FALSE)
}
# Remove cluster 9 from the 1st dataset object 1st_dataset_removed_7_9_removed_16.with_reInt.rds
removed.clusters <- c(9)
cluster.resolution <- 1

if (file.exists(file.path(path.to.01.output, sprintf("1st_dataset_removed_7_9_removed_16_%s.%s.res%s", 
                                                     paste(removed.clusters, collapse =  "_"), re.integrate.status, cluster.resolution))) == FALSE){
  s.obj.1st.dataset.tmp <- readRDS(file.path(path.to.01.output, 
                                             "1st_dataset_removed_7_9_removed_16.without_reInt.res1",
                                             "1st_dataset_removed_7_9_removed_16.without_reInt.res1.rds"))
  
  remove_clusters_from_an_object(s.obj.1st.dataset.tmp,
                                 clusters.to.be.removed = c(9),
                                 with.re.integration = TRUE,
                                 save.dataset.name = "1st_dataset_removed_7_9_removed_16",
                                 path.to.output = path.to.01.output,
                                 cluster.resolution = cluster.resolution,
                                 generate.html = FALSE)
}



tmp.s.obj <- readRDS("/media/outdir/BSimons/THESIS_OUTPUT_20231026/01_output/1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1/1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1.rds")

