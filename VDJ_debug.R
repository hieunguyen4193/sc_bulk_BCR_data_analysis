gc()
rm(list = ls())

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
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.04.output <- file.path(outdir, "VDJ_output", "04_output")
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

thres <- 0.85

#####---------------------------------------------------------------------------#####
##### GENERATE A BIG DATAFRAME CONTAINING ALL CLONES FROM ALL DATASETS
#####---------------------------------------------------------------------------#####
# install.packages("https://github.com/igraph/rigraph/archive/refs/tags/v2.1.1.tar.gz", type = "source", repos = NULL)

all.clone.files <- Sys.glob(file.path(outdir, "VDJ_output", "*", sprintf("VDJ_output_%s", thres), "preprocessed_files", "clonesets*.split_clones.xlsx" ))

dataset.origin <- list(
  `1st_2nd_BSimons_Datasets` = "sc",
  `220701_etc_biopsies` = "bulk",
  `240805_BSimons` = "sc",
  `240826_BSimons` = "bulk",
  `241002_BSimons` = "sc"
)

names(all.clone.files) <- to_vec(
  for (item in all.clone.files) str_replace(str_replace(basename(item), "clonesets_", ""), ".split_clones.xlsx", "")
)

mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")
count.mid.in.mouse <- table(mid.metadata$population, mid.metadata$mouse) %>% colSums()

#####----------------------------------------------------------------------#####
##### PREPARE LIST OF MID W.R.T EACH CASE, ALL SAMPLES, BIOPSY ONLY OR YPF+-
#####----------------------------------------------------------------------#####
yfp.mids <- list()
for (mouse.id in names(count.mid.in.mouse[count.mid.in.mouse >= 4])){
  yfp.mids[[mouse.id]] <- list(all = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
                               pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$X,
                               neg = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)$X,
                               biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$population == "biopsy")$X)   
}
tmpdf <- readxl::read_excel(all.clone.files[["220701_etc_biopsies"]])
tmpdf <- subset(tmpdf, tmpdf$id %in% yfp.mids[["m42"]]$all)
subset.tmpdf <- subset(tmpdf, tmpdf$VJ.len.combi == "IGHV1-52_IGHJ4_60")
seqs <- unique(subset.tmpdf$aaSeqCDR3)
# cluster.output <- assign_clusters_to_sequences(seqs, threshold = 0.15)$res
threshold <- 0.15

distance_matrix <- compute_distance_matrix(seqs)
max_length <- max(nchar(seqs))
normalized_matrix <- distance_matrix / max_length
graph <- graph_from_adjacency_matrix(normalized_matrix < threshold, mode = "undirected", diag = FALSE)
cl <- cluster_optimal(graph)
clus.res <- data.frame(seq = seqs, cluster = cl$membership)
