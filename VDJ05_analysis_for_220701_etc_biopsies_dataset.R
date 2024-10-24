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
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output")
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
thres <- 0.85
clonedf <- read.csv(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"), row.names = "X")
clonedf <- subset(clonedf, clonedf$num_mutation != "region_not_covered-skip")

clonedf$num_mutation <- as.numeric(clonedf$num_mutation)
mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")
count.mid.in.mouse <- table(mid.metadata$population, mid.metadata$mouse) %>% colSums()

if (file.exists(file.path(path.to.05.output, "sample_list_based_on_YFP.rds")) == FALSE){
  yfp.mids <- list()
  for (mouse.id in names(count.mid.in.mouse[count.mid.in.mouse >= 4])){
    yfp.mids[[mouse.id]] <- list(all = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
                                 pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$X,
                                 neg = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)$X,
                                 biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$population == "biopsy")$X)   
  }
  saveRDS(yfp.mids, file.path(path.to.05.output, "sample_list_based_on_YFP.rds"))
} else {
  yfp.mids <- readRDS(file.path(path.to.05.output, "sample_list_based_on_YFP.rds"))
}

mouse.id <- "m53"
tmpdf <- subset(clonedf, clonedf$id %in% yfp.mids[[mouse.id]]$all)

