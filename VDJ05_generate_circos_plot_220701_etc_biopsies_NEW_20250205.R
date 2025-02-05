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
library(circlize)
#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")
yfp.mids <- list()
sample.list <- list()
for (mouse.id in unique(mid.metadata$mouse)){
  sample.list[[mouse.id]] <- subset(mid.metadata, mid.metadata$mouse == mouse.id)$X
  yfp.mids[[mouse.id]] <- list(
    all_w_biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id)$X,
    all_yfp = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
    pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$X,
    neg = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)$X,
    biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$population == "biopsy")$X)   
}

PROJECT <- "220701_etc_biopsies"
path.to.VDJ.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> NEW DATA
#####----------------------------------------------------------------------#####
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output_NEW_20250205", PROJECT)
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, "circos_plot"), showWarnings = FALSE, recursive = TRUE)

path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
dir.create(file.path(path.to.05.output, "single_MID_clone_df"), showWarnings = FALSE, recursive = TRUE)

clone.obj <- run_preprocessing_all_bulk_VDJ_data(path.to.mid.output = path.to.mid.output,
                                                 path.to.save.output = path.to.VDJ.output,
                                                 PROJECT = PROJECT,
                                                 thres = thres, 
                                                 thres.dis = thres.dis,
                                                 savefile = savefile,
                                                 verbose = verbose,
                                                 rerun = rerun,
                                                 define.clone.clusters = define.clone.clusters)   


full.clonedf <- clone.obj$clonesets
full.clonedf <- subset(full.clonedf, full.clonedf$aaSeqCDR3 != "region_not_covered")

for (mid in unique(mid.metadata$X)){
  if (file.exists(file.path(path.to.05.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid))) == FALSE){
    print(sprintf("Working on sample MID %s", mid))
    clonedf <- subset(full.clonedf, full.clonedf$id == mid)
    input.circos <- subset(clonedf, select = c(V.gene, 
                                               J.gene, 
                                               aaSeqCDR3, 
                                               nSeqCDR3, 
                                               uniqueMoleculeCount)) %>%
      rowwise() %>%
      mutate(VJnt = sprintf("%s_%s_%s", V.gene, J.gene, nSeqCDR3)) %>%
      group_by(VJnt) %>%
      summarise(cloneCount = sum(uniqueMoleculeCount)) %>% 
      arrange(desc(cloneCount))
    colnames(input.circos) <- c("id", "cloneCount")
    input.circos <- input.circos %>% rowwise() %>%
      mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
      mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
      mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]])
    write.table(input.circos, file.path(path.to.05.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid)), quote = FALSE, sep = "\t", row.names = FALSE)  
  }
}

all.mid.files <- Sys.glob(file.path(path.to.05.output, "single_MID_clone_df", "*.simplified.csv"))
names(all.mid.files) <- to_vec(
  for (item in all.mid.files){
    str_replace(basename(item), ".simplified.csv", "")
  }
)

#####----------------------------------------------------------------------#####
##### GENERATE CIRCOS PLOT WITH FABIO SCRIPT
#####----------------------------------------------------------------------#####
source(file.path(path.to.main.src, "circos.R"))
count.mice <- table(mid.metadata$mouse)

#####----------------------------------------------------------------------#####
##### GENERATE CIRCOS PLOT WITH HIEU SCRIPT
#####----------------------------------------------------------------------#####
plot.mice <- count.mice[count.mice >= 2] %>% names()