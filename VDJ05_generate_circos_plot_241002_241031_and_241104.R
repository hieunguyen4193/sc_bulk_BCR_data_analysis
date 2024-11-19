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

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
bulk.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240826_BSimons/240829 sample sheet.xlsx")
list.of.PROJECT <- c("241002_BSimons", "241031_BSimons", "241104_BSimons")
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### read the Seurat object of single cell dataset
#####----------------------------------------------------------------------#####
s.obj <- list()
s.obj.metadata <- list()
list.of.ht <- list()
for (input.s.obj.name in c("241002_BSimons", "241104_BSimons")){
  s.obj[[input.s.obj.name]] <- readRDS(path.to.all.s.obj[[input.s.obj.name]])
  s.obj.metadata[[input.s.obj.name]] <- s.obj[[input.s.obj.name]]@meta.data %>% rownames_to_column("barcode")
  list.of.ht[[input.s.obj.name]] <- list()
  for (sample.id in unique(s.obj.metadata[[input.s.obj.name]]$name)){
    list.of.ht[[input.s.obj.name]][[sample.id]] <- unique(subset(s.obj.metadata[[input.s.obj.name]], 
                                                                 s.obj.metadata[[input.s.obj.name]]$name == sample.id)$HTO_classification)
  }
}
#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> NEW DATA
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.05.output, "all_data.rds")) == FALSE){
  all.data <- list()
  for (PROJECT in list.of.PROJECT){
    path.to.VDJ.output <- file.path( outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
    path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
    
    path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
    path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
    
    if (PROJECT == "241031_BSimons"){
      clone.obj <- run_preprocessing_all_bulk_VDJ_data(path.to.mid.output = path.to.mid.output,
                                                       path.to.save.output = path.to.save.output,
                                                       PROJECT = PROJECT,
                                                       thres = thres, 
                                                       thres.dis = thres.dis,
                                                       savefile = savefile,
                                                       verbose = verbose,
                                                       rerun = rerun,
                                                       define.clone.clusters = define.clone.clusters) 
    } else if (PROJECT %in% c("241002_BSimons", "241104_BSimons")){
      clone.obj <- run_preprocessing_all_sc_data( path.to.VDJ.output = path.to.VDJ.output, 
                                                  path.to.save.output = path.to.save.output, 
                                                  PROJECT = PROJECT,
                                                  thres = thres, 
                                                  thres.dis = thres.dis,
                                                  savefile = savefile,
                                                  rerun = rerun,
                                                  define.clone.clusters =  define.clone.clusters)
    }
    
    full.clonedf <- clone.obj$clonesets
    full.clonedf <- subset(full.clonedf, full.clonedf$aaSeqCDR3 != "region_not_covered")
    if (PROJECT %in% c("241002_BSimons", "241104_BSimons")){
      full.clonedf <- full.clonedf %>% rowwise() %>%
        mutate(barcode.full = sprintf("%s_%s", id, barcode)) %>%
        # keep only cells that have hashtag information and remain in the 
        # filtered seurat object data.
        mutate(HT = ifelse(barcode.full %in% s.obj.metadata[[PROJECT]]$barcode, 
                           subset(s.obj.metadata[[PROJECT]], 
                                  s.obj.metadata[[PROJECT]]$barcode == barcode.full)$HTO_classification,
                           NA)) %>%
        subset(is.na(HT) == FALSE) %>%
        mutate(id.origin = id) %>%
        mutate(id = sprintf("%s_%s", id, HT))
    }
    
    dir.create(file.path(path.to.save.output, "single_MID_clone_df"), showWarnings = FALSE, recursive = TRUE)
    
    ##### split the full clone dataframe to smaller dataframe for each MID/each single cell sample hashtag
    for (mid in unique(full.clonedf$id)){
      if (file.exists(file.path(path.to.save.output, 
                                "single_MID_clone_df", 
                                sprintf("%s.simplified.csv", mid))) == FALSE){
        
        print(sprintf("Working on sample MID %s", mid))
        clonedf <- subset(full.clonedf, full.clonedf$id == mid)
        input.circos <- subset(clonedf, select = c(V.gene, 
                                                   J.gene, 
                                                   aaSeqCDR3, 
                                                   nSeqCDR3, 
                                                   uniqueMoleculeCount)) %>%
          rowwise() %>%
          mutate(VJnt = sprintf("%s_%s_%s", V.gene, J.gene, nSeqCDR3))
        input.circos <- data.frame(table(input.circos$VJnt))
        colnames(input.circos) <- c("id", "cloneCount")
        input.circos <- input.circos %>% rowwise() %>%
          mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
          mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
          mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]]) %>%
          arrange(desc(cloneCount))
        write.table(input.circos, 
                    file.path(path.to.save.output, 
                              "single_MID_clone_df", 
                              sprintf("%s.simplified.csv", mid)), 
                    quote = FALSE, 
                    sep = "\t", 
                    row.names = FALSE) 
        all.data[[sprintf("%s_%s", PROJECT, mid)]] <- input.circos
      } else {
        print(sprintf("File exists at %s, reading in ...", 
                      file.path(path.to.save.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid))))
        all.data[[sprintf("%s_%s", PROJECT, mid)]] <- read.csv(file.path(path.to.save.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid)), sep = "\t")
      }
    }
    
    ##### if the dataset is a single cell dataset, get clonedf for all hashtag in a mouse sample. 
    if (PROJECT %in% c("241002_BSimons", "241104_BSimons")){
      for (mid in unique(full.clonedf$id.origin)){
        if (file.exists(file.path(path.to.save.output, 
                                  "single_MID_clone_df", 
                                  sprintf("%s.simplified.csv", mid))) == FALSE){
          
          print(sprintf("Working on sample MID %s", mid))
          clonedf <- subset(full.clonedf, full.clonedf$id.origin == mid)
          input.circos <- subset(clonedf, select = c(V.gene, 
                                                     J.gene, 
                                                     aaSeqCDR3, 
                                                     nSeqCDR3, 
                                                     uniqueMoleculeCount)) %>%
            rowwise() %>%
            mutate(VJnt = sprintf("%s_%s_%s", V.gene, J.gene, nSeqCDR3))
          input.circos <- data.frame(table(input.circos$VJnt))
          colnames(input.circos) <- c("id", "cloneCount")
          input.circos <- input.circos %>% rowwise() %>%
            mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
            mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
            mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]]) %>%
            arrange(desc(cloneCount))
          write.table(input.circos, 
                      file.path(path.to.save.output, 
                                "single_MID_clone_df", 
                                sprintf("%s.simplified.csv", mid)), 
                      quote = FALSE, 
                      sep = "\t", 
                      row.names = FALSE) 
          all.data[[sprintf("%s_%s", PROJECT, mid)]] <- input.circos
        } else {
          print(sprintf("File exists at %s, reading in ...", 
                        file.path(path.to.save.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid))))
          all.data[[sprintf("%s_%s", PROJECT, mid)]] <- read.csv(file.path(path.to.save.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid)), sep = "\t")
        }
      }
    }
  }
  saveRDS(all.data, file.path(path.to.05.output, "all_data.rds"))
} else {
  print(sprintf("All samples data exists, reading in ..."))
  all.data <- readRDS(file.path(path.to.05.output, "all_data.rds"))
}
