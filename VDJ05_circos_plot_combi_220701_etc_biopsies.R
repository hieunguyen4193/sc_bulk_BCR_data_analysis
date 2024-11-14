gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))
source(file.path(path.to.main.src, "VDJ_generate_circos_plots.R"))
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
library(circlize)
#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")
yfp.mids <- list()
for (mouse.id in unique(mid.metadata$mouse)){
  yfp.mids[[mouse.id]] <- list(
    all_w_biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id)$X,
    all_yfp = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
    pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$X,
    neg = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)$X,
    biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$population == "biopsy")$X)   
}

PROJECT <- "220701_etc_biopsies"
path.to.VDJ.output <- file.path( outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

#####----------------------------------------------------------------------#####
##### READ CLONE DATA
#####----------------------------------------------------------------------#####
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", PROJECT)
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, "circos_plot"), showWarnings = FALSE, recursive = TRUE)

path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

clone.obj <- run_preprocessing_all_bulk_VDJ_data(path.to.mid.output = path.to.mid.output,
                                                 path.to.save.output = path.to.save.output,
                                                 PROJECT = PROJECT,
                                                 thres = thres, 
                                                 thres.dis = thres.dis,
                                                 savefile = savefile,
                                                 verbose = verbose,
                                                 rerun = rerun,
                                                 define.clone.clusters = define.clone.clusters)   

count.mice <- table(mid.metadata$mouse)
all.mice <- count.mice[count.mice >= 2] %>% names()

#####----------------------------------------------------------------------#####
##### generate clonesets and assign clones to clusters for each set 
#####----------------------------------------------------------------------#####
# for (mouse.id in all.mice){
#   for (input.case in setdiff(names(yfp.mids[[mouse.id]]), c("biopsy"))  ){
#     selected.mids <- yfp.mids[[mouse.id]][[input.case]]
#     
#     clonesets <- subset(clone.obj$clonesets, clone.obj$clonesets$id %in% selected.mids)
#     if (file.exists(file.path(path.to.05.output, sprintf("%s_%s_%s.clonesets.csv", PROJECT, mouse.id, input.case))) == FALSE){
#       new.clonesets <- data.frame()
#       for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
#         tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
#         seqs <- unique(tmpdf$aaSeqCDR3)
#         print(sprintf("VJ.len.combi: %s, num seqs: %s", input.VJ.combi, length(seqs)))
#         if (length(seqs) >= 2){
#           cluster.output <- assign_clusters_to_sequences(seqs = seqs, threshold = thres.dis)$res
#           tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
#             tmpdf$aaSeqCDR3, function(x){
#               return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
#             }
#           ))    
#         } else {
#           tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
#         }
#         new.clonesets <- rbind(new.clonesets, tmpdf)
#       }
#       
#       clonedf <- table(new.clonesets$VJcombi_CDR3_0.85) %>% data.frame()
#       colnames(clonedf) <- c("clone", "count")
#       clonedf <- clonedf %>% arrange(desc(count))
#       
#       #####----------------------------------------------------------------------#####
#       ##### PLOT CIRCOS 
#       #####----------------------------------------------------------------------#####
#       write.csv(new.clonesets, file.path(path.to.05.output, sprintf("%s_%s_%s.clonesets.csv", PROJECT, mouse.id, input.case)))
#     } else {
#       new.clonesets <- read.csv(file.path(path.to.05.output, sprintf("%s_%s_%s.clonesets.csv", PROJECT, mouse.id, input.case)))
#     }
#   }
# }

#####----------------------------------------------------------------------#####
##### MAIN RUN
#####----------------------------------------------------------------------#####
# mouse.id <- "m12"
# input.case <- "all_w_biopsy"
# clone.def <- "VJnt"

for (mouse.id in all.mice){
  for (input.case in c("all_w_biopsy",
                       "all_yfp",
                       "neg",
                       "pos")){
    for (clone.def in c("VJnt",
                        "VJaa",
                        sprintf("VJcombi_CDR3_%s", thres))){
      selected.mids <- yfp.mids[[mouse.id]][[input.case]]
      new.clonesets <- read.csv(file.path(path.to.05.output, sprintf("%s_%s_%s.clonesets.csv", PROJECT, mouse.id, input.case))) %>%
        subset(aaSeqCDR3 != "region_not_covered") %>%
        rowwise() %>%
        mutate(VJnt = sprintf("%s_%s", VJ.combi, nSeqCDR3)) %>%
        mutate(VJaa = sprintf("%s_%s", VJ.combi, aaSeqCDR3))
      
      populations <- to_vec(
        for (item in subset(mid.metadata, mid.metadata$X %in% selected.mids)$population %>% unique()){
          sprintf("%s_%s", mouse.id, item)
        }
      )
      names(populations) <- unlist(lapply(
        populations, function(x){
          subset(mid.metadata, mid.metadata$population == str_split(x, "_")[[1]][[2]] & mid.metadata$mouse == mouse.id)$X
        }
      ))
      
      input.clonesets <- new.clonesets
      path.to.save.svg <- file.path(path.to.05.output, "circos_plot", mouse.id, input.case, clone.def)
      svg.name <- sprintf("%s_mouse_%s_input_%s.%s.svg", PROJECT, mouse.id, input.case, clone.def)
      clone.def <- clone.def
      data.type <- "bulk"
      populations <- populations
      if (file.exists(file.path(path.to.save.svg, svg.name)) == FALSE){
        generate_circos_plot(input.clonesets = new.clonesets,
                             path.to.save.svg = path.to.save.svg,
                             svg.name = svg.name,
                             clone.def = clone.def,
                             data.type = data.type,
                             populations = populations)        
      } else {
        print(sprintf("File %s exists", file.path(path.to.save.svg, svg.name)))
      }
    }
  }
}




