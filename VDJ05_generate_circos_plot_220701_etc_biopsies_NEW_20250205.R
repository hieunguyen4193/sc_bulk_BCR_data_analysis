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

all.input.files <- Sys.glob(file.path(path.to.05.output, "single_MID_clone_df", "*.simplified.csv"))
names(all.input.files) <- to_vec(
  for (item in all.input.files){
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

count.clonedf <- data.frame()
for (mouse.id in plot.mice){
  print(sprintf("Working on mouse ID: %s", mouse.id))
  selected.mids <- subset(mid.metadata, mid.metadata$mouse == mouse.id)$X
  input.files <- all.input.files[selected.mids]
  
  filter.clone <- FALSE
  filter.clone.cutoff <- NA
  source(file.path(path.to.main.src, "circos_helper.R"))
  thres.dis <- 0.15
  thres <- 0.85
  
  clonesets <- read_tsv(input.files, id = "fileName") %>%
    rowwise() %>%
    mutate(fileName = basename(fileName) %>% str_replace(".simplified.csv", "")) %>%
    mutate(bestVHit = str_split(bestVHit, "[*]")[[1]][[1]]) %>%
    mutate(bestJHit = str_split(bestJHit, "[*]")[[1]][[1]]) %>%
    mutate(len = nchar(nSeqCDR3)) %>%
    mutate(VJ.len.combi = sprintf("%s_%s_%s", bestVHit, bestJHit, len ))
  
  colnames(clonesets) <- c("fileName", "id", "cloneCount", "bestVHit", "bestJHit", "seq", "len", "VJ.len.combi")
  ##### Group sequences + Gene usages to clones
  new.clonesets <- data.frame()
  for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
    tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
    seqs <- unique(tmpdf$seq)
    print(sprintf("VJ.len.combi: %s, num seqs: %s", input.VJ.combi, length(seqs)))
    if (length(seqs) >= 2){
      cluster.output <- assign_clusters_to_sequences(seqs = seqs, threshold = thres.dis)$res
      tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
        tmpdf$seq, function(x){
          return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
        }
      ))    
    } else {
      tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
    }
    new.clonesets <- rbind(new.clonesets, tmpdf)
  }
  
  tmp.count.clonedf <- data.frame(mouse = c(mouse.id), 
                                  numClone = c(length(unique(new.clonesets$VJcombi_CDR3_0.85))),
                                  totalAbundance = c(sum(new.clonesets$cloneCount)))
  count.clonedf <- rbind(count.clonedf, tmp.count.clonedf)
}

