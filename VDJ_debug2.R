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

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####

path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
thres <- 0.85
thres.dis <- 0.15
PROJECT <- "220701_etc_biopsies"
ref.gene <- "IMGT"
ref.gene.config <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/ref_gene_config.R"

#####----------------------------------------------------------------------#####
##### PATHS
#####----------------------------------------------------------------------#####
path.to.save.fasta <- file.path(outdir, "VDJ_output", "03_output", "FASTA_output", PROJECT,  sprintf("VDJ_output_%s", thres))
dir.create(path.to.save.fasta, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### METADATA
#####----------------------------------------------------------------------#####
mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")
count.mid.in.mouse <- table(mid.metadata$population, mid.metadata$mouse) %>% colSums()
if (file.exists(file.path(path.to.save.fasta, "sample_list_based_on_YFP.rds")) == FALSE){
  yfp.mids <- list()
  for (mouse.id in names(count.mid.in.mouse[count.mid.in.mouse >= 4])){
    yfp.mids[[mouse.id]] <- list(all = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
                                 pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$X,
                                 neg = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)$X,
                                 biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$population == "biopsy")$X)   
  }
  saveRDS(yfp.mids, file.path(path.to.save.fasta, "sample_list_based_on_YFP.rds"))
} else {
  yfp.mids <- readRDS(file.path(path.to.save.fasta, "sample_list_based_on_YFP.rds"))
}

#####----------------------------------------------------------------------#####
##### READ THE CLONE DATA
#####----------------------------------------------------------------------#####
path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

clone.output <- run_preprocessing_all_bulk_VDJ_data(
  path.to.mid.output = path.to.mid.output, 
  path.to.save.output = path.to.save.output,
  PROJECT = PROJECT,
  thres = thres, 
  thres.dis = thres.dis,
  savefile = TRUE,
  verbose = TRUE,
  rerun = FALSE,
  define.clone.clusters = FALSE
)  

clonesets <- clone.output$clonesets

input.case <- "all"
mouse.id <- "m11"
clonesets <- subset(clonesets, clonesets$id %in% yfp.mids[[mouse.id]][[input.case]])

new.clonesets <- data.frame()
for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
  tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
  seqs <- unique(tmpdf$aaSeqCDR3)
  print(sprintf("VJ.len.combi: %s, num seqs: %s", input.VJ.combi, length(seqs)))
  if (length(seqs) >= 2){
    cluster.output <- assign_clusters_to_sequences(seqs = seqs, threshold = thres.dis)$res
    tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
      tmpdf$aaSeqCDR3, function(x){
        return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
      }
    ))    
  } else {
    tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
  }
  new.clonesets <- rbind(new.clonesets, tmpdf)
}
clonesets <- new.clonesets
print("IMPORTANT NOTE: THE CLONE CLUSTERS ARE GENERATED BASED ON SAMPLES/MIDS IN THIS SELECTED MOUSE/SAMPLE ONLY; NOT THE AGGREGATED TABLE.")
old.clonesets <- readxl::read_excel("/media/hieunguyen/GSHD_HN01/raw_data/mixcr_gctree_archived_output/mixcr_pipeline_output/data_analysis/05_output/CDR3_0.15/m11/all_YFP_MIDs/clonesets_m11.split_clones_0.15.xlsx")
old.clonesets <- old.clonesets %>% rowwise() %>%
  mutate(VJcombi_CDR3_0.15 = str_replace_all(VJcombi_CDR3_0.15, "[*]01", "")) %>%
  mutate(VJcombi_CDR3_0.15 = str_replace_all(VJcombi_CDR3_0.15, "[*]02", "")) %>%
  mutate(VJcombi_CDR3_0.15 = str_replace_all(VJcombi_CDR3_0.15, "[*]03", "")) %>%
  mutate(VJcombi_CDR3_0.15 = str_replace_all(VJcombi_CDR3_0.15, "[*]04", ""))

# setdiff(old.clonesets$VJcombi_CDR3_0.15, clonesets$VJcombi_CDR3_0.85)
# setdiff(clonesets$VJcombi_CDR3_0.85, old.clonesets$VJcombi_CDR3_0.15)

source(ref.gene.config) # path to the configuration file for the input BCR reference genes
if (ref.gene == "10x"){
  ref.fasta <- readDNAStringSet(ref.genes$`10x`)
} else if (ref.gene == "IMGT"){
  s.V.genes <- readDNAStringSet(ref.genes$IMGT$V.gene)
  names(s.V.genes) <- lapply(names(s.V.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    x <- str_split(x, "[*]")[[1]][[1]]
    return(x)
  })
  s.J.genes <- readDNAStringSet(ref.genes$IMGT$J.gene)
  names(s.J.genes) <- lapply(names(s.J.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    x <- str_split(x, "[*]")[[1]][[1]]
    return(x)
  })
}

save_fasta <- TRUE
path.to.save.output <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1/VDJ_output/test"
all.saved.files <- Sys.glob(file.path(path.to.save.output, "*"))
all.saved.files <- to_vec(
  for (item in all.saved.files){
    str_replace(basename(item), ".fasta", "")
  }
)

path.to.old.files <- "/media/hieunguyen/GSHD_HN01/raw_data/mixcr_gctree_archived_output/mixcr_pipeline_output/data_analysis/05_output/CDR3_0.15/m11/all_YFP_MIDs"
old.saved.files <- Sys.glob(file.path(path.to.old.files, "*.fasta"))
old.saved.files <- to_vec(
  for (item in old.saved.files){
    str_replace(basename(item), ".aln.fasta", "") %>% str_replace("m11_all_YFP_", "") %>% 
      str_replace_all("-01", "") %>% 
      str_replace_all("-02", "") %>% 
      str_replace_all("-03", "") %>% 
      str_replace_all("-04", "")
  }
)

setdiff(clonesets$VJcombi_CDR3_0.85, all.saved.files)

setdiff(old.saved.files, all.saved.files)

setdiff(all.saved.files, clonesets$VJcombi_CDR3_0.85)

for (i in setdiff(old.saved.files, all.saved.files)){
  print(subset(old.clonesets, clonesets$VJcombi_CDR3_0.85 == i) %>% nrow())   
}

for (i in setdiff(old.saved.files, all.saved.files)){
  print(subset(new.clonesets, clonesets$VJcombi_CDR3_0.85 == i) %>% nrow()) 
}

