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

mouse.id <- "m11"

path.to.input <- file.path(path.to.05.output, "input", mouse.id)
dir.create(path.to.input, showWarnings = FALSE, recursive = TRUE)

path.to.clonesets <- file.path(outdir, 
                               "VDJ_output/03_output/FASTA_output/220701_etc_biopsies/VDJ_output_0.85", 
                               mouse.id, 
                               "all_samples_including_biopsy",
                               "clonesets.csv")

clonedf <- read.csv(path.to.clonesets) %>% 
  rowwise() %>%
  mutate(tmp.V.gene = str_split(str_split(VJcombi_CDR3_0.85, "_")[[1]][[1]], "[*]")[[1]][[1]] ) %>%
  mutate(tmp.J.gene = str_split(str_split(VJcombi_CDR3_0.85, "_")[[1]][[2]], "[*]")[[1]][[1]]) %>%
  mutate(tmp.len = str_split(VJcombi_CDR3_0.85, "_")[[1]][[3]]) %>%
  mutate(tmp.i = str_split(VJcombi_CDR3_0.85, "_")[[1]][[4]]) %>%
  mutate(VJcombi_CDR3_0.85 = sprintf("%s_%s_%s_%s", tmp.V.gene, tmp.J.gene, tmp.len, tmp.i))

all.samples <- unique(clonedf$id)

for (input.MID in all.samples){
  clonedf <- subset(clonedf, clonedf$id == input.MID) 
  clonedf <- subset(clonedf, select = c(VJcombi_CDR3_0.85, 
                                        uniqueMoleculeCount, 
                                        tmp.V.gene, 
                                        tmp.J.gene, 
                                        tmp.len)) 
  new.clonedf <- clonedf  %>%
    group_by(VJcombi_CDR3_0.85) %>%
    summarise(cloneCount = sum(uniqueMoleculeCount)) 
  colnames(new.clonedf) <- c("id", "cloneCount")
  new.clonedf <- new.clonedf %>% 
    rowwise() %>%
    mutate(bestVHit = unique(subset(clonedf, clonedf$VJcombi_CDR3_0.85 == id)$tmp.V.gene)[[1]]) %>%
    mutate(bestJHit = unique(subset(clonedf, clonedf$VJcombi_CDR3_0.85 == id)$tmp.J.gene)[[1]]) %>%
    mutate(nSeqCDR3 = unique(subset(clonedf, clonedf$VJcombi_CDR3_0.85 == id)$tmp.len)[[1]])
  
  write.table(new.clonedf, 
              file.path(path.to.input, 
                        sprintf("%s.simplified.csv", input.MID)), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE) 
}


