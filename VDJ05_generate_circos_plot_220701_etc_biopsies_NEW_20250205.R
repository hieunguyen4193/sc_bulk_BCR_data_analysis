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

mid.metadata <- mid.metadata %>% rowwise() %>%
  mutate(sample_type = str_replace(str_replace(population, "Ly6c[+]", ""), "Ly6c[-]", ""))

count.mice <- table(mid.metadata$mouse)
plot.mice <- count.mice[count.mice >= 2] %>% names()
for (mouse.id in plot.mice){
  path.to.input <- file.path(path.to.05.output, "input", mouse.id)
  dir.create(path.to.input, showWarnings = FALSE, recursive = TRUE)
  
  path.to.pool.samplelist <- file.path(path.to.05.output, "pooled_sample_list", mouse.id)
  dir.create(path.to.pool.samplelist, showWarnings = FALSE, recursive = TRUE)
  
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
    if (file.exists(file.path(path.to.input, 
                              sprintf("%s.simplified.csv", input.MID))) == FALSE){
      new.clonedf <- clonedf  %>% 
        subset(id == input.MID) %>%
        group_by(VJcombi_CDR3_0.85) %>%
        summarise(cloneCount = sum(uniqueMoleculeCount)) 
      colnames(new.clonedf) <- c("clone", "cloneCount")
      new.clonedf <- new.clonedf %>% 
        rowwise() %>%
        mutate(bestVHit = unique(subset(clonedf, clonedf$VJcombi_CDR3_0.85 == clone)$tmp.V.gene)[[1]]) %>%
        mutate(bestJHit = unique(subset(clonedf, clonedf$VJcombi_CDR3_0.85 == clone)$tmp.J.gene)[[1]]) %>%
        mutate(nSeqCDR3 = unique(subset(clonedf, clonedf$VJcombi_CDR3_0.85 == clone)$tmp.len)[[1]])
      colnames(new.clonedf) <- c("id", "cloneCount", "bestVHit", "bestJHit", "nSeqCDR3")
      
      write.table(new.clonedf, 
                  file.path(path.to.input, 
                            sprintf("%s.simplified.csv", input.MID)), 
                  quote = FALSE, 
                  sep = "\t", 
                  row.names = FALSE) 
    } else {
      print(sprintf("File %s exists, reading in ...", file.path(path.to.input, 
                                                                sprintf("%s.simplified.csv", input.MID))))
    }
    all.input.files <- Sys.glob(file.path(path.to.input, "*.simplified.csv"))
    
    input.metadata <- data.frame(
      path = all.input.files,
      SampleID = to_vec(for (item in all.input.files){
        str_replace(basename(item), ".simplified.csv", "") 
      }))
    
    input.metadata <- input.metadata %>% rowwise() %>%
      mutate(sample.type = subset(mid.metadata, mid.metadata$X == SampleID)$sample_type)
    
    all.input.files <- input.metadata$path
    names(all.input.files) <- input.metadata$SampleID
    
    sample.list <- list()
    for (input.sample.type in unique(input.metadata$sample.type)){
      sample.list[[input.sample.type]] <- subset(input.metadata, input.metadata$sample.type == input.sample.type)$SampleID
      tmpdf <- read_tsv(all.input.files[sample.list[[input.sample.type]] ])
      df <- tmpdf %>%
        group_by(id) %>%
        summarise(new.cloneCount = sum(cloneCount))
      df <- df %>% rowwise() %>%
        mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
        mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
        mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]])
      colnames(df) <- c("id", "cloneCount", "bestVHit", "bestJHit", "nSeqCDR3")
      write.table(df, 
                  file.path(path.to.pool.samplelist, sprintf("%s_%s.csv", mouse.id, input.sample.type)), 
                  quote = FALSE, 
                  sep = "\t", 
                  row.names = FALSE) 
    }
  }
}

#####----------------------------------------------------------------------------#####
##### RUN CIRCOS PLOT
#####----------------------------------------------------------------------------#####



