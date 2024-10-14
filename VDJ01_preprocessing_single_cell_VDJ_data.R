gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/00_helper_functions.R")
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
new.pkgs <- c("svglite")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }   
}
new.pkgs <- c("msa")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() ==  FALSE){
    BiocManager::install(pkg, update = FALSE)
  }
}
library("msa")

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis"
PROJECT <- "1st_2nd_BSimons_Datasets"
thres <- 0.85

path.to.VDJ.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres))
path.to.save.output <- file.path(path.to.VDJ.output, "preprocessed_files")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

all.VDJ.files <- Sys.glob(file.path(path.to.save.output, "*.xlsx"))
names(all.VDJ.files) <- unlist(lapply(all.VDJ.files, function(x){
  str_replace(basename(x), ".xlsx", "")
}))
all.VDJ.files <- all.VDJ.files[names(all.VDJ.files) != "GF_7w_14w"]

#####----------------------------------------------------------------------#####
##### pre-processing steps for the input VDJ files
#####----------------------------------------------------------------------#####
all.clonedf <- data.frame()
for (file in all.VDJ.files){
  tmpdf <- readxl::read_xlsx(file)  
  all.clonedf <- rbind(all.clonedf, tmpdf)
}

# Keep IGH chain only, to compare with the BULK data
all.clonedf <- subset(all.clonedf, all.clonedf$chain == "IGH") %>%
  rowwise() %>%
  mutate(VJseq.combi = sprintf("%s_%s_%s_%s", v_gene, j_gene, cdr3, cdr3_nt))
  
new.all.clonedf <- data.frame(
  VJseq.combi.tmp = unique(all.clonedf$VJseq.combi)) %>%
  rowwise() %>%
  mutate(CDR3aa = str_split(VJseq.combi.tmp, "_")[[1]][[3]]) %>%
  mutate(V.gene = str_split(VJseq.combi.tmp, "_")[[1]][[1]]) %>%
  mutate(J.gene = str_split(VJseq.combi.tmp, "_")[[1]][[2]]) %>%
  mutate(CDR3nt = str_split(VJseq.combi.tmp, "_")[[1]][[4]]) %>%
  mutate(cloneSize = nrow(subset(all.clonedf, all.clonedf$VJseq.combi == VJseq.combi.tmp))) %>% 
  mutate(CDR3aa.length = nchar(CDR3aa)) %>% 
  mutate(CDR3nt.length = nchar(CDR3nt)) %>% 
  mutate(samples = paste(unique(subset(all.clonedf, all.clonedf$VJseq.combi == VJseq.combi.tmp)$id), collapse = ",")
)

