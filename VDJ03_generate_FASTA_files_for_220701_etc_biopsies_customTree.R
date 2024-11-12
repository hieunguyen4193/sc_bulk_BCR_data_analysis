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
    yfp.mids[[mouse.id]] <- list(
      all_samples_including_biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id)$X,
      all = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
      pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$
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

re_define_clone_cluster <- TRUE
rerun <- FALSE
savefile <- TRUE
verbose <- TRUE
save_fasta <- TRUE
define.clone.clusters <- FALSE

clone.output <- run_preprocessing_all_bulk_VDJ_data(
  path.to.mid.output = path.to.mid.output, 
  path.to.save.output = path.to.save.output,
  PROJECT = PROJECT,
  thres = thres, 
  thres.dis = thres.dis,
  savefile = savefile,
  verbose = verbose,
  rerun = rerun,
  define.clone.clusters = define.clone.clusters
)  

clonesets <- clone.output$clonesets

for (mouse.id in names(yfp.mids)){
  for (input.case in c("all", "neg", "pos", "biopsy")){
    if (input.case == "biopsy"){
      save.folder.name <- input.case
    } else {
      save.folder.name <- sprintf("%s_YFP", input.case)
    }
    input.clonesets <- subset(clonesets, clonesets$id %in% yfp.mids[[mouse.id]][[input.case]])
    ##### This script might not run when in RSTUDIO, run in BASH command line. 
    output.all.fasta <- generate_fasta(
       clonesets = input.clonesets, 
       mouse.id = mouse.id,
       path.to.save.output = file.path(path.to.save.fasta, mouse.id, input.case), 
       ref.gene = ref.gene, 
       ref.gene.config = ref.gene.config,
       PROJECT = PROJECT,
       thres = 0.85,
       thres.dis = 0.15,
       save_fasta = save_fasta,
       re_define_clone_cluster = re_define_clone_cluster)
  }
}
