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
PROJECT <- "240826_BSimons"
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
mid.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240826_BSimons/240829 sample sheet.xlsx")

sample.list <- list()
for (mouse.id in unique(mid.metadata$mouse)){
  sample.list[[mouse.id]] <- list(
    all = subset(mid.metadata, mid.metadata$mouse == mouse.id)$MID,
    without_colon_sample = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$organ != "colon")$MID
  )
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

for (mouse.id in names(sample.list)){
  for (input.case in c("all", "without_colon_sample")){
    input.clonesets <- subset(clonesets, clonesets$id %in% sample.list[[mouse.id]][[input.case]])
    ##### This script might not run when in RSTUDIO, run in BASH command line. 
    output.all.fasta <- generate_fasta(clonesets = input.clonesets, 
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
