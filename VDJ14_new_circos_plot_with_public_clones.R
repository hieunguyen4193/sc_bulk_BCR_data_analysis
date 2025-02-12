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
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

sc.projects <- c("240805_BSimons_filterHT_cluster_renamed")
bulk.projects <- c("240826_BSimons")
list.of.PROJECT <- c(sc.projects, bulk.projects)

thres <- 0.85 

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
path.to.14.output <- file.path(outdir, "VDJ_output", "14_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.14.output, showWarnings = FALSE, recursive = TRUE)

path.to.input.circos.plots <- file.path(path.to.14.output, "input")
dir.create(path.to.input.circos.plots, showWarnings = FALSE, recursive = TRUE)

meta.data.name <- "without_hashtags"

all.input.files <- Sys.glob(file.path(path.to.05.output, 
                                      sprintf("VJcombi_CDR3_%s", thres), 
                                      meta.data.name,
                                      "*.simplified.csv"))

input.metadata <- data.frame(
  path = all.input.files,
  SampleID = to_vec(for (item in all.input.files){
    str_replace(basename(item), ".simplified.csv", "") 
  }),
  PROJECT = to_vec(for (item in all.input.files){
    str_split(item, "/")[[1]][[8]]
  })
) 

public.clonedf <- read.csv(file.path(path.to.main.src, "public_clones.csv")) %>%
  rowwise() %>% 
  mutate(VJ = sprintf("%s_%s", V_gene, J_gene))

df <- read.csv(all.input.files[[1]], sep = "\t") %>%
  rowwise() %>%
  mutate(VJ.gene = sprintf("%s_%s", bestVHit, bestJHit))
