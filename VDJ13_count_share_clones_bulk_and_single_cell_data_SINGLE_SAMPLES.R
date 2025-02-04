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

thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

circos.group.type <- "VJaa"

source(file.path(path.to.main.src, "convert_sampleID_to_locationName.R"))

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
bulk.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240826_BSimons/240829 sample sheet.xlsx")
# sc.projects <- c("240805_BSimons_filterHT_cluster_renamed")
# bulk.projects <- c("240826_BSimons")

sc.projects <- c("241002_241104_BSimons")
bulk.projects <- c("241031_BSimons")

list.of.PROJECT <- c(sc.projects, bulk.projects)

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

outputdir <- file.path(path.to.05.output,
                       sprintf("VJcombi_CDR3_%s", thres),
                       "circos_plot")

path.to.13.output <- file.path(outdir, "VDJ_output", "13_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.13.output, showWarnings = FALSE, recursive = TRUE)

all.files <- Sys.glob(file.path(outputdir, "*", "m*_circos.csv"))
all.files <- to_vec(for (item in all.files) if (grepl("hashtag", basename(item)) == TRUE) item)

input.file <- all.files[[1]]
# for (input.file in all.files){
  mouse.id <- str_split(basename(input.file), "_")[[1]][[1]]
  countdf <- read.csv(input.file)
  
  all.samples <- to_vec(
    for (item in colnames(countdf)){
      if (grepl("M", item) == TRUE | grepl("P", item) == TRUE | grepl("MID", item) == TRUE){
        str_replace(str_replace(str_replace(item, "_Count", ""), "_Freq", ""), "_accumFreq", "")
      }
    }
  ) %>% unique()
# }