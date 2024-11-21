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
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

# circos.group.type <- "VJnt"
circos.group.type <- "VJaa"

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
bulk.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240826_BSimons/240829 sample sheet.xlsx")
sc.projects <- c("240805_BSimons")
bulk.projects <- c("240826_BSimons")
list.of.PROJECT <- c(sc.projects, bulk.projects)

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, circos.group.type), showWarnings = FALSE, recursive = TRUE)

meta.data <- readxl::read_excel(file.path(path.to.05.output, "all_data_metadata.xlsx"))
exclude.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
meta.data.splitted <- subset(meta.data, meta.data$SampleID %in% exclude.samples == FALSE)
meta.data.non.splitted <- subset(meta.data, grepl("_", meta.data$SampleID) == FALSE)

all.input.files <- Sys.glob(file.path(outdir, "VDJ_output",
                                      "*",
                                      sprintf("VDJ_output_%s", thres),
                                      "preprocessed_files",
                                      circos.group.type,
                                      "*"))

input.metadata <- data.frame(
  path = all.input.files,
  SampleID = to_vec(for (item in all.input.files){
    str_replace(basename(item), ".simplified.csv", "") 
  }),
  PROJECT = to_vec(for (item in all.input.files){
    str_split(item, "/")[[1]][[8]]
  })
) %>%
  subset(PROJECT %in% list.of.PROJECT)

all.input.files <- input.metadata$path
names(all.input.files) <- input.metadata$SampleID

mouse.id <- "m1"
selected.mids <- subset(meta.data.non.splitted, meta.data.non.splitted$mouse == mouse.id)$SampleID
input.files <- all.input.files[selected.mids]

fileAliases <- to_vec(
  for (item in names(input.files)){
    sprintf("%s (%s)", item, subset(meta.data.non.splitted, meta.data.non.splitted$SampleID == item)$organ)
  }
)
saveFileName <- sprintf("%s_circos.svg", mouse.id)
outputdir <- file.path(path.to.05.output,
                       circos.group.type,
                       "circos_plot")
filter.clone <- FALSE
filter.clone.cutoff <- NA
source(file.path(path.to.main.src, "circos_helper.R"))

if (file.exists(file.path(outputdir, saveFileName)) == FALSE){
  generate_circos(
    input = input.files,
    fileAliases = fileAliases,
    saveFileName = saveFileName,
    outputdir = outputdir,
    filter.clone = filter.clone,
    filter.clone.cutoff = filter.clone.cutoff
  )
}

cloneCountdf <- data.frame()
for (input.mid in names(input.files)){
  print(sprintf("reading in clones from sample %s", input.mid))
  tmp.clonedf <- read_tsv(input.files[[input.mid]])
  print(dim(tmp.clonedf))
  if (filter.clone == TRUE){
    tmp.clonedf <- subset(tmp.clonedf, tmp.clonedf$cloneCount > filter.clone.cutoff)
  }
  totalCloneCount <- tmp.clonedf$cloneCount %>% sum()
  tmp.clonedf <- tmp.clonedf %>%
    rowwise() %>%
    mutate(Freq = cloneCount/totalCloneCount) %>%
    arrange(Freq)
  accum.sum <- 0
  all.accum.sum <- list()
  for (i in seq(1, nrow(tmp.clonedf))){
    accum.sum <- tmp.clonedf[i, ][["Freq"]] + accum.sum
    all.accum.sum <- c(all.accum.sum, c(accum.sum))
  }
  tmp.clonedf$accum.Freq <- unlist(all.accum.sum)
  tmp.clonedf <- subset(tmp.clonedf, select = c(id, Freq, accum.Freq))
  tmp.clonedf$SampleID <- input.mid
  cloneCountdf <- rbind(cloneCountdf, tmp.clonedf)
}
