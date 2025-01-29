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

# circos.group.type <- "VJnt"
circos.group.type <- "VJaa"

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
bulk.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241031_BSimons/241031_sample_sheet.xlsx") %>%
  rowwise() %>% 
  mutate(MID = sprintf("MID%s", MID))
sc.projects <- c("241002_241104_BSimons")
bulk.projects <- c("241031_BSimons")
list.of.PROJECT <- c(sc.projects, bulk.projects)

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, circos.group.type), showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
#### read the Seurat object of single cell dataset
######----------------------------------------------------------------------#####
s.obj <- list()
sc.meta.data <- list()
list.of.ht <- list()
for (PROJECT in sc.projects){
  s.obj[[PROJECT]] <- readRDS(path.to.all.s.obj[[PROJECT]])
  sc.meta.data[[PROJECT]] <- s.obj[[PROJECT]]@meta.data %>% rownames_to_column("barcode")
  list.of.ht[[PROJECT]] <- list()
  for (sample.id in unique(sc.meta.data[[PROJECT]]$name)){
    list.of.ht[[PROJECT]][[sample.id]] <- unique(subset(sc.meta.data[[PROJECT]], sc.meta.data[[PROJECT]]$name == sample.id)$HTO_classification)
  }
}

yfp.sample.list <- list()
for (mouse.id in c("m3", "m7")){
  yfp.sample.list[[mouse.id]] <- list()
  for (i in unique(bulk.metadata$YFP)){
    yfp.sample.list[[mouse.id]][[i]] <- subset(bulk.metadata, bulk.metadata$YFP == i & bulk.metadata$mouse == mouse.id)$MID
  }  
}

injected.samples <- list(
  m3 = "PP3_HT2",
  m7 = "PP7_HT1"
)

#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> NEW DATA
#####----------------------------------------------------------------------#####
print(sprintf("All samples data exists, reading in ..."))
all.data <- readRDS(file.path(path.to.05.output, circos.group.type, "all_data.rds"))
meta.data <- readxl::read_excel(file.path(path.to.05.output, "all_data_metadata.xlsx"))

##### generate circos plot for all hashtags
exclude.samples <- c("PP3", "PP7")

meta.data.splitted.or.not <- list(
  with_hashtags = subset(meta.data, meta.data$SampleID %in% exclude.samples == FALSE),
  without_hashtags = subset(meta.data, grepl("_", meta.data$SampleID) == FALSE)
)

meta.data.name <- "with_hashtags"

  tmp.metadata <- meta.data.splitted.or.not[[meta.data.name]]
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
  
  all.input.files <- input.metadata$path
  names(all.input.files) <- input.metadata$SampleID
  
for (mouse.id in c("m3", "m7")){
  yfp.samples <- yfp.sample.list[[mouse.id]]$`YFP+`
  nonyfp.samples <- yfp.sample.list[[mouse.id]]$`YFP-`
  # for (mouse.id in c("m3", "m7")){
  #   selected.mids <- subset(tmp.metadata, tmp.metadata$mouse == mouse.id)$SampleID
  
  tmp.yfpdf <- read_tsv(all.input.files[yfp.samples])
  tmp.nonyfpdf <- read_tsv(all.input.files[nonyfp.samples])
  
  yfpdf <- tmp.yfpdf %>%
    group_by(id) %>%
    summarise(new.cloneCount = sum(cloneCount))
  yfpdf <- yfpdf %>% rowwise() %>%
    mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
    mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
    mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]])
  colnames(yfpdf) <- c("id", "cloneCount", "bestVHit", "bestJHit", "nSeqCDR3")
  
  nonyfpdf <- tmp.nonyfpdf %>%
    group_by(id) %>%
    summarise(new.cloneCount = sum(cloneCount))
  nonyfpdf <- nonyfpdf %>% rowwise() %>%
    mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
    mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
    mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]])
  colnames(nonyfpdf) <- c("id", "cloneCount", "bestVHit", "bestJHit", "nSeqCDR3")
  
  dir.create(file.path(path.to.05.output, 
                       sprintf("VJcombi_CDR3_%s", thres), 
                       meta.data.name, 
                       "MERGE_YFP_SAMPLES"), showWarnings = FALSE, recursive = TRUE)
  write.table(yfpdf, 
              file.path(path.to.05.output, 
                        sprintf("VJcombi_CDR3_%s", thres), 
                        meta.data.name, 
                        "MERGE_YFP_SAMPLES",
                        sprintf("YFP_%s.simplified.csv", mouse.id)), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE) 
  write.table(nonyfpdf, 
              file.path(path.to.05.output, 
                        sprintf("VJcombi_CDR3_%s", thres), 
                        meta.data.name, 
                        "MERGE_YFP_SAMPLES",
                        sprintf("NONYFP_%s.simplified.csv", mouse.id)), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE) 
  
  new.input.files <- list(
    injected_PP = all.input.files[[injected.samples[[mouse.id]]]],
    YFP_positive = file.path(path.to.05.output, 
                             sprintf("VJcombi_CDR3_%s", thres), 
                             meta.data.name, 
                             "MERGE_YFP_SAMPLES",
                             sprintf("YFP_%s.simplified.csv", mouse.id)),
    YFP_negative = file.path(path.to.05.output, 
                             sprintf("VJcombi_CDR3_%s", thres), 
                             meta.data.name, 
                             "MERGE_YFP_SAMPLES",
                             sprintf("NONYFP_%s.simplified.csv", mouse.id))
  )
  
  fileAliases <- c("injected PP", "YFP+", "YFP-")
  
  names(fileAliases) <- names(new.input.files)
  saveFileName <- sprintf("%s_circos_MERGE_YFP_INJECTED_PP.svg", mouse.id)
  outputdir <- file.path(path.to.05.output,
                         sprintf("VJcombi_CDR3_%s", thres),
                         "circos_plot")
  filter.clone <- FALSE
  filter.clone.cutoff <- NA
  source(file.path(path.to.main.src, "circos_helper.R"))
  
  if (file.exists(file.path(outputdir, saveFileName)) == FALSE){
    generate_circos(
      input.files = new.input.files,
      fileAliases = fileAliases,
      saveFileName = saveFileName,
      outputdir = outputdir,
      filter.clone = filter.clone,
      filter.clone.cutoff = NULL,
      group.to.highlight1 = NULL,
      group.to.highlight2 = NULL,
      linkColor1 = "#FF000080",
      linkColor2 = "#FF000080",
      ordered.samples = NULL
    )
  }
}
