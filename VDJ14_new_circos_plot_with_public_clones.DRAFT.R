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

path.to.tree.05.output <- file.path(outdir, "GEX_output", "05_output", sc.projects[[1]])
path.to.tree.01.output <- file.path(outdir, "tree_analysis", bulk.projects[[1]], "01_output")

match.bulkdf <- read.csv(file.path(path.to.tree.01.output, "match_public_clonedf.csv"))
match.scdf <- readxl::read_excel(file.path(path.to.tree.05.output, "240805_BSimons_filterHT_cluster_renamed_public_clone_full.xlsx")) %>%
  subset(select = -c(`...1`)) %>% 
  subset(select = c(barcode, VJcombi_CDR3_0.85, dist_to_public_clone))

path.to.input.circos.plots <- file.path(path.to.14.output, "input")
dir.create(path.to.input.circos.plots, showWarnings = FALSE, recursive = TRUE)

meta.data.name <- "without_hashtags"

all.input.files <- Sys.glob(file.path(path.to.05.output, 
                                      sprintf("VJcombi_CDR3_%s", thres), 
                                      meta.data.name,
                                      "*.simplified.csv"))
names(all.input.files) <- to_vec(
  for (item in all.input.files){
    str_replace(basename(item), ".simplified.csv", "")
  }
)
input.metadata <- data.frame(
  path = all.input.files,
  SampleID = to_vec(for (item in all.input.files){
    str_replace(basename(item), ".simplified.csv", "") 
  }),
  PROJECT = to_vec(for (item in all.input.files){
    str_split(item, "/")[[1]][[8]]
  })
) 

meta.data <- readxl::read_excel(file.path(path.to.05.output, "all_data_metadata.xlsx"))

mouse.id <- "m1"

tmp.metadata <- subset(meta.data, meta.data$mouse  == mouse.id) %>%
  subset(grepl("_HT", SampleID) == FALSE)

all.samples <- tmp.metadata$SampleID
bulk.samples <- to_vec(
  for (item in all.samples){
    if (grepl("MID", item) == TRUE){
      item
    }
  }
)

single.cell.samples <- setdiff(all.samples, bulk.samples)

if (file.exists(file.path(path.to.14.output, "bulkdf.csv")) == FALSE){
  bulkdf <- read_tsv(all.input.files[bulk.samples])[, c("id", "cloneCount")]
  bulkdf <- bulkdf %>%
    group_by(id) %>%
    summarise(sum.cloneCount = sum(cloneCount))
  colnames(bulkdf) <- c("id", "cloneCount")
  bulkdf <- bulkdf %>% rowwise() %>%
    mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
    mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
    mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]])
  
  write.table(bulkdf, 
              file.path(path.to.14.output, "bulkdf.csv"), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE) 
} else {
  bulkdf <- read.csv(file.path(path.to.14.output, "bulkdf.csv"), sep = "\t")
}

##### prepare data for input to circos plot
input.data <- list()
for (sample.id in single.cell.samples){
  input.data[[sample.id]] <- read.csv(all.input.files[[sample.id]], sep = "\t")
}
# input.data[["bulk"]] <- bulkdf

#####----------------------------------------------------------------------#####
##### CUSTOM: GENERATE CIRCOS PLOT INPUT FOR REAL CLONES AND PUBLIC CLONES
#####----------------------------------------------------------------------#####

##### input params
filter.clone = FALSE
filter.clone.cutoff = NULL
linkColor1 = "#FF000080"
linkColor2 = "lightgray"
group.to.highlight1 = NULL
group.to.highlight2 = NULL
ordered.samples = NULL

##### main function
cloneCountdf <- data.frame()
for (input.mid in names(input.data)){
  tmp.clonedf <- input.data[[input.mid]]
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
  tmp.clonedf <- subset(tmp.clonedf, select = c(id, cloneCount, Freq, accum.Freq))
  tmp.clonedf$SampleID <- input.mid
  tmp.clonedf <- rbind(data.frame(cloneCount = 0, id = "START", SampleID = input.mid, Freq = 0, accum.Freq = 0), tmp.clonedf)
  cloneCountdf <- rbind(cloneCountdf, tmp.clonedf)
}

count.clone.in.samples <- table(cloneCountdf$SampleID)
exclude.samples <- count.clone.in.samples[count.clone.in.samples <= 1] %>% names()
cloneCountdf <- subset(cloneCountdf, cloneCountdf$SampleID %in% exclude.samples == FALSE)
if (is.null(ordered.samples) == TRUE){
  cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = names(input.data))
} else {
  print(sprintf("Samples are ordered based on input order %s", paste(ordered.samples, collapse = ", ")))
  cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = ordered.samples)
  cloneCountdf <- cloneCountdf[order(cloneCountdf$SampleID), ]
}

all.plotdf <- list()
all.unique.clones <- c()
for (sample.id in levels(cloneCountdf$SampleID)){
  tmp.plotdf <- subset(cloneCountdf, cloneCountdf$SampleID == sample.id) %>%
    subset(select = c(id, cloneCount, Freq, accum.Freq))
  colnames(tmp.plotdf) <- c("id", 
                            sprintf("%s_Count", sample.id),
                            sprintf("%s_Freq", sample.id),
                            sprintf("%s_accumFreq", sample.id))
  all.plotdf[[sample.id]] <- tmp.plotdf
  all.unique.clones <- c(all.unique.clones, unique(tmp.plotdf$id)) %>% unique()
}

plotdf <- data.frame(id = all.unique.clones)
for (sample.id in names(all.plotdf)){
  plotdf <- merge(plotdf, all.plotdf[[sample.id]], by.x = "id", by.y = "id", all.x = TRUE)
}
plotdf[is.na(plotdf)] <- 0

##### at this point we have plotdf

subset(plotdf, plotdf$id %in% unique(match.scdf$VJcombi_CDR3_0.85))
