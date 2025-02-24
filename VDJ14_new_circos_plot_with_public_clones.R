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

path.to.tree.05.output <- file.path(outdir, 
                                    "GEX_output", 
                                    "05_output", 
                                    sc.projects[[1]], 
                                    "input_circos_public_clones")
all.sc.files <- Sys.glob(file.path(path.to.tree.05.output, "*.csv"))
names(all.sc.files) <- to_vec(for (item in all.sc.files){
  str_replace(basename(item), "_publicClone.csv", "")
})

path.to.tree.01.output <- file.path(outdir, 
                                    "tree_analysis",
                                    bulk.projects[[1]], 
                                    "01_output", 
                                    "input_circos_public_clones")
all.bulk.files <- Sys.glob(file.path(path.to.tree.01.output, "*.csv"))
names(all.bulk.files) <- to_vec(for (item in all.bulk.files){
  str_replace(basename(item), "_publicClone.csv", "")
})

filter.clone = FALSE
filter.clone.cutoff = NULL
linkColor1 = "#FF000080"
linkColor2 = "lightgray"
group.to.highlight1 = NULL
group.to.highlight2 = NULL
ordered.samples = NULL

fileAliases <- c("M3", "P3", "m3")
names(fileAliases) <- c("M3", "P3", "m3")

input.files <- c(all.sc.files[c("M3", "P3")], all.bulk.files[c("m3")])

cloneCountdf <- data.frame()
for (input.mid in names(input.files)){
  tmpdf <- read.csv(input.files[[input.mid]]) 
  if (input.mid %in% all.bulk.files){
    tmpdf <- tmpdf[, c("cloneid", "clone_sum_abundance", "match_public_clone")]
  }
  colnames(tmpdf) <- c("id", "cloneCount", "match_public_clone")
  tmp.clonedf <- tmpdf
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
  tmp.clonedf <- subset(tmp.clonedf, select = c(id, cloneCount, Freq, accum.Freq, match_public_clone))
  tmp.clonedf$SampleID <- input.mid
  tmp.clonedf <- rbind(data.frame(cloneCount = 0, id = "START", SampleID = input.mid, Freq = 0, accum.Freq = 0, match_public_clone = "no"), tmp.clonedf)
  cloneCountdf <- rbind(cloneCountdf, tmp.clonedf)
}

count.clone.in.samples <- table(cloneCountdf$SampleID)
exclude.samples <- count.clone.in.samples[count.clone.in.samples <= 1] %>% names()
cloneCountdf <- subset(cloneCountdf, cloneCountdf$SampleID %in% exclude.samples == FALSE)
keep.samples <- setdiff(names(input.files), exclude.samples)
if (is.null(ordered.samples) == TRUE){
  cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = keep.samples)
} else {
  print(sprintf("Samples are ordered based on input order %s", paste(ordered.samples, collapse = ", ")))
  cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = ordered.samples)
  cloneCountdf <- cloneCountdf[order(cloneCountdf$SampleID), ]
}

new.fileAliases <- to_vec(
  for (item in keep.samples){
    fileAliases[[item]]
  }
)
fileAliases <- new.fileAliases
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

# tmp.clonedf.plot <- tmp.clonedf[, c("id", "Freq", "accum.Freq")]
# colnames(tmp.clonedf.plot) <- c("id", sprintf("%s_Freq", input.mid), sprintf("%s_accum.Freq", input.mid))
# matched.public.clones <- subset(tmp.clonedf, tmp.clonedf$match_public_clone == "yes")$id
# for (p.clone in matched.public.clones){
#   tmp.clonedf.plot[[sprintf("%s_Freq", p.clone)]] <- to_vec(
#     for (item in tmp.clonedf.plot$id){
#       if (item == p.clone){
#         return(1)
#       } else {
#         return(0)
#       }
#     }
#   )
#   tmp.clonedf.plot[[sprintf("%s_accum.Freq", p.clone)]] <- to_vec(
#     for (item in tmp.clonedf.plot$id){
#       if (item == p.clone){
#         return(1)
#       } else {
#         return(0)
#       }
#     }
#   )
# }
# 
# plotdf <- tmp.clonedf.plot
