gc()
rm(list = ls())
library(tidyr)
library(dplyr)
library(comprehenr)
library(stringr)
source("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/circos.R")

# path.to.input <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets/test_circos_inputs"
# path.to.save.output <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets/test_circos_outputs"

path.to.input <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets/hieudata_"

all.files <- Sys.glob(file.path(path.to.input, "*"))
names(all.files) <- to_vec(
  for (item in basename(all.files)){
    str_replace(item, "_clones.csv", "")
  }
)

mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")

for (mouse.id in unique(mid.metadata$mouse)){
  path.to.save.output <- sprintf("/media/hieunguyen/HNSD01/storage/all_BSimons_datasets/hieudata_%s_", mouse.id)
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  list.of.samples <- tolower(subset(mid.metadata, mid.metadata$mouse == mouse.id)$X)
  names(list.of.samples) <- to_vec(
    for (item in list.of.samples){
      subset(mid.metadata, mid.metadata$X == toupper(item))$population
    }
  )
  files <- all.files[list.of.samples]
  circos(files = files,
         fileAliases = names(list.of.samples),
         saveFolder = path.to.save.output,
         cutoff = 1.0,
         sort = FALSE,
         countColors = c("#FFFFFFFF", "#0000FFFF"),
         linkColors = rep("#FF000080", ifelse(length(files) == 1, 1, length(combn(length(files), 2)))),
         showLinks = rep(TRUE, length(combn(length(files), 2))))
}

