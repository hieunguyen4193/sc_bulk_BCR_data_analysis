gc()
rm(list = ls())

outdir <- "/home/hieu/outdir"
PROJECT <- "mixcr_pipeline_output"

path.to.mid.metadata <- file.path(outdir, "FT_output", "mid_labels.csv") # <<<<< change here, input mid_labels.csv file. 
mid.metadata <- read.csv(path.to.mid.metadata, sep = ";")

path.to.project.src <- "/home/hieu/src/BCRTree_release/gctree" # <<<<< change here to current working source dir
path.to.save.output <- file.path(path.to.project.src, "sample_groups") 
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

##### mouse based output 
dir.create(file.path(path.to.save.output, "mouse_based_output"), showWarnings = FALSE, recursive = TRUE)
for (mouse.id in unique(mid.metadata$mouse)){
  tmpdf <- subset(mid.metadata, mid.metadata$mouse == mouse.id)
  write.table(tmpdf, sep = "\t", col.names = FALSE, row.names = FALSE,
              file.path(path.to.save.output, "mouse_based_output", sprintf("mouse_%s.tsv", mouse.id)))
}

dir.create(file.path(path.to.save.output, "mouse_YFP_based_output"), showWarnings = FALSE, recursive = TRUE)
for (mouse.id in unique(mid.metadata$mouse)){
  tmpdf <- subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)
  write.table(tmpdf, sep = "\t", col.names = FALSE, row.names = FALSE,
              file.path(path.to.save.output, "mouse_YFP_based_output", sprintf("mouse_%s.tsv", mouse.id)))
}

dir.create(file.path(path.to.save.output, "mouse_YFPpos_based_output"), showWarnings = FALSE, recursive = TRUE)
for (mouse.id in unique(mid.metadata$mouse)){
  tmpdf <- subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)
  write.table(tmpdf, sep = "\t", col.names = FALSE, row.names = FALSE,
              file.path(path.to.save.output, "mouse_YFPpos_based_output", sprintf("mouse_%s.tsv", mouse.id)))
}

dir.create(file.path(path.to.save.output, "mouse_YFPneg_based_output"), showWarnings = FALSE, recursive = TRUE)
for (mouse.id in unique(mid.metadata$mouse)){
  tmpdf <- subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)
  write.table(tmpdf, sep = "\t", col.names = FALSE, row.names = FALSE,
              file.path(path.to.save.output, "mouse_YFPneg_based_output", sprintf("mouse_%s.tsv", mouse.id)))
}