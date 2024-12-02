gc()
rm(list = ls())

install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.tar.gz", type = "source", repos = NULL)
new.pkgs <- c("APackOfTheClones", "svglite", "car", "ggpubr", "ggthemes")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
library(ggpubr)
library(ggthemes)
library(APackOfTheClones)

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

dataset.name <- "1st_dataset"
save.dev <- "png"

print(sprintf("Working on dataset %s", dataset.name))
s.obj <- readRDS(path.to.all.s.obj[[dataset.name]])
DefaultAssay(s.obj) <- "RNA"

path.to.02.output <- file.path(outdir, "GEX_output", "02_output", dataset.name)
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)
reduction.name <- "INTE_UMAP"

s.obj <- RunAPOTC(seurat_obj = s.obj, reduction_base = reduction.name, clonecall = "CTstrict")

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("cell.barcode")
clonedf <- data.frame(meta.data$CTstrict %>% table())
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count))
for (sampleid in unique(s.obj$name)){
  clonedf[[sampleid]] <- unlist(lapply(clonedf$clone, function(x){
    nrow(subset(meta.data, meta.data$name == sampleid & meta.data$CTstrict == x))
  }))
}

plot.top <- 20
for (plot.sampleid in unique(s.obj$name)){
  dir.create(file.path(path.to.02.output, 
                       dataset.name, 
                       "topClones_in_each_sample", 
                       plot.sampleid), showWarnings = FALSE, recursive = TRUE)
  
  tmp.clonedf <- clonedf[, c("clone", plot.sampleid)]
  colnames(tmp.clonedf) <- c("clone", "count")
  tmp.clonedf <- subset(tmp.clonedf, tmp.clonedf$count != 0)
  
  top20.clones <- head(tmp.clonedf, plot.top) %>% pull(clone)
  split.top5.clones <- split(top20.clones, seq(1,4))
  
  colors <- tableau_color_pal(palette = "Tableau 20")(20)
  split.colors <- split(colors, seq(1,4))
  
  apotc.clone.plot <- list()
  for (i in seq(1,4)){
    tmp.plot <- vizAPOTC(s.obj, clonecall = "CTaa", 
                         verbose = FALSE, 
                         reduction_base = reduction.name, 
                         show_labels = TRUE, 
                         legend_position = "top_right", 
                         legend_sizes = 2) %>% 
      showCloneHighlight(clonotype =  as.character(split.top5.clones[[i]]), 
                         fill_legend = TRUE,
                         color_each = split.colors[[i]], 
                         default_color = "#808080") 
    apotc.clone.plot[[i]] <- tmp.plot
  }
  library("gridExtra")
  merge.plot <- arrangeGrob(ggarrange(apotc.clone.plot[[1]], apotc.clone.plot[[2]],
                                      apotc.clone.plot[[3]], apotc.clone.plot[[4]], 
                                      ncol = 2, nrow = 2 ))
  
  ggsave(plot = merge.plot, filename = sprintf("%s_top%sClones.%s", plot.sampleid, plot.top, save.dev),
         path = file.path(path.to.02.output, dataset.name, "topClones_in_each_sample", plot.sampleid),
         device = save.dev, width = 20, height = 14
  )
}
