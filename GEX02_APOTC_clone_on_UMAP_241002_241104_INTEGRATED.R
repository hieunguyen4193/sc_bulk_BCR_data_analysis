gc()
rm(list = ls())
# install.packages("dplyr")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.tar.gz", type = "source", repos = NULL)
new.pkgs <- c("APackOfTheClones", "svglite", "car", "ggpubr", "ggthemes", "dplyr")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
library(ggpubr)
library(ggthemes)
library(APackOfTheClones)
library("gridExtra")

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
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.all.s.obj <- path.to.all.s.obj[setdiff(names(path.to.all.s.obj), c("BonnData"))]

sc.projects.with.ht <- c("240805_BSimons_filterHT_cluster_renamed", "241002_241104_BSimons")
name_or_sampleHT <- "sample_ht"

if (name_or_sampleHT == "sample_ht"){
  path.to.all.s.obj <- path.to.all.s.obj[sc.projects.with.ht]
}

clone.name <- "VJcombi_CDR3_0.85"

dataset.name <- "241002_241104_BSimons"

topN <- 5
save.dev <- "svg"
# save.dev <- "tiff"
get.YFP.clone.only <- TRUE
##### sample.list for each mouse:

print(sprintf("Working on dataset %s", dataset.name))

s.obj <- readRDS(path.to.all.s.obj[[dataset.name]])

if (dataset.name %in% sc.projects.with.ht){
  meta.data <- s.obj@meta.data %>% 
    rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(sample_ht = sprintf("%s_%s", name, HTO_classification)) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$sample_ht, col.name = "sample_ht")      
}

DefaultAssay(s.obj) <- "RNA"

path.to.02.output <- file.path(outdir, 
                               "GEX_output", 
                               sprintf("02_output_241002_241104_BSimons"), 
                               dataset.name, 
                               clone.name,
                               sprintf("top%s", topN),
                               save.dev)
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)
if (dataset.name %in% c("241104_BSimons", "241002_BSimons")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

Idents(s.obj) <- "seurat_clusters"
s.obj <- RunAPOTC(seurat_obj = s.obj, reduction_base = reduction.name, clonecall = clone.name)

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("cell.barcode")
clonedf <- data.frame(meta.data[[clone.name]] %>% table())
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count))

for (sampleid in unique(s.obj@meta.data[[name_or_sampleHT]])){
  print(sprintf("Generating clonedf, working on sample: %s", sampleid))
  clonedf[[sampleid]] <- unlist(lapply(clonedf$clone, function(x){
    nrow(subset(meta.data, meta.data[[name_or_sampleHT]] == sampleid & meta.data[[clone.name]] == x))
  }))
}

if (get.YFP.clone.only == TRUE){
  count.yfp.exprs <- GetAssayData(object = s.obj, slot = "counts", assay = "RNA")["YFP", ]
  yfp.cells <- count.yfp.exprs[count.yfp.exprs != 0] %>% names()
  non.yfp.cells <- count.yfp.exprs[count.yfp.exprs == 0] %>% names()
  
  print("Generating data for YFP clones only, clones that have at least 1 cell YFP+")
  path.to.02.output <- file.path(path.to.02.output, "YFP_clones_only")
  dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)
  
  clonedf <- clonedf %>% rowwise() %>%
    mutate(clone.size = length(row.names(subset(s.obj@meta.data, s.obj@meta.data$VJcombi_CDR3_0.85 == clone)))) %>%
    mutate(YFP.clone = 
             length(intersect(
               row.names(subset(s.obj@meta.data, s.obj@meta.data$VJcombi_CDR3_0.85 == clone)), yfp.cells
             ))
    ) %>%
    mutate(non.YFP.clone = 
             length(intersect(
               row.names(subset(s.obj@meta.data, s.obj@meta.data$VJcombi_CDR3_0.85 == clone)), non.yfp.cells
             ))
    )
  clonedf <- subset(clonedf, clonedf$YFP.clone != 0)
}

sample.list <- list(
  PP3 = c("PP3_HT1", "PP3_HT2", "PP3_HT3"),
  PP7 = c("PP7_HT1", "PP7_HT2", "PP7_HT3")
)

for (sample.list.name in names(sample.list)){
  if (file.exists(file.path(path.to.02.output, sprintf("APOTC_%s.%s", sample.list.name, save.dev))) == FALSE){
    plot.colors <- tableau_color_pal(palette = "Tableau 20")(20)
    plot.clonedf <- data.frame()
    for (sample.id in sample.list[[sample.list.name]]){
      tmpdf <- clonedf[, c("clone", sample.id)]
      colnames(tmpdf) <- c("clone", "SampleID")
      tmpdf <- tmpdf %>% arrange(desc(SampleID)) %>% head(topN)
      tmp.plot.clonedf <- data.frame(clone = tmpdf$clone)
      tmp.plot.clonedf$SampleID <- sample.id
      plot.clonedf <- rbind(plot.clonedf, tmp.plot.clonedf)
    }
    plot.clonedf$color <- head(plot.colors, nrow(plot.clonedf))
    
    tmp.plot <- vizAPOTC(s.obj, clonecall = clone.name, 
                         verbose = FALSE, 
                         reduction_base = reduction.name, 
                         show_labels = TRUE, 
                         repulsion_strength = 5,
                         legend_position = "top_right", 
                         legend_sizes = 2) %>% 
      showCloneHighlight(clonotype =  as.character(plot.clonedf$clone), 
                         fill_legend = TRUE,
                         color_each = plot.clonedf$color, 
                         default_color = "lightgray") 
    ggsave(plot = tmp.plot, filename = sprintf("APOTC_%s.%s", sample.list.name, save.dev), path = path.to.02.output, dpi = 300, width = 14, height = 10)        
  } else {
    print(sprintf("File %s exists, not generating a new plot ...", 
                  file.path(path.to.02.output, sprintf("APOTC_%s.%s", sample.list.name, save.dev))))
  }
}

if (get.YFP.clone.only == TRUE){
  print("Generate APOTC for YFP clones only for 2 mice, color by mice basis ...")
  plot.clonedf <- data.frame()
  for (sample.list.name in names(sample.list)){
    for (sample.id in sample.list[[sample.list.name]]){
      tmpdf <- clonedf[, c("clone", sample.id)]
      colnames(tmpdf) <- c("clone", "SampleID")
      tmpdf <- tmpdf %>% arrange(desc(SampleID)) %>% subset(SampleID != 0)
      tmp.plot.clonedf <- data.frame(clone = tmpdf$clone)
      tmp.plot.clonedf$SampleID <- sample.id
      plot.clonedf <- rbind(plot.clonedf, tmp.plot.clonedf)
    }    
  }
  plot.colors <- tableau_color_pal(palette = "Tableau 10")(2)
  plot.clonedf$color <- to_vec(
    for (item in plot.clonedf$SampleID){
      if (grepl("PP3", item) == TRUE){
        plot.colors[[1]]
      } else {
        plot.colors[[2]]
      }
    }
  )
  tmp.plot <- vizAPOTC(s.obj, clonecall = clone.name, 
                       verbose = FALSE, 
                       reduction_base = reduction.name, 
                       show_labels = TRUE, 
                       repulsion_strength = 5,
                       legend_position = "top_right", 
                       legend_sizes = 2) %>% 
    showCloneHighlight(clonotype =  as.character(plot.clonedf$clone), 
                       fill_legend = TRUE,
                       color_each = plot.clonedf$color, 
                       default_color = "lightgray") +
    ggtitle(sprintf("PP3: %s, PP7: %s", plot.colors[[1]], plot.colors[[2]]))
  ggsave(plot = tmp.plot, filename = sprintf("APOTC_2_mice_YFP_clones_only.%s", save.dev), path = path.to.02.output, dpi = 300, width = 14, height = 10)        
}
