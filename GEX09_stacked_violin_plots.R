gc()
rm(list = ls())

new.pkgs <- c("svglite", "dplyr")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
library(patchwork)
library(scales)
library(ggthemes)
library(Seurat)
library("gridExtra")
library(ggplot2)
library(dplyr)
library(tidyverse)
path.to.s.obj <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1/240805_BSimons/data_analysis/07_output/quantile_0.85/all_samples/s8_output/240805_BSimons.renamedClusters.output.s8.rds"
dataset.name <- "240805_BSimons_filterHT_cluster_renamed"
save.dev <- "svg"

print(sprintf("Working on dataset %s", dataset.name))

s.obj <- readRDS(path.to.s.obj)
DefaultAssay(s.obj) <- "RNA"
Idents(s.obj) <- "seurat_clusters"

if (dataset.name %in% c("241104_BSimons", "241002_BSimons")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

path.to.save.output <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1/tmp"
path.to.save.output <- file.path(path.to.save.output, save.dev)

dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

##### table color - cluster
modify_vlnplot_from_seurat<- function(obj, 
                                      feature, 
                                      colors,
                                      pt.size = 0, 
                                      plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                                      cluster.order = NULL,
                                      plot.for.CCR7 = FALSE) {
  if (plot.for.CCR7 == TRUE){
    Idents(obj) <- "ccr7_cluster"
  } else {
    Idents(obj) <- "seurat_clusters"
  }
  if (is.null(cluster.order) == TRUE){
    p <- VlnPlot(obj, features = feature, pt.size = pt.size)
  } else {
    tmp.obj <- subset(obj, seurat_clusters %in% cluster.order)
    tmp.obj$new.ident <- factor(tmp.obj$seurat_clusters, levels = cluster.order )
    Idents(s.obj) <- "new.ident"
    p <- VlnPlot(tmp.obj, features = feature, pt.size = pt.size, group.by = "new.ident")
  }
  
  p <- p + scale_fill_manual(values = colors) +
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max <- function(p){
  ymax <- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
generate_vlnplot <- function(obj, 
                             features,
                             colors,
                             pt.size = 0,
                             cluster.order = NULL,
                             plot.for.CCR7 = FALSE) {
  
  plot_list <- purrr::map(features, function(x) modify_vlnplot_from_seurat(obj = obj,
                                                                           feature = x, 
                                                                           cluster.order = cluster.order, 
                                                                           colors = colors, 
                                                                           plot.for.CCR7 = plot.for.CCR7))
  
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(), axis.ticks.x = element_line())
  
  ymaxs <- purrr::map_dbl(plot_list, extract_max)
  plot_list <- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                             scale_y_continuous(breaks = c(y)) + 
                             expand_limits(y = y))
  
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

##### stacked violin plot

gene.list <- list(
  `240805_BSimons_filterHT_cluster_renamed` = list(  group1 = c("Nr4a1", "Klf4", "Cd86", "Cd83", "Ccr6"),
                                                     group2= c("Foxp1", "Cxcr4", "Cxcr5", "Nt5e"))
)


all.clusters <- sort(unique(s.obj$seurat_clusters))

cluster.order <- c(9,8,7,6,5,1,2,3,4)

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("barcode")
meta.data$seurat_clusters <- factor(meta.data$seurat_clusters, levels = cluster.order)
meta.data <- meta.data %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, col.name = "seurat_clusters", metadata = meta.data$seurat_clusters)

plot.clusters <- tableau_color_pal(palette = "Tableau 20")(length(all.clusters)) 
names(plot.clusters) <- cluster.order

p.umap <- DimPlot(object = s.obj, reduction = reduction.name, label = TRUE, label.box = TRUE, 
                  group.by = "seurat_clusters", cols = plot.clusters)

ggsave(plot = p.umap, filename = sprintf("UMAP.%s", save.dev), 
       path = path.to.save.output, device = save.dev, width = 14, height = 10, dpi = 300)  

for (g in names(gene.list[[dataset.name]])){
  p <- generate_vlnplot(obj = s.obj, features = gene.list[[dataset.name]][[g]], pt.size = 0, cluster.order = cluster.order, 
                        colors = plot.clusters)
  ggsave(plot = p, filename = sprintf("violin_plot_group_%s.svg", g), 
         path = path.to.save.output, device = save.dev, width = 4, height = 10, dpi = 300)  
}
