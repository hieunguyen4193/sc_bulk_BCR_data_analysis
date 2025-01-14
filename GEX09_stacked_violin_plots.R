gc()
rm(list = ls())

library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)

path.to.input.s.obj <- "/media/hieunguyen/HD0/outdir/CRC1382/FHager_datasets/20230513/220907_FH_v0.1/data_analysis/01_output/220907_FH_v0.1.CCR7_cluster.rds"
s.obj <- readRDS(path.to.input.s.obj)

DimPlot(object = s.obj, reduction = "RNA_UMAP", label = TRUE, label.box = TRUE, repel = TRUE)

##### table color - cluster
all.colors <- data.frame(color = hue_pal()(length(unique(s.obj@meta.data$seurat_clusters))))
all.colors$cluster <- seq(0, length(unique(s.obj@meta.data$seurat_clusters))-1)

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
features <- c("Xcr1", "Rab7b", "Clec9a", "Irf8", "Tlr12", "Ppt1", "Cadm1", "Agpat3", "Tlr3", "Plpp1")
path.to.save.output <- "/home/hieunguyen/CRC1382/outdir/tmp" 
cluster.order <- c(3,1,4,5)
p <- generate_vlnplot(obj = s.obj, features = features, pt.size = 0, cluster.order = cluster.order, 
                      colors = unlist(lapply(cluster.order, function(x){
                        return(subset(all.colors, all.colors$cluster == x)$color)
                      })))
# ggsave(plot = p, filename = "test.svg", path = path.to.save.output, device = "svg", width = 4, height = 10, dpi = 300)



