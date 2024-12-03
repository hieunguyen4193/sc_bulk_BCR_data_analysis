gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
input.dataset <- "241002_BSimons"
vdj.output <- readRDS(path.to.all.VDJ.output[[input.dataset]])

for (input.dataset in names(path.to.all.s.obj)){
  print(sprintf("Working on dataset %s", input.dataset))
  s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
  DefaultAssay(s.obj) <- "RNA"
  vdj.output <- vdj.output[unique(s.obj$name)]
  vdjdf <- data.frame()
  for (sample.id in names(vdj.output)){
    vdjdf <- rbind(vdjdf, vdj.output[[sample.id]])
  }
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
  if (length(intersect(meta.data$barcode, vdjdf$barcode)) == 0){
    vdjdf$barcode <- unlist(lapply(
      vdjdf$barcode, function(x){
        sampleid <- subset(vdjdf, vdjdf$barcode == x)$sample
        str_replace(x, sprintf("%s_%s", sampleid, sampleid), sampleid)
      }
    ))
  }
  
  vdjdf <- vdjdf %>% rowwise() %>%
    mutate(C.gene = ifelse(is.na(IGH) == FALSE, str_split(IGH, "[.]")[[1]][[4]], NA) )
  vdjdf <- subset(vdjdf, select = c(barcode, C.gene))
  meta.data <- merge(meta.data, vdjdf, by.x = "barcode", by.y = "barcode", all.x = TRUE)
  meta.data <- meta.data %>% column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, col.name = "C.gene", metadata = meta.data$C.gene)
  c.gene.barcodes <- list()
  for (g in setdiff(unique(s.obj$C.gene), c(NA, "NA"))){
    c.gene.barcodes[[g]] <- subset(meta.data, meta.data$C.gene == g) %>% row.names()
  }
  
  c.gene.umap <- DimPlot(object = s.obj, 
                         label = TRUE, label.box = TRUE, 
                         cells.highlight = c.gene.barcodes, 
                         cols.highlight = hue_pal()(length(names(c.gene.barcodes))))
  ggsave(plot = c.gene.umap, filename = sprintf("all_C_genes_marked_cells.svg"), 
         path = file.path(path.to.01.output, "svg", "module_scores"), device = "svg", width = 14, height = 10)
  
  for (g in names(c.gene.barcodes)){
    c.gene.umap <- DimPlot(object = s.obj, cells.highlight = c.gene.barcodes[[g]], 
                           cols.highlight = "red", label = TRUE, label.box = TRUE) + theme(legend.position = "none") + ggtitle(g)
    ggsave(plot = c.gene.umap, filename = sprintf("all_C_genes_marked_cells_%s.svg", g), 
           path = file.path(path.to.01.output, "svg", "module_scores"), device = "svg", width = 14, height = 10)
  }  
  
  count.C.gene.clusters <- table(s.obj$seurat_clusters, s.obj$C.gene) %>% as.data.frame() %>% pivot_wider(names_from = "Var1", values_from = "Freq")
  colnames(count.C.gene.clusters) <- c(c("Gene"), to_vec(for (i in seq(1, length(unique(s.obj$seurat_clusters)))) sprintf( 
    "cluster_%s", i)) )
  write.csv(count.C.gene.clusters, file.path(path.to.01.output, "svg", "module_scores", "count_C_gene_in_each_cluster.csv"))
  path.to.01.output <- file.path(outdir, "GEX_output", "01_output", input.dataset)
  dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(path.to.01.output, "svg", "module_scores"), showWarnings = FALSE, recursive = TRUE)
  
  tmpdf <- readxl::read_excel(file.path(path.to.project.src, "module_score_Bcells.xlsx"))
  module.gene.list <- list() 
  for (g in colnames(tmpdf)){
    tmp <- tmpdf[[g]] %>% unique()
    g <- str_replace_all(g, " ", "_")
    module.gene.list[[g]] <- tmp[is.na(tmp) == FALSE]
  }
  
  for (input.list in names(module.gene.list)){
    DefaultAssay(s.obj) <- "RNA"
    s.obj <- AddModuleScore(object = s.obj, features = list(module.gene.list[[input.list]]), name = sprintf("%s_", input.list), ctrl = 50)
  }
  
  fake.module.gene.list <- to_vec(
    for (item in names(module.gene.list)){
      sprintf("%s_1", item)
    }
  )
  library(viridis)
  if ("svglite" %in% installed.packages() == FALSE){
    install.packages("svglite")
  }
  if (input.dataset %in% c("241002_BSimons", "241104_BSimons")){
    reduction.name <- "RNA_UMAP"
  } else {
    reduction.name <- "INTE_UMAP"
  }

  DefaultAssay(s.obj) <- "RNA"
  Idents(s.obj) <- "seurat_clusters"
  
  feature.plot <- FeaturePlot(object = s.obj, reduction = reduction.name, ncol = 3, features = fake.module.gene.list, pt.size = 1) &
    scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
  
  violin.plot <- VlnPlot(object = s.obj, features = fake.module.gene.list, pt.size = 0)
  
  heatmapdf <- s.obj@meta.data[, c("seurat_clusters", fake.module.gene.list)] %>%
    group_by(seurat_clusters) %>%
    summarise(across(everything(), mean)) %>%
    column_to_rownames("seurat_clusters")
  
  colors <- rev(RColorBrewer::brewer.pal(n = length(fake.module.gene.list), name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)
  
  heatmapdf.scaled <- (heatmapdf - rowMeans(heatmapdf))/rowSds(as.matrix(heatmapdf))
  colnames(heatmapdf.scaled) <- to_vec(
    for (item in colnames(heatmapdf.scaled)){
      str_split(item, "_")[[1]][[1]]
    }
  )
  heatmap.plot <- heatmapdf.scaled %>% rownames_to_column("cluster") %>% 
    pivot_longer(!cluster, names_to = "signature", values_to = "z_score") %>%
    ggplot(aes(x = cluster, y = signature, fill = z_score)) + geom_tile() + 
    scale_fill_distiller(palette = "RdBu") + 
    theme(axis.text = element_text(size = 22))
  
  ggsave(plot = feature.plot, filename = sprintf("feature_plot_module_scores.svg"), 
         path = file.path(path.to.01.output, "svg", "module_scores"), device = "svg", width = 20, height = 14)
  ggsave(plot = violin.plot, filename = sprintf("violin_plot_module_scores.svg"), 
         path = file.path(path.to.01.output, "svg", "module_scores"), device = "svg", width = 20, height = 14)
  ggsave(plot = heatmap.plot, filename = sprintf("heatmap_module_scores.svg"), 
         path = file.path(path.to.01.output, "svg", "module_scores"), device = "svg", width = 10, height = 10)
}

