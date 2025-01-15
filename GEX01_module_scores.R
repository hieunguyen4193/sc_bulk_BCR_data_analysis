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

library(viridis)
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
all.module.genedf <- list(
  old.version = file.path(path.to.project.src, "module_score_Bcells.xlsx"),
  ver05122024 = file.path(path.to.project.src, "module_score_Bcells.05122024.xlsx"),
  ver07012025 = file.path(path.to.project.src, "module_score_07012025.xlsx"),
  ver09012025 = file.path(path.to.project.src, "module_score_09012025.xlsx"),
  ver09012025_noBcells = file.path(path.to.project.src, "module_score_09012025_noBcell.xlsx"),
  ver09012025_2 = file.path(path.to.project.src, "module_score_09012025_Bcells_to_others.xlsx")
)

module.score.version <- "ver09012025_noBcells"
module.genedf <- readxl::read_excel(all.module.genedf[[module.score.version]])

# save.figure.type <- "tiff"
save.figure.type <- "svg"

for (input.dataset in names(path.to.all.s.obj)){
  path.to.01.output <- file.path(outdir, "GEX_output", sprintf("01_output_%s", module.score.version), input.dataset)
  if (file.exists(file.path(path.to.01.output, sprintf("finished_%s_%s", input.dataset, save.figure.type))) == FALSE){
    if (input.dataset %in% c("241002_BSimons", "241104_BSimons", "BonnData")){
      reduction.name <- "RNA_UMAP"
    } else {
      reduction.name <- "INTE_UMAP"
    }

    dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(path.to.01.output, save.figure.type, "module_scores"), showWarnings = FALSE, recursive = TRUE)
    
    print(sprintf("Working on dataset %s", input.dataset))
    s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
    
    input.colors <- tableau_color_pal(palette = "Tableau 20")(length(unique(s.obj$seurat_clusters)))
    umap.p <- DimPlot(object = s.obj, 
                      reduction = reduction.name,
                      label = TRUE, 
                      label.box = TRUE, 
                      repel = TRUE, 
                      group.by = "seurat_clusters", 
                      cols = input.colors)
    ggsave(plot = umap.p, filename = sprintf("UMAP.%s", save.figure.type), 
           path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 14, height = 10)
    
    DefaultAssay(s.obj) <- "RNA"
    if (input.dataset == "BonnData"){
      vdj.output <- readRDS(path.to.all.VDJ.output[[input.dataset]])
    } else if (input.dataset %in% names(all.integration.cases)){
      vdj.output <- c()
      for (tmp.dataset in all.integration.cases[[input.dataset]]){
        tmp.vdj.output <- readRDS(path.to.all.VDJ.output[[tmp.dataset]])
        if (tmp.dataset == "1st_dataset"){
          tmp.vdj.output <- tmp.vdj.output[dataset1_sample_list]
        } else if (tmp.dataset == "2nd_dataset"){
          tmp.vdj.output <- tmp.vdj.output[dataset2_sample_list]
        }
        vdj.output <- c(vdj.output, tmp.vdj.output)
      }
    } else {
      vdj.output <- readRDS(path.to.all.VDJ.output[[input.dataset]])
      vdj.output <- vdj.output[unique(s.obj$name)]
    }
    
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
    
    c.gene.umap <- DimPlot(object = s.obj, reduction = reduction.name, 
                           label = TRUE, label.box = TRUE, 
                           cells.highlight = c.gene.barcodes, 
                           cols.highlight = hue_pal()(length(names(c.gene.barcodes))))
    ggsave(plot = c.gene.umap, filename = sprintf("all_C_genes_marked_cells.%s", save.figure.type), 
           path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 14, height = 10)
    
    for (g in names(c.gene.barcodes)){
      c.gene.umap <- DimPlot(object = s.obj, cells.highlight = c.gene.barcodes[[g]], reduction = reduction.name,
                             cols.highlight = "red", label = TRUE, label.box = TRUE) + theme(legend.position = "none") + ggtitle(g)
      ggsave(plot = c.gene.umap, filename = sprintf("all_C_genes_marked_cells_%s.%s", g, save.figure.type), 
             path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 14, height = 10)
    }  
    
    count.C.gene.clusters <- table(s.obj$seurat_clusters, s.obj$C.gene) %>% as.data.frame() %>% pivot_wider(names_from = "Var1", values_from = "Freq")
    colnames(count.C.gene.clusters) <- c(c("Gene"), to_vec(for (i in seq(1, length(unique(s.obj$seurat_clusters)))) sprintf( 
      "cluster_%s", i)) )
    write.csv(count.C.gene.clusters, file.path(path.to.01.output, save.figure.type, "module_scores", "count_C_gene_in_each_cluster.csv"))
    
    module.gene.list <- list() 
    single.module.genes <- c()
    for (g in colnames(module.genedf)){
      tmp <- module.genedf[[g]] %>% unique()
      g <- str_replace_all(g, " ", "_")
      module.gene.list[[g]] <- tmp[is.na(tmp) == FALSE]
      single.module.genes <- c(single.module.genes, module.gene.list[[g]])
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
    DefaultAssay(s.obj) <- "RNA"
    Idents(s.obj) <- "seurat_clusters"
    
    feature.plot <- FeaturePlot(object = s.obj, reduction = reduction.name, label = TRUE, ncol = 3, features = fake.module.gene.list, pt.size = 1) &
      scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")
    
    violin.plot <- VlnPlot(object = s.obj, features = fake.module.gene.list, pt.size = 0)
    
    heatmapdf <- s.obj@meta.data[, c("seurat_clusters", fake.module.gene.list)] %>%
      group_by(seurat_clusters) %>%
      summarise(across(everything(), mean)) %>%
      column_to_rownames("seurat_clusters")
    
    colors <- rev(RColorBrewer::brewer.pal(n = length(fake.module.gene.list), name = "RdBu"))
    colors.use <- grDevices::colorRampPalette(colors = colors)(100)
    
    # heatmapdf.scaled <- (heatmapdf - colMins(as.matrix(heatmapdf)))/(colMaxs(as.matrix(heatmapdf)) - colMins(as.matrix(heatmapdf))) 
    
    heatmapdf <- heatmapdf %>% t() %>% as.data.frame()
    heatmapdf.scaled <- (heatmapdf - rowMeans(heatmapdf))/rowSds(as.matrix(heatmapdf))
    
    colnames(heatmapdf.scaled) <- to_vec(
      for (item in colnames(heatmapdf.scaled)){
        str_replace(item, "_1", "")
      }
    )
    for (c in row.names(heatmapdf.scaled)){
      heatmap.plot <- heatmapdf.scaled[c, ] %>% rownames_to_column("cluster") %>% 
        pivot_longer(!cluster, names_to = "signature", values_to = "z_score") %>%
        ggplot(aes(x = cluster, y = signature, fill = z_score)) + geom_tile() + 
        scale_fill_distiller(palette = "RdBu") + 
        theme(axis.text.x = element_text(size = 22, angle = 90),
              axis.text.y = element_text(size = 22))
      ggsave(plot = heatmap.plot, filename = sprintf("heatmap_module_scores.%s.%s", c, save.figure.type), 
             path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 10, height = 10)
    }

    ggsave(plot = feature.plot, filename = sprintf("feature_plot_module_scores.%s", save.figure.type), 
           path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 20, height = 14)
    ggsave(plot = violin.plot, filename = sprintf("violin_plot_module_scores.%s", save.figure.type), 
           path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 20, height = 14)
    
    
    ##### single gene plots
    single.module.genes <- intersect(single.module.genes, row.names(s.obj))
    # gex.mat <- GetAssayData(object = s.obj, slot = "data", assay = "RNA")[single.module.genes, ] %>%
    #   t() %>% data.frame() %>%
    #   rownames_to_column("barcode")
    # meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
    # gex.mat <- merge(gex.mat, meta.data, by.x = "barcode", by.y = "barcode")
    # input.heatmap <- gex.mat[, c(single.module.genes, "seurat_clusters")] %>% 
    #   group_by(seurat_clusters) %>% summarise_all("mean") %>%
    #   column_to_rownames("seurat_clusters") %>% 
    #   t() %>% 
    #   data.frame()
    # colnames(input.heatmap) <- c(seq(1, length(unique(s.obj$seurat_clusters))) )
    # 
    # heatmapdf.scaled <- (input.heatmap - rowMeans(input.heatmap))/rowSds(as.matrix(input.heatmap))
    # heatmap.plot.single.gene <- heatmapdf.scaled %>% rownames_to_column("cluster") %>% 
    #   pivot_longer(!cluster, names_to = "signature", values_to = "z_score") %>%
    #   ggplot(aes(x = cluster, y = signature, fill = z_score)) + geom_tile() + 
    #   scale_fill_distiller(palette = "RdBu") + 
    #   theme(axis.text = element_text(size = 22, angle = 90))
    
    heatmap.plot.single.gene <- DoHeatmap(object = s.obj, 
                                          assay = "RNA", 
                                          features = single.module.genes, 
                                          group.by = "seurat_clusters", draw.lines = TRUE
    ) +
      scale_fill_gradient(low = "white", high = "red")
    ggsave(plot = heatmap.plot.single.gene, filename = sprintf("heatmap_all_module_genes.%s", save.figure.type), 
           path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 20, height = 15)
    split.single.module.genes <- split(single.module.genes, ceiling(seq_along(single.module.genes) / 9))
    for (g in names(split.single.module.genes)){
      input.genes <- split.single.module.genes[[g]]
      feature.plot <- FeaturePlot(object = s.obj, 
                                  reduction = reduction.name, label = TRUE, 
                                  ncol = 3,
                                  features = input.genes, 
                                  pt.size = 1) &
        scale_color_gradient(low = "lightgray", 
                             high = "#FF0000", 
                             na.value = "lightgray")
      violin.plot <- VlnPlot(object = s.obj, features = input.genes, pt.size = 0)
      ggsave(plot = feature.plot, filename = sprintf("feature_plot_module_single_gene_%s.%s", g, save.figure.type), 
             path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 20, height = 14)
      ggsave(plot = violin.plot, filename = sprintf("violin_plot_module_single_gene_%s.%s", g, save.figure.type), 
             path = file.path(path.to.01.output, save.figure.type, "module_scores"), device = save.figure.type, width = 20, height = 14)
    }
    write.csv(data.frame(status = c("finished")), 
              file.path(path.to.01.output, sprintf("finished_%s_%s", input.dataset, save.figure.type)))    
  } else {
    print(sprintf("%s module score exists", input.dataset))
  }
}


