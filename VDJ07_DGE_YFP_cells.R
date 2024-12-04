gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
# source(file.path(path.to.main.src, "VDJ_helper_functions.R"))
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
library(circlize)
#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/GSHD_HN01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
sc.projects <-  c("241002_BSimons", "241104_BSimons", "2nd_dataset")

for (PROJECT in sc.projects){
  print(sprintf("Working on project %s", PROJECT))
  path.to.07.output <- file.path(outdir, "VDJ_output", "07_output", PROJECT)
  dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- readRDS(path.to.all.s.obj[[PROJECT]])
  DefaultAssay(s.obj) <- "RNA"
  all.genes <- row.names(s.obj)
  
  all.cells <- GetAssayData(object = s.obj, assay = "RNA", slot = "data")["YFP", ]
  
  yfp.cells <- all.cells[all.cells != 0] %>% names()
  non.yfp.cells <- all.cells[all.cells == 0] %>% names()
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(YFP = ifelse(barcode %in% yfp.cells, "yes", "no")) %>%
    column_to_rownames("barcode")
  
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$YFP, col.name = "YFP")
  
  yfp.de.markers.raw <- FindMarkers(object = s.obj, group.by = "YFP", ident.1 = "yes", ident.2 = "no", assay = "RNA", test.use = "wilcox",
                                    features = setdiff(all.genes, "YFP"))
  saveRDS(yfp.de.markers.raw, file.path(path.to.07.output, "yfp.de.markers.raw.rds"))
  yfp.de.markers <- yfp.de.markers.raw %>% 
    rownames_to_column("Gene") %>%
    subset(p_val_adj <= 0.05) %>% 
    rowwise() %>%
    mutate(abs.avg_log2FC = abs(avg_log2FC)) %>%
    arrange(desc(avg_log2FC))
  
  input.df <- yfp.de.markers.raw %>%
    rownames_to_column("Gene") %>%
    rowwise() %>%
    mutate(sig = ifelse(p_val_adj <= 0.05, "sig", "non-sig")) %>%
    mutate(show.gene.name = ifelse(sig == "sig", Gene, NA)) %>%
    mutate(abs.avg_log2FC = abs(avg_log2FC))
  
  cutoff.adjp <- 0.05
  
  write.csv(input.df, file.path(path.to.07.output, sprintf("input_VolcanoPlot_YFP_%s_cells_vs_nonYFP_%s_cells.svg", length(yfp.cells), length(non.yfp.cells))))
  volcano.plot <- ggplot(data=input.df, 
                         aes(x=avg_log2FC, y=-log10(p_val_adj), col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
    geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
    geom_text_repel() +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    xlim(c(-max(input.df$abs.avg_log2FC), max(input.df$abs.avg_log2FC)))
  
  ggsave(plot = volcano.plot, 
         filename = sprintf("volcano_plot_YFP_%s_cells_vs_nonYFP_%s_cells.svg", length(yfp.cells), length(non.yfp.cells)), 
         path = path.to.07.output, 
         dpi = 300, 
         width = 14, 
         height = 10)
}
