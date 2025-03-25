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
if ("ggthemes" %in% installed.packages() == FALSE){
  install.packages("ggthemes")
}
library(ggthemes)
#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

path.to.12.output <- file.path(outdir, "GEX_output", "12_output")
dir.create(path.to.12.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.12.output, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.12.output, "volcano_plots"), showWarnings = FALSE, recursive = TRUE)

sc.dataset.with.hto <- c("240805_BSimons",
                         "241104_BSimons",
                         "240805_BSimons_filterHT",
                         "240805_BSimons_filterHT_cluster_renamed",
                         "241002_BSimons",
                         "240805_BSimons_filterHT_cluster",
                         "241002_241104_BSimons")

input.dataset <- "240805_BSimons_filterHT_cluster_renamed"
s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
# DimPlot(object = s.obj, reduction = "INTE_UMAP", label = TRUE, label.box = TRUE, group.by = "seurat_clusters")

cluster.id1 <- 6
to.test.clusters <- setdiff(unique(s.obj$seurat_clusters), c(cluster.id1))

de.markers <- list()
raw.de.markers <- list()

for (cluster.id2 in to.test.clusters){
  tmp <- FindMarkers(object = s.obj, 
                     ident.1 = cluster.id1, 
                     ident.2 = cluster.id2,
                     group.by = "seurat_clusters", 
                     assay = "RNA", 
                     test.use = "wilcox", 
                     min.pct = 0.5)
  tmp.raw <- tmp %>% rownames_to_column("Gene")
  tmp <- tmp %>% rownames_to_column("Gene") %>% subset(p_val_adj <= 0.05) %>%
    rowwise() %>%
    mutate(abs.avg_log2FC = abs(avg_log2FC)) %>%
    arrange(desc(avg_log2FC))
  
  de.markers[[sprintf("%s_vs_%s", cluster.id1, cluster.id2)]] <- tmp
  raw.de.markers[[sprintf("%s_vs_%s", cluster.id1, cluster.id2)]] <- tmp.raw
}

input.df <- raw.de.markers$`6_vs_5` 

input.df <- input.df %>% rowwise() %>%
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-300, p_val_adj))
top10up.genes <- input.df %>% subset(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC)) %>% head(10)
top10down.genes <- input.df %>% subset(p_val_adj <= 0.05) %>% arrange(desc(avg_log2FC)) %>% tail(10)
input.df <- input.df %>%
  rowwise() %>%
  mutate(sig = ifelse(p_val_adj <= 0.05, "Sig.", "not Sig.")) %>%
  mutate(show.gene.name = ifelse(Gene %in% c(top10up.genes$Gene, top10down.genes$Gene), Gene, NA))


cutoff.adjp <- 0.05
volcano.plot <- ggplot(data=input.df, 
                       aes(x=avg_log2FC, 
                           y=-log10(p_val_adj), 
                           col = sig, label=Gene)) + 
  geom_point(size = 3) + geom_label_repel(label = input.df$show.gene.name, size = 8) + 
  scale_color_manual(values=c("#c0d2f0", "#f28095")) +
  geom_vline(xintercept=c(-1, 1), col="#9a9fa6") +
  geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6") +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 12), axis.text=element_text(size=12)) + 
  theme_bw()
volcano.plot

