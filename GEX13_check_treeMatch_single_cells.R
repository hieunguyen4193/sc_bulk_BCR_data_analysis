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

# input.dataset <- "240805_BSimons_filterHT_cluster_renamed"
input.dataset <- "241002_241104_BSimons"

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

path.to.13.output <- file.path(outdir, "GEX_output", "13_output", input.dataset)
dir.create(path.to.13.output, showWarnings = FALSE, recursive = TRUE)

sc.dataset.with.hto <- c("240805_BSimons",
                         "241104_BSimons",
                         "240805_BSimons_filterHT",
                         "240805_BSimons_filterHT_cluster_renamed",
                         "241002_BSimons",
                         "240805_BSimons_filterHT_cluster",
                         "241002_241104_BSimons")



s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])

if (input.dataset == "240805_BSimons_filterHT_cluster_renamed"){
  match.scdf <- readxl::read_excel(file.path(outdir, "tree_analysis",
  "240826_BSimons", "04_output", 
  "240826_BSimons_240805_BSimons_filterHT_cluster_renamed",
  "240826_BSimons_match_barcode_scdf.xlsx"))
} else if (input.dataset == "241002_241104_BSimons"){
  match.scdf <- readxl::read_excel(file.path(outdir, "tree_analysis",
                                      "241031_BSimons", "04_output", 
                                      "241031_BSimons_241002_241104_BSimons",
                                      "241031_BSimons_match_barcode_scdf.xlsx")) 
}
match.scdf <- subset(match.scdf, select = -c(`...1`))

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("barcode") %>%
  rowwise() %>%
  mutate(matchTree = ifelse(barcode %in% match.scdf$barcode, "match", "no"))
  
p <- DimPlot(object = s.obj, reduction = "INTE_UMAP", 
        cells.highlight = subset(meta.data, meta.data$matchTree == "match")$barcode, 
        label = TRUE, label.box = TRUE, repel = TRUE, group.by = "seurat_clusters") + 
  theme(legend.position = "none") +
  ggtitle("UMAP highlighting cells which match tree sequences")
ggsave(plot = p, filename = "UMAP_highlight_matchTree_cells.svg",
       path = path.to.13.output, device = "svg", width = 14, height = 10, dpi = 300)

countdf <- table(meta.data$seurat_clusters, meta.data$matchTree) %>%
  data.frame() %>% pivot_wider(names_from = "Var1", values_from = "Freq") %>% data.frame()
colnames(countdf) <- to_vec(for (i in colnames(countdf)){
  if (i == "Var2"){
    "Group"
  } else {
    str_replace(i, "X", "cluster_")
  }
})
writexl::write_xlsx(countdf, file.path(path.to.13.output, "count_match_cells_in_clusters.xlsx"))
  
countdf <- countdf %>% column_to_rownames("Group")
writexl::write_xlsx(countdf/rowSums(countdf), file.path(path.to.13.output, "pct_match_cells_in_clusters.xlsx"))