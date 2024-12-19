gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
library(circlize)
library(ggpubr)
#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

PROJECT <- "220701_etc_biopsies"
path.to.main.output <- file.path(outdir, "tree_analysis", PROJECT)
path.to.08.output <- file.path(path.to.main.output, "08_output")

path.to.save.figures <- file.path(path.to.08.output, "figures_nodedf")
dir.create(path.to.save.figures, showWarnings = FALSE, recursive = TRUE)

mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep = ";")
if (file.exists(file.path(path.to.08.output, "nodedf.addMetadata.csv")) == FALSE){
  nodedf <- read.csv(file.path(path.to.08.output, "nodedf.csv")) %>%
    rowwise() %>%
    mutate(age = unique(subset(mid.metadata, mid.metadata$mouse == mouseid)$age)[[1]]) %>%
    mutate(day = unique(subset(mid.metadata, mid.metadata$mouse == mouseid)$day)[[1]]) %>%
    mutate(age_day = sprintf("%s_%s", age, day))
  write.csv(nodedf, file.path(path.to.08.output, "nodedf.addMetadata.csv"))  
} else {
  nodedf <- read.csv(file.path(path.to.08.output, "nodedf.addMetadata.csv"))
}

my_comparisons <- list(
  c("3w_d0", "3w_d90"),
  c("3w_d0", "8w_d0"),
  c("3w_d0", "8w_d90")
)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                    symbols = c("****", "***", "**", "*", "ns"))

summarydf <- table(nodedf$cloneID) %>% data.frame()
colnames(summarydf) <- c("cloneID", "count")
summarydf <- summarydf %>% rowwise() %>%
  mutate(mouseid= str_split(cloneID, "_")[[1]][[1]])%>%
  mutate(age = unique(subset(mid.metadata, mid.metadata$mouse == mouseid)$age)[[1]]) %>%
  mutate(day = unique(subset(mid.metadata, mid.metadata$mouse == mouseid)$day)[[1]]) %>%
  mutate(age_day = sprintf("%s_%s", age, day)) %>% arrange(desc(count))
write.csv(summarydf, file.path(path.to.08.output, "count_node_in_trees.csv"))

count.tree.plots <- summarydf %>% ggplot(aes(x = age_day, y = count)) + 
  geom_boxplot() +
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90, size = 25),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 25)) +
  stat_compare_means(method = "t.test", comparisons = my_comparisons, symnum.args = symnum.args)
ggsave(plot = count.tree.plots, filename = sprintf("number_of_trees_each_mouse_age_day.svg"), 
       path = path.to.save.figures, dpi = 300, width = 14, height = 10)

plotdf <- nodedf[, c("mouseid", "seqid", "cloneID", "population", "abundance", 
                     "topo_dist_to_root", "dist_to_root", 
                     "topo_dist_to_furthest_child_node", "dist_to_furthest_child_node", 
                     "age", "day", "age_day")]
colnames(plotdf) <- c("mouseid", "seqid", "cloneID", "population", "abundance",
                      "topo_rootness", "rootness", 
                      "topo_leafness", "leafness",
                      "age", "day", "age_day")

#####----------------------------------------------------------------------#####
##### All nodes
#####----------------------------------------------------------------------#####
for (input.feat in c("rootness", "topo_rootness", "leafness", "topo_leafness")){
  dir.create(file.path(path.to.save.figures, "all_nodes"), showWarnings = FALSE, recursive = TRUE)
  p <- plotdf %>% ggplot(aes_string(x = "age_day", y = input.feat)) + 
    geom_boxplot() + 
    theme_pubr() + 
    theme(axis.text.x = element_text(angle = 90, size = 25),
          axis.text.y = element_text(size = 25),
          axis.title = element_text(size = 25)) +
    stat_compare_means(method = "t.test", comparisons = my_comparisons, symnum.args = symnum.args)
  ggsave(plot = p, filename = sprintf("boxplot_%s_all_nodes.svg", input.feat), 
         path = file.path(path.to.save.figures, "all_nodes"), dpi = 300, width = 14, height = 10)
}

#####----------------------------------------------------------------------#####
##### filter node
#####----------------------------------------------------------------------#####
for (node.type in unique(plotdf$population)){
  dir.create(file.path(path.to.save.figures, node.type), showWarnings = FALSE, recursive = TRUE)
  for (input.feat in c("rootness", "topo_rootness", "leafness", "topo_leafness")){
    p <- plotdf %>% 
      subset(population == node.type) %>% 
      ggplot(aes_string(x = "age_day", y = input.feat)) + 
      geom_boxplot() + 
      theme_pubr() + 
      theme(axis.text.x = element_text(angle = 90, size = 25),
            axis.text.y = element_text(size = 25),
            axis.title = element_text(size = 25)) +
      stat_compare_means(method = "t.test", comparisons = my_comparisons, symnum.args = symnum.args)
    ggsave(plot = p, filename = sprintf("boxplot_%s_%s.svg", input.feat, node.type), 
           path = file.path(path.to.save.figures, node.type), dpi = 300, width = 14, height = 10)
  }
}
