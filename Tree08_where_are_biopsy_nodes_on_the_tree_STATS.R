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
mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep = ";")
nodedf <- read.csv(file.path(path.to.08.output, "nodedf.csv"))
nodedf <- merge(nodedf, subset(mid.metadata, select = c(mouse, age, day)), 
                by.x = "mouseid", by.y = "mouse", all.x = TRUE) %>%
  rowwise() %>%
  mutate(age_day = sprintf("%s_%s", age, day))

mouse.id <- "m31"

sub.nodedf <- subset(nodedf, nodedf$mouseid == mouse.id)

sub.nodedf %>% ggplot(aes(x = population, y = dist_to_root)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)")

sub.nodedf %>% ggplot(aes(x = population, y = dist_to_deepest)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to deepest node")

nodedf %>% subset(population == "biopsy") %>% 
  ggplot(aes(x = age, y = dist_to_deepest)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to deepest node")

nodedf %>% subset(population == "biopsy") %>% 
  ggplot(aes(x = day, y = dist_to_root)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)")

nodedf %>% subset(population == "biopsy") %>% 
  ggplot(aes(x = age, y = dist_to_root)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)")

nodedf %>% subset(population == "biopsy") %>% 
  ggplot(aes(x = age_day, y = dist_to_root)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)")

nodedf %>% subset(population == "biopsy") %>% 
  ggplot(aes(x = age_day, y = dist_to_deepest)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)")

