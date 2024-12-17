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

my_comparisons <- list(
  c("3w_d0", "3w_d90"),
  c("3w_d0", "8w_d0"),
  c("3w_d0", "8w_d90"),
  c("3w_d0", "9w_d90")
)
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), 
                    symbols = c("****", "***", "**", "*", "ns"))

p1 <- nodedf %>% subset(population == "biopsy") %>% 
  ggplot(aes(x = age_day, y = dist_to_deepest)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to deepest node") +
  stat_compare_means(comparisons = my_comparisons)

p2 <- nodedf %>% subset(population == "biopsy") %>% 
  ggplot(aes(x = age_day, y = dist_to_root)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", 
                     symnum.args = symnum.args)


p1 <- nodedf %>%  
  ggplot(aes(x = age_day, y = dist_to_deepest)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to deepest node") +
  stat_compare_means(comparisons = my_comparisons)

p2 <- nodedf %>% 
  ggplot(aes(x = age_day, y = dist_to_root)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", 
                     symnum.args = symnum.args)


nodedf %>% subset(cloneID == "m11_IGHV1-78-01_IGHJ3-01_24_3.aln") %>% ggplot(aes(x = dist_to_root, y = dist_to_deepest)) + geom_point()
p1 + p2


p1 <- nodedf  %>% subset(mixed_node == "no") %>%  
  ggplot(aes(x = population, y = dist_to_deepest)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to deepest node") +
  stat_compare_means(comparisons = my_comparisons)

p2 <- nodedf %>% subset(mixed_node == "no") %>%
  ggplot(aes(x = population, y = dist_to_root)) + 
  geom_boxplot() + 
  theme_pubr() + 
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("Distance to root node (GL)") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test", 
                     symnum.args = symnum.args)
p1 + p2

summarydf <- table(nodedf$cloneID) %>% data.frame()
colnames(summarydf) <- c("cloneID", "count")
summarydf <- summarydf %>% rowwise() %>%
  mutate(mouseid= str_split(cloneID, "_")[[1]][[1]])%>%
  mutate(age = unique(subset(mid.metadata, mid.metadata$mouse == mouseid)$age)[[1]]) %>%
  mutate(day = unique(subset(mid.metadata, mid.metadata$mouse == mouseid)$day)[[1]]) %>%
  mutate(age_day = sprintf("%s_%s", age, day))



table(subset(summarydf, summarydf$count >= 50)$age_day)

