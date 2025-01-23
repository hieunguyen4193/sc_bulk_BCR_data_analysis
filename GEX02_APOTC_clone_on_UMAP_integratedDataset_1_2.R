gc()
rm(list = ls())
# install.packages("dplyr")
# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-1.tar.gz", type = "source", repos = NULL)
new.pkgs <- c("APackOfTheClones", "svglite", "car", "ggpubr", "ggthemes", "dplyr")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }
}
library(ggpubr)
library(ggthemes)
library(APackOfTheClones)
library("gridExtra")

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.all.s.obj <- path.to.all.s.obj[setdiff(names(path.to.all.s.obj), c("BonnData"))]

sc.projects.with.ht <- c("240805_BSimons_filterHT_cluster_renamed")

# name_or_colonization <- "colonization"
name_or_colonization <- "name"

if (name_or_colonization == "sample_ht"){
  path.to.all.s.obj <- path.to.all.s.obj[sc.projects.with.ht]
}
clone.name <- "VJcombi_CDR3_0.85"
dataset.name <- "Dataset1_2"
save.dev <- "svg"
topN <- 5
print(sprintf("Working on dataset %s", dataset.name))

s.obj <- readRDS(path.to.all.s.obj[[dataset.name]])

if (dataset.name %in% sc.projects.with.ht){
  meta.data <- s.obj@meta.data %>% 
    rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(sample_ht = sprintf("%s_%s", name, HTO_classification)) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$sample_ht, col.name = "sample_ht")      
}

DefaultAssay(s.obj) <- "RNA"
path.to.02.output <- file.path(outdir, 
                               "GEX_output", 
                               sprintf("02_output_Dataset1_2_%s", name_or_colonization), 
                               dataset.name, 
                               clone.name,
                               sprintf("top%s", topN),
                               save.dev)
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

if (dataset.name %in% c("241104_BSimons", "241002_BSimons")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}
Idents(s.obj) <- "seurat_clusters"
s.obj <- RunAPOTC(seurat_obj = s.obj, reduction_base = reduction.name, clonecall = clone.name)

colonizations <- list( MM9_Ecoli = c("17_MM9_Ecoli", 
                                     "20_MM9_Ecoli"),
                       MM9_Ecoli_SPF = c("21_MM9_Ecoli_SPF"),
                       MM9 = c("MM9_S2",
                                   "MM9_S4"),
                       MM9_SPF = c( "MM9_SPF_S3",
                                    "MM9_SPF_S9"),
                       dataset2 = c("Sample_132",
                                    "Sample_133")
                       )
colonizationdf <- data.frame(SampleID = unlist(colonizations)) %>%
  rownames_to_column("colonization") %>%
  rowwise() %>%
  mutate(colonization = ifelse(colonization != "MM9_Ecoli_SPF", 
                               paste(str_split(colonization, "")[[1]][1:nchar(colonization)- 1], collapse = ""), 
                               "MM9_Ecoli_SPF"))

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") %>%
  rowwise() %>%
  mutate(colonization = subset(colonizationdf, colonizationdf$SampleID == name)$colonization) %>%
  column_to_rownames("barcode")

meta.data <- meta.data[row.names(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, col.name = "colonization", metadata = meta.data$colonization)

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("cell.barcode")
clonedf <- data.frame(meta.data[[clone.name]] %>% table())
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count))

top.clones <- list()
all.top.clones <- c()
for (sampleid in unique(s.obj@meta.data[[name_or_colonization]])){
  print(sprintf("Generating clonedf, working on sample: %s", sampleid))
  clonedf[[sampleid]] <- unlist(lapply(clonedf$clone, function(x){
    nrow(subset(meta.data, meta.data[[name_or_colonization]] == sampleid & meta.data[[clone.name]] == x))
  }))
  top.clones[[sampleid]] <- clonedf[order(clonedf[[sampleid]], decreasing = TRUE),] %>% head(topN) %>% pull(clone) %>% as.character()
  all.top.clones <- c(all.top.clones, top.clones[[sampleid]])
}

plot.clonedf <- data.frame(clone = unique(all.top.clones))
if (length(plot.clonedf$clone) <= 20){
  colors <- tableau_color_pal(palette = "Tableau 20")(length(plot.clonedf$clone))  
} else {
  colors <- c(tableau_color_pal(palette = "Tableau 20")(20), hue_pal()(length(plot.clonedf$clone) - 20))
}

plot.clonedf$color <- colors

tmp.plot <- vizAPOTC(s.obj, clonecall = clone.name, 
                     verbose = FALSE, 
                     reduction_base = reduction.name, 
                     show_labels = TRUE, 
                     repulsion_strength = 5,
                     legend_position = "top_right", 
                     legend_sizes = 2) %>% 
  showCloneHighlight(clonotype =  as.character(plot.clonedf$clone), 
                     fill_legend = TRUE,
                     color_each = plot.clonedf$color, 
                     default_color = "lightgray") 

ggsave(plot = tmp.plot, 
       filename = sprintf("APOTC.%s", save.dev), 
       path = path.to.02.output, dpi = 300, width = 14, height = 10)    

for (c in names(colonizations)){
  p.umap <- DimPlot(object = subset(s.obj, colonization == c), reduction = reduction.name, label = TRUE, label.box = TRUE)
  tmp.plot <- vizAPOTC(subset(s.obj, colonization == c), clonecall = clone.name, 
                       verbose = FALSE, 
                       reduction_base = reduction.name, 
                       show_labels = TRUE, 
                       repulsion_strength = 5,
                       legend_position = "top_right", 
                       legend_sizes = 2) %>% 
  showCloneHighlight(clonotype =  as.character(top.clones[[c]]), 
                       fill_legend = TRUE,
                       color_each = subset(plot.clonedf, plot.clonedf$clone %in% top.clones[[c]])$color, 
                       default_color = "lightgray") 
  
  ggsave(plot = tmp.plot, 
         filename = sprintf("APOT_%s.%s", c, save.dev), 
         path = path.to.02.output, dpi = 300, width = 14, height = 10)    
  ggsave(plot = p.umap, 
         filename = sprintf("UMAP_colonization_%s.%s", c, save.dev), 
         path = path.to.02.output, dpi = 300, width = 14, height = 10)    
}

