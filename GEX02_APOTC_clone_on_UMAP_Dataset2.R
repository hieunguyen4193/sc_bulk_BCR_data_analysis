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
clone.name <- "VJcombi_CDR3_0.85"
dataset.name <- "2nd_dataset"
save.dev <- "svg"
topN <- 5
path.to.02.output <- file.path(outdir, 
                               "GEX_output", 
                               sprintf("02_output_Dataset2"), 
                               dataset.name, 
                               clone.name,
                               sprintf("top%s", topN),
                               save.dev)
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

sc.projects.with.ht <- c("240805_BSimons_filterHT_cluster_renamed")
reduction.name <- "INTE_UMAP"

print(sprintf("Working on dataset %s", dataset.name))

s.obj <- readRDS(path.to.all.s.obj[[dataset.name]])

count.yfp.gex <- GetAssayData(object = s.obj, assay = "RNA", slot = "count")["YFP", ]

yfp.cells <- count.yfp.gex[count.yfp.gex != 0] %>% names()
non.yfp.cells <- count.yfp.gex[count.yfp.gex == 0] %>% names()

#####----------------------------------------------------------------------#####
##### PREPARATION
#####----------------------------------------------------------------------#####

clonedf <- table(s.obj$VJcombi_CDR3_0.85) %>% as.data.frame() 
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count)) %>%
  rowwise() %>%
  mutate(Sample_132 = nrow(subset(s.obj@meta.data, s.obj@meta.data$name == "Sample_132" & VJcombi_CDR3_0.85 == clone))) %>%
  mutate(Sample_133 = nrow(subset(s.obj@meta.data, s.obj@meta.data$name == "Sample_133" & VJcombi_CDR3_0.85 == clone)))
  
#####----------------------------------------------------------------------#####
##### PLOTTING TOP CLONES IN EACH MOUSE
#####----------------------------------------------------------------------#####
topN <- 5
top.clones <- list()
all.top.clones <- c()
for (sample.id in unique(s.obj$name)){
  top.clones[[sample.id]] <- head(clonedf[order(clonedf[[sample.id]], decreasing = TRUE), ], topN)$clone %>% as.character()
  all.top.clones <- c(all.top.clones, top.clones[[sample.id]])
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


#####----------------------------------------------------------------------#####
##### PLOTTING YFP CLONES FROM 2 MICE
#####----------------------------------------------------------------------#####
dir.create(file.path(path.to.02.output, "YFP_clones"), showWarnings = FALSE, recursive = TRUE)
clonedf <- clonedf %>% rowwise() %>%
  mutate(clone.size = length(row.names(subset(s.obj@meta.data, s.obj@meta.data$VJcombi_CDR3_0.85 == clone)))) %>%
  mutate(YFP.clone = 
    length(intersect(
      row.names(subset(s.obj@meta.data, s.obj@meta.data$VJcombi_CDR3_0.85 == clone)), yfp.cells
    ))
  ) %>%
  mutate(non.YFP.clone = 
           length(intersect(
             row.names(subset(s.obj@meta.data, s.obj@meta.data$VJcombi_CDR3_0.85 == clone)), non.yfp.cells
           ))
  )

clonedf$inSample <- unlist(lapply(clonedf$clone, function(x){
  tmpdf <- subset(clonedf, clonedf$clone == x)
  if (tmpdf$Sample_132 == 0 & tmpdf$Sample_133 != 0){
    return("Sample_133")
  } else if (tmpdf$Sample_132 != 0 & tmpdf$Sample_133 == 0){
    return("Sample_132")
  } else if (tmpdf$Sample_132 != 0 & tmpdf$Sample_133 != 0){
    return("shared_clone")
  }
}))

for (sample.set in list(c("Sample_132"), 
                        c("Sample_133"), 
                        c("Sample_132", "Sample_133"),
                        c("shared_clone"))){
  yfp.clones <- subset(clonedf, 
                         clonedf$YFP.clone != 0 &
                         clonedf$inSample %in% sample.set) %>% arrange(desc(YFP.clone))
  topN <- 10
  
  plot.clonedf <- data.frame(clone = unique(head(yfp.clones, topN))$clone)
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
         filename = sprintf("APOTC.%s.%s", paste(sample.set, collapse = "_"), save.dev), 
         path = file.path(path.to.02.output, "YFP_clones"), dpi = 300, width = 14, height = 10)    
  
}

##### plot APOTC with color by sample Sample_132 or Sample_133
sample.set <- c("Sample_132", "Sample_133")
colors <- tableau_color_pal(palette = "Tableau 10")(2)
names(colors) <- unique(s.obj$name)

yfp.clones <- subset(clonedf, 
                     clonedf$YFP.clone != 0 &
                       clonedf$inSample %in% sample.set) %>% arrange(desc(YFP.clone)) %>%
  rowwise() %>%
  mutate(color = colors[[inSample]])

plot.clonedf <- head(yfp.clones, topN) %>% subset( select = c(clone, color, inSample)) 

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
       filename = sprintf("APOTC.%s.color_by_mouseID.%s", paste(sample.set, collapse = "_"), save.dev), 
       path = path.to.02.output, dpi = 300, width = 14, height = 10) 
