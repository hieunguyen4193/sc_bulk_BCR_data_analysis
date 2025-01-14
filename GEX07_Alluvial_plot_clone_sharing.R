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
# install.packages("https://cran.r-project.org/src/contrib/ggplot2_3.5.1.tar.gz", repos = NULL, type = "source")

library(ggalluvial)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library("gridExtra")

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
# source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

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

name_or_sampleHT <- "name" # sample_ht
clone.name <- "VJcombi_CDR3_0.85"
dataset.name <- "240805_BSimons_filterHT_cluster_renamed"
save.dev <- "svg"
num.top <- 5

if (name_or_sampleHT == "sample_ht"){
  path.to.all.s.obj <- path.to.all.s.obj[sc.projects.with.ht]
}
path.to.07.output <- file.path(outdir, 
                               "GEX_output", 
                               sprintf("07_output_%s", name_or_sampleHT), 
                               dataset.name, 
                               clone.name,
                               sprintf("top%s", num.top),
                               save.dev
                               )
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

##### read seurat object 
s.obj <- readRDS(path.to.all.s.obj[[dataset.name]])

if (dataset.name %in% sc.projects.with.ht){
  meta.data <- s.obj@meta.data %>% 
    rownames_to_column("barcode") %>%
    rowwise() %>%
    mutate(sample_ht = sprintf("%s_%s", name, HTO_classification)) %>%
    mutate(mouseid = sprintf("m%s", str_replace(name, "M", ""))) %>%
    mutate(mouseid = sprintf("%s", str_replace(mouseid, "P", ""))) %>%
    mutate(mouseid = sprintf("%s", str_replace(mouseid, "PP", ""))) %>%
    column_to_rownames("barcode")
  meta.data <- meta.data[row.names(s.obj@meta.data), ]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$sample_ht, col.name = "sample_ht")    
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$mouseid, col.name = "mouseid")    
}

DefaultAssay(s.obj) <- "RNA"

#####----------------------------------------------------------------------#####
##### Generate clonedf
#####----------------------------------------------------------------------#####
meta.data <- s.obj@meta.data %>% 
  rownames_to_column("cell.barcode")
clonedf <- data.frame(meta.data[[clone.name]] %>% table())
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count))

selected.top.clones <- list()
all.selected.top.clones <- c()
for (sampleid in unique(s.obj@meta.data[[name_or_sampleHT]])){
  print(sprintf("Generating clonedf, working on sample: %s", sampleid))
  clonedf[[sampleid]] <- unlist(lapply(clonedf$clone, function(x){
    nrow(subset(meta.data, meta.data[[name_or_sampleHT]] == sampleid & meta.data[[clone.name]] == x))
  }))
  selected.top.clones[[sampleid]] <- clonedf[order(clonedf[[sampleid]], decreasing = TRUE),] %>% 
    head(num.top) %>% 
    pull(clone) %>%
    as.character()
  all.selected.top.clones <- c(all.selected.top.clones, selected.top.clones[[sampleid]])
}

#####----------------------------------------------------------------------#####
##### prepare the data frame for plotting
#####----------------------------------------------------------------------#####
for (mouse.id in unique(s.obj$mouseid)){
  print(sprintf("Working on mouse ID: %s", mouse.id))
  all.plot.samples <- unique(subset(s.obj@meta.data, s.obj@meta.data$mouseid == mouse.id)[[name_or_sampleHT]]) %>% 
    sort()
  
  p.sample <- all.plot.samples[grepl("P", all.plot.samples)]
  m.sample <- all.plot.samples[grepl("M", all.plot.samples)]
  
  plotdf <- clonedf[, c("clone", all.plot.samples)] %>% 
    subset(clone %in% unlist(selected.top.clones[all.plot.samples])) 
  
  for (sampleid in all.plot.samples){
    plotdf[[sampleid]] <- to_vec(for (item in plotdf[[sampleid]]){
      item / sum(plotdf[[sampleid]])
    })  
  }
  
  row.names(plotdf) <- NULL
  plotdf$check.share <- plotdf %>% column_to_rownames("clone") %>% as.matrix() %>% rowProds()
  
  tmp.plotdf <- plotdf %>% column_to_rownames("clone")
  share.clones <- tmp.plotdf[rowSums(tmp.plotdf != 0) >= 2, ] %>% row.names() %>% as.character()
  plotdf <- subset(plotdf, select = -c(check.share))
  
  ##### plot alluvial plots
  plotdf.pivot <- plotdf %>% pivot_longer(!clone, names_to = "SampleID", values_to = "Count")
  colnames(plotdf.pivot) <- c("Clone", "SampleID", "Count", "CloneID")
  plotdf.pivot$Clone <- as.character(plotdf.pivot$Clone)
  plotdf.pivot <- plotdf.pivot %>% rowwise() %>%
    mutate(CloneID = ifelse(Clone %in% as.character(share.clones), Clone, NA))
  
  # is_lodes_form(plotdf.pivot, key = SampleID, value = Count, id = clone, silent = TRUE)
  unique.clones <- unique(plotdf.pivot$Clone)
  if (length(unique.clones) <= 20){
    colors <- tableau_color_pal(palette = "Tableau 20")(length(unique.clones))  
  } else {
    colors <- c(tableau_color_pal(palette = "Tableau 20")(20), hue_pal()(length(unique.clones) - 20))
  }
  names(colors) <- unique.clones
  
  if (name_or_sampleHT == "name"){
    clone.orders <- plotdf[order(plotdf[[p.sample]]), ]$clone %>% as.character()
    colors <- colors[clone.orders]
    plotdf.pivot$Clone <- factor(plotdf.pivot$Clone, levels = clone.orders)  
  }
  
  plotdf.pivot$SampleID <- factor(plotdf.pivot$SampleID, levels = c(p.sample, m.sample))
  
  allu.plot <- ggplot(plotdf.pivot,
                      aes(x = SampleID, stratum = Clone, alluvium = Clone, y = Count, fill = Clone)) +
    geom_alluvium(aes(fill = CloneID), width=.5) +
    geom_stratum(width=.5) +
    # geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_pubr() + 
    scale_fill_manual(values = colors, na.value = "lightgray") + 
    theme(legend.position = "bottom")
  ggsave(plot = allu.plot, filename = sprintf("alluvial_plot_%s.%s", mouse.id, save.dev),
         path = file.path(path.to.07.output), device = save.dev, width = 14, height = 10, dpi = 300)
}

