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
# name_or_sampleHT <- "sample_ht"
name_or_sampleHT <- "name"

if (name_or_sampleHT == "sample_ht"){
  path.to.all.s.obj <- path.to.all.s.obj[sc.projects.with.ht]
}
clone.name <- "VJcombi_CDR3_0.85"
dataset.name <- "240805_BSimons_filterHT_cluster_renamed"
save.dev <- "png"

##### sample.list for each mouse:
if (grepl("240805_BSimons", dataset.name) == TRUE){
  sample.list <- list(
    M_samples = c("M1", "M2", "M3"),
    P_samples = c("P1", "P2", "P3")
  )
  mouse.list <- list(
    m1 = c("M1", "P1"),
    m2 = c("M2", "P2"),
    m3 = c("M3", "P3")
  )
} else if (grepl("241002_BSimons", dataset.name) == TRUE){
  sample.list <- list(
    m3 = c("PP3")
  )
} else if (grpel("241104_BSimons", dataset.name) == TRUE){
  sample.list <- list(
    m7 = c("PP7")
  )
}
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
path.to.07.output <- file.path(outdir, 
                               "GEX_output", 
                               sprintf("07_output"), 
                               dataset.name, 
                               clone.name)
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)
if (dataset.name %in% c("241104_BSimons", "241002_BSimons")){
  reduction.name <- "RNA_UMAP"
} else {
  reduction.name <- "INTE_UMAP"
}

Idents(s.obj) <- "seurat_clusters"

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("cell.barcode")
clonedf <- data.frame(meta.data[[clone.name]] %>% table())
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count))

for (sampleid in unique(s.obj@meta.data[[name_or_sampleHT]])){
  print(sprintf("Generating clonedf, working on sample: %s", sampleid))
  clonedf[[sampleid]] <- unlist(lapply(clonedf$clone, function(x){
    nrow(subset(meta.data, meta.data[[name_or_sampleHT]] == sampleid & meta.data[[clone.name]] == x))
  }))
}

#####----------------------------------------------------------------------#####
##### Clone sharing between samples
#####----------------------------------------------------------------------#####
mouse.id <- "m1"

selected.clones <- c()
top.clones <- list()
for (sample.id in mouse.list[[mouse.id]]){
  tmpdf <- clonedf[, c("clone", sample.id)]
  colnames(tmpdf) <- c("clone", "SampleID")
  tmpdf <- tmpdf %>% arrange(desc(SampleID)) %>% head(5)
  selected.clones <- c(selected.clones, as.character(tmpdf$clone))
  top.clones[[sample.id]] <- as.character(tmpdf$clone)
}

plotdf <- clonedf[, c("clone", mouse.list[[mouse.id]])] %>% subset(clone %in% selected.clones)

##### scale to percentage
for (sample.id in mouse.list[[mouse.id]]){
  plotdf[[sample.id]] <- to_vec(for (item in plotdf[[sample.id]]){
    item/sum(plotdf[[sample.id]])
  })
}

p.sample <- mouse.list[[mouse.id]][grepl("P", mouse.list[[mouse.id]])]
m.sample <- mouse.list[[mouse.id]][grepl("M", mouse.list[[mouse.id]])]

share.clones <- subset(plotdf, 
                         plotdf[[colnames(plotdf)[[2]]]] != 0 &
                         plotdf[[colnames(plotdf)[[3]]]] != 0 )$clone

plotdf.pivot <- plotdf %>% pivot_longer(!clone, names_to = "SampleID", values_to = "Count")
colnames(plotdf.pivot) <- c("Clone", "SampleID", "Count", "CloneID")
plotdf.pivot$Clone <- as.character(plotdf.pivot$Clone)
plotdf.pivot <- plotdf.pivot %>% rowwise() %>%
  mutate(CloneID = ifelse(Clone %in% as.character(share.clones), Clone, NA))

# is_lodes_form(plotdf.pivot, key = SampleID, value = Count, id = clone, silent = TRUE)
clone.orders <- plotdf[order(plotdf[[p.sample]]), ]$clone %>% as.character()
unique.clones <- unique(plotdf.pivot$Clone)
colors <- tableau_color_pal(palette = "Tableau 20")(length(unique.clones))
names(colors) <- unique.clones

colors <- colors[clone.orders]

plotdf.pivot$Clone <- factor(plotdf.pivot$Clone, levels = clone.orders)
plotdf.pivot$SampleID <- factor(plotdf.pivot$SampleID, levels = c(p.sample, m.sample))

ggplot(plotdf.pivot,
       aes(x = SampleID, stratum = Clone, alluvium = Clone, y = Count, fill = Clone)) +
  geom_alluvium(aes(fill = CloneID), width=.5) +
  geom_stratum(width=.5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_pubr() + 
  scale_fill_manual(values = colors, na.value = "lightgray") + 
  theme(legend.position = "bottom")
