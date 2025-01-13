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
#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
thres <- 0.85
thres.dis <- 0.15

PROJECT <- "241002_BSimons"
path.to.VDJ.output <- file.path( outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

#####----------------------------------------------------------------------#####
##### READ CLONE DATA
#####----------------------------------------------------------------------#####
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", PROJECT)
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)

clone.obj <- run_preprocessing_all_sc_data(path.to.VDJ.output = path.to.VDJ.output, 
                                                   path.to.save.output = path.to.save.output, 
                                                   PROJECT = PROJECT,
                                                   thres = thres, 
                                                   thres.dis = thres.dis,
                                                   savefile = TRUE,
                                                   rerun = FALSE,
                                                   define.clone.clusters =  FALSE) 

#####----------------------------------------------------------------------#####
##### READ GENE EXPRESSION DATA
#####----------------------------------------------------------------------#####
s.obj <- readRDS(path.to.all.s.obj[[PROJECT]])
meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")

if (PROJECT %in% c("241002_BSimons", "241104_BSimons")){
  yfp.exprs <- GetAssayData(object = s.obj, assay = "RNA", slot = "data")["YFP", ]
  all.cells <- colnames(s.obj)
  yfp.cells <- yfp.exprs[yfp.exprs != 0] %>% names()
  non.yfp.cells <- setdiff(all.cells, yfp.cells)
  meta.data <- meta.data %>% rowwise() %>%
    mutate(YFP = ifelse(barcode %in% yfp.cells, "YES", "NO"))
  
  to.add.metadata <- meta.data %>% column_to_rownames("barcode")
  to.add.metadata <- to.add.metadata[row.names(s.obj@meta.data),]
  s.obj <- AddMetaData(object = s.obj, col.name = "YFP", metadata = to.add.metadata$YFP)
  if ("INTE_UMAP" %in% names(s.obj@reductions)){
    yfp.plot <- DimPlot(object = s.obj, reduction = "INTE_UMAP", cells.highlight = yfp.cells, cols.highlight = "red") +
      theme(legend.position = "none")
  } else {
    yfp.plot <- DimPlot(object = s.obj, reduction = "RNA_UMAP", cells.highlight = yfp.cells, cols.highlight = "red") +
      theme(legend.position = "none")
  }
}

clonesets <- clone.obj$clonesets %>% rowwise() %>%
  mutate(barcode_full = sprintf("%s_%s", id, barcode)) %>%
  mutate(in_GEX_data = ifelse(barcode_full %in% colnames(s.obj), "YES", "NO")) 

clonesets.filtered <- clonesets %>%
  subset(in_GEX_data == "YES") %>%
  mutate(ht = subset(meta.data, meta.data$barcode == barcode_full)$HTO_classification) %>%
  mutate(id_hashtag = sprintf("%s_%s", id, ht))

new.clonesets <- data.frame()
for (input.VJ.combi in unique(clonesets.filtered$VJ.len.combi)){
  tmpdf <- subset(clonesets.filtered, clonesets.filtered$VJ.len.combi == input.VJ.combi)
  seqs <- unique(tmpdf$aaSeqCDR3)
  print(sprintf("VJ.len.combi: %s, num seqs: %s", input.VJ.combi, length(seqs)))
  if (length(seqs) >= 2){
    cluster.output <- assign_clusters_to_sequences(seqs = seqs, threshold = thres.dis)$res
    tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
      tmpdf$aaSeqCDR3, function(x){
        return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
      }
    ))    
  } else {
    tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
  }
  new.clonesets <- rbind(new.clonesets, tmpdf)
}

clonedf <- table(new.clonesets$VJcombi_CDR3_0.85) %>% data.frame()
colnames(clonedf) <- c("clone", "count")
clonedf <- clonedf %>% arrange(desc(count))
for (sample.id in sort(unique(clonesets.filtered$id_hashtag))){
  clonedf[[sample.id]] <- unlist(lapply(
    clonedf$clone, function(x){
      subset(new.clonesets, new.clonesets[[sprintf("VJcombi_CDR3_%s", thres)]] == x &
               new.clonesets$id_hashtag == sample.id) %>% nrow()
    }
  ))
}

#####----------------------------------------------------------------------#####
##### PLOT CIRCOS 
#####----------------------------------------------------------------------#####
write.csv(new.clonesets, file.path(path.to.05.output, sprintf("%s.clonesets.xlsx", PROJECT)))
source(file.path(path.to.main.src, "VDJ_generate_circos_plots.R"))
dir.create(file.path(path.to.05.output, "circos_plot"), showWarnings = FALSE, recursive = TRUE)
generate_circos_plot(input.clonesets = new.clonesets, 
                     path.to.save.svg = file.path(path.to.05.output, "circos_plot"),
                     svg.name = sprintf("%s_all_samples_hashtags.svg", PROJECT))
