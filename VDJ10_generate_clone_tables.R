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

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

input.dataset <- "Dataset1_2"
for (input.dataset in names(path.to.all.s.obj)){
  path.to.10.output <- path.to.05.output <- file.path(outdir, "VDJ_output", "10_output", input.dataset)
  dir.create(path.to.10.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
  
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
  
  clonedf <- table(meta.data$VJcombi_CDR3_0.85) %>% data.frame()
  colnames(clonedf) <- c("clone", "Freq")
  clonedf <- clonedf %>% rowwise() %>%
    mutate(V.gene = str_split(clone, "_")[[1]][[1]]) %>%
    mutate(J.gene = str_split(clone, "_")[[1]][[2]]) %>%
    mutate(CDR3_seqs = paste(
      subset(meta.data, meta.data$VJcombi_CDR3_0.85 == clone)$aaSeqCDR3 %>% unique(), collapse = "_"
    ))
  
  writexl::write_xlsx(clonedf, file.path(path.to.10.output, sprintf("%s.cloneInfo.xlsx", input.dataset)))
  
  count.celldf <- table(meta.data$seurat_clusters, meta.data$name) %>% data.frame() 
  colnames(count.celldf) <- c("seurat_clusters", "SampleID", "Freq")
  count.celldf <- count.celldf %>% pivot_wider(names_from = "seurat_clusters", values_from = "Freq")
  writexl::write_xlsx(count.celldf, file.path(path.to.10.output, sprintf("count_celldf_%s.xlsx", input.dataset)))
}


  
