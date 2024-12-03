gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
# source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

asssign_mutation_to_group <- function(x){
  if (is.na(x) == TRUE){
    return(NA)
  }
  if (x <= 3){
    return("low")
  } else if (x <= 15){
    return("intermediate")
  } else if (x > 15){
    return("high")
  } 
}

sample.list <- list(
  `1st_dataset` = c("17_MM9_Ecoli",
                    "20_MM9_Ecoli",
                    "21_MM9_Ecoli_SPF",
                    "MM9_S2",
                    "MM9_S4",
                    "MM9_SPF_S3",      
                    "MM9_SPF_S9"),
  `2nd_dataset` = c("Sample_132", 
                    "Sample_133")
)
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

for (input.dataset in names(path.to.all.s.obj)){
  path.to.04.output <- file.path(outdir, "VDJ_output", "04_output")
  
  path.to.06.output <- file.path(outdir, "VDJ_output", "06_output", input.dataset)
  s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
  dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
  if (input.dataset %in% c("1st_dataset", "2nd_dataset")){
    clonedf <- read.csv(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"), row.names = "X") %>%
      subset(id == sample.list[[input.dataset]]) %>%
      rowwise() %>%
      mutate(group_mutation = asssign_mutation_to_group(num_mutation)) %>%
      mutate(barcode_full = sprintf("%s_%s", id, barcode))
  } else {
    clonedf <- read.csv(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"), row.names = "X") %>%
      subset(dataset.name == input.dataset) %>%
      rowwise() %>%
      mutate(group_mutation = asssign_mutation_to_group(num_mutation)) %>%
      mutate(barcode_full = sprintf("%s_%s", id, barcode))
  }
  
  thres <- 0.85
  
  #####---------------------------------------------------------------------------#####
  ##### GENERATE A BIG DATAFRAME CONTAINING ALL CLONES FROM ALL DATASETS
  #####---------------------------------------------------------------------------#####
  all.clone.files <- Sys.glob(file.path(outdir, 
                                        "VDJ_output", 
                                        "*", 
                                        sprintf("VDJ_output_%s", thres), 
                                        "preprocessed_files", 
                                        "clonesets*.split_clones.xlsx" ))
  
  names(all.clone.files) <- to_vec(
    for (item in all.clone.files) str_replace(str_replace(basename(item), "clonesets_", ""), ".split_clones.xlsx", "")
  )
  
  
  meta.data <- s.obj@meta.data %>% 
    rownames_to_column("input_barcode")
  if (length(intersect(meta.data$input_barcode, clonedf$barcode_full)) == 0){
    print("rename barcode")
    meta.data <- meta.data %>%
      rowwise() %>%
      mutate(input_barcode = str_replace(input_barcode, sprintf("%s_%s", name, name), name))
  }
  meta.data <- meta.data %>%
    rowwise() %>%
    mutate(num_mutation = ifelse(
      nrow(subset(clonedf, clonedf$barcode_full == input_barcode)) == 1.0,
      subset(clonedf, clonedf$barcode_full == input_barcode)$num_mutation %>% as.numeric(),
      NA
    )) %>%
    mutate(group_mutation = asssign_mutation_to_group(num_mutation)) %>%
    column_to_rownames("input_barcode")
  
  meta.data <- meta.data[row.names(s.obj@meta.data),]
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$group_mutation, col.name = "group_mutation")
  s.obj <- AddMetaData(object = s.obj, metadata = meta.data$num_mutation, col.name = "num_mutation")
  
  if (file.exists(file.path(path.to.06.output, "DGE_high_low_mutation_rate.rds")) == FALSE){
    dge.mutation <- list()
    dge.mutation.raw <- list()
    
    dge.mutation.raw[["high_vs_low"]] <- FindMarkers(object = s.obj, 
                                                     assay = "RNA",
                                                     ident.1 = "high", 
                                                     ident.2 = "low",
                                                     group.by = "group_mutation") 
    dge.mutation[["high_vs_low"]] <- dge.mutation.raw[["high_vs_low"]] %>%
      subset(p_val_adj <= 0.05) %>%
      rownames_to_column("Gene") %>%
      rowwise() %>%
      mutate(abs_avg_log2FC = abs(avg_log2FC))
    
    dge.mutation.raw[["high_vs_intermediate"]] <- FindMarkers(object = s.obj, 
                                                              assay = "RNA",
                                                              ident.1 = "high", 
                                                              ident.2 = "intermediate",
                                                              group.by = "group_mutation") 
    dge.mutation[["high_vs_intermediate"]] <- dge.mutation.raw[["high_vs_intermediate"]] %>%
      subset(p_val_adj <= 0.05) %>%
      rownames_to_column("Gene") %>%
      rowwise() %>%
      mutate(abs_avg_log2FC = abs(avg_log2FC))
    
    dge.mutation.raw[["intermediate_vs_low"]] <- FindMarkers(object = s.obj, 
                                                             assay = "RNA",
                                                             ident.1 = "intermediate", 
                                                             ident.2 = "low",
                                                             group.by = "group_mutation") 
    dge.mutation[["intermediate_vs_low"]] <- dge.mutation.raw[["intermediate_vs_low"]] %>%
      subset(p_val_adj <= 0.05) %>%
      rownames_to_column("Gene") %>%
      rowwise() %>%
      mutate(abs_avg_log2FC = abs(avg_log2FC))
    
    saveRDS(dge.mutation.raw, file.path(path.to.06.output, "DGE_high_low_mutation_rate.raw.rds"))
    saveRDS(dge.mutation, file.path(path.to.06.output, "DGE_high_low_mutation_rate.rds"))
    
  } else {
    dge.mutation <- readRDS(file.path(path.to.06.output, "DGE_high_low_mutation_rate.rds"))
  }
  
  for (g in names(dge.mutation.raw)){
    input.df <- dge.mutation.raw[[g]] %>%
      rownames_to_column("Gene") %>%
      rowwise() %>%
      mutate(sig = ifelse(p_val_adj <= 0.05, "sig", "non-sig")) %>%
      mutate(show.gene.name = ifelse(sig == "sig", Gene, NA)) %>%
      mutate(abs.avg_log2FC = abs(avg_log2FC))
    
    cutoff.adjp <- 0.05
    
    write.csv(input.df, file.path(path.to.06.output, sprintf("inputVolcanoPlot_%s.csv", g) ))
    volcano.plot <- ggplot(data=input.df, 
                           aes(x=avg_log2FC, y=-log10(p_val_adj), col=sig, label=show.gene.name)) + 
      geom_point() + 
      scale_color_manual(values=c("#c0d2f0", "#f28095")) +
      theme_minimal() +
      geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
      geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
      geom_text_repel() +
      theme_bw() + 
      theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
      xlim(c(-max(input.df$abs.avg_log2FC), max(input.df$abs.avg_log2FC)))
    
    ggsave(plot = volcano.plot, 
           filename = sprintf("VolcanoPlot_%s.csv", g), 
           path = path.to.06.output, 
           dpi = 300, 
           width = 14, 
           height = 10)
  }
  
  for (n in names(dge.mutation)){
    writexl::write_xlsx(dge.mutation[[n]] %>% 
                          arrange(desc(abs_avg_log2FC)) %>% 
                          arrange(desc(abs_avg_log2FC)), 
                        file.path(path.to.06.output, sprintf("DGE_test_result_%s.xlsx", n)))
  }
}


