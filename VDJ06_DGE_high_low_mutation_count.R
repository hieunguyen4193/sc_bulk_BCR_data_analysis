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
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

input.dataset <- "1st_2nd_BSimons_Datasets"
input.dataset.1or2 <- "1st_dataset"
path.to.04.output <- file.path(outdir, "VDJ_output", "04_output")

if (input.dataset == "1st_2nd_BSimons_Datasets"){
  sample.lists <- list(
    `1st_dataset` = c(
      "17_MM9_Ecoli",
      "20_MM9_Ecoli",
      "21_MM9_Ecoli_SPF",
      "MM9_S2",
      "MM9_S4",
      "MM9_SPF_S3",
      "MM9_SPF_S9"),
    `2nd_dataset` = c("Sample_132",
                      "Sample_133" )
  )
  s.obj <- readRDS(path.to.all.s.obj[[input.dataset.1or2]])
  path.to.06.output <- file.path(outdir, "VDJ_output", "06_output", input.dataset.1or2)
  dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
  clonedf <- read.csv(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"), row.names = "X") %>%
    subset(dataset.name == input.dataset) %>%
    rowwise() %>%
    mutate(group_mutation = asssign_mutation_to_group(num_mutation)) %>%
    mutate(barcode_full = sprintf("%s_%s_%s", id, id, barcode)) %>%
    subset(id %in% sample.lists[[input.dataset.1or2]]) 
  
} else{
  path.to.06.output <- file.path(outdir, "VDJ_output", "06_output", input.dataset)
  dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
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
all.clone.files <- Sys.glob(file.path(outdir, "VDJ_output", "*", sprintf("VDJ_output_%s", thres), "preprocessed_files", "clonesets*.split_clones.xlsx" ))

dataset.origin <- list(
  `1st_2nd_BSimons_Datasets` = "sc",
  `220701_etc_biopsies` = "bulk",
  `240805_BSimons` = "sc",
  `240826_BSimons` = "bulk",
  `241002_BSimons` = "sc",
  `241031_BSimons` = "bulk",
  `241104_BSimons` = "sc"
)

names(all.clone.files) <- to_vec(
  for (item in all.clone.files) str_replace(str_replace(basename(item), "clonesets_", ""), ".split_clones.xlsx", "")
)

meta.data <- s.obj@meta.data %>% 
  rownames_to_column("input_barcode")
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
  dge.mutation[["high_vs_low"]] <- FindMarkers(object = s.obj, 
                                               assay = "RNA",
                                               ident.1 = "high", 
                                               ident.2 = "low",
                                               group.by = "group_mutation") %>% 
    subset(p_val_adj <= 0.05) %>%
    rownames_to_column("Gene") %>%
    rowwise() %>%
    mutate(abs_avg_log2FC = abs(avg_log2FC))
  dge.mutation[["intermediate_vs_low"]] <- FindMarkers(object = s.obj, 
                                                       assay = "RNA",
                                                       ident.1 = "intermediate", 
                                                       ident.2 = "low",
                                                       group.by = "group_mutation") %>% 
    subset(p_val_adj <= 0.05) %>%
    rownames_to_column("Gene") %>%
    rowwise() %>%
    mutate(abs_avg_log2FC = abs(avg_log2FC))
  dge.mutation[["high_vs_intermediate"]] <- FindMarkers(object = s.obj, 
                                                        assay = "RNA",
                                                        ident.1 = "high", 
                                                        ident.2 = "intermediate",
                                                        group.by = "group_mutation") %>% 
    subset(p_val_adj <= 0.05) %>%
    rownames_to_column("Gene") %>%
    rowwise() %>%
    mutate(abs_avg_log2FC = abs(avg_log2FC))
  
  saveRDS(dge.mutation, file.path(path.to.06.output, "DGE_high_low_mutation_rate.rds"))
} else {
  dge.mutation <- readRDS(file.path(path.to.06.output, "DGE_high_low_mutation_rate.rds"))
}

for (n in names(dge.mutation)){
  writexl::write_xlsx(dge.mutation[[n]] %>% 
                        arrange(desc(abs_avg_log2FC)) %>% 
                        arrange(desc(abs_avg_log2FC)), 
                      file.path(path.to.06.output, sprintf("DGE_test_result_%s.xlsx", n)))
}

