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

countdf <- data.frame()
for (input.dataset in names(path.to.all.s.obj)){
  path.to.03.output <- file.path(outdir, "VDJ_output", "03_output")
  path.to.12.output <- path.to.05.output <- file.path(outdir, "VDJ_output", "12_output", input.dataset)
  dir.create(path.to.12.output, showWarnings = FALSE, recursive = TRUE)
  
  s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
  meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
  
  dataset.name <- input.dataset
  clone.case <- "single_cell"

  if (input.dataset %in% c("240805_BSimons",
                           "241104_BSimons",
                           "240805_BSimons_filterHT",
                           "240805_BSimons_filterHT_cluster_renamed",
                           "241002_BSimons",
                           "240805_BSimons_filterHT_cluster",
                           "241002_241104_BSimons")){
    for (mouse.id in unique(meta.data$name)){
      tmp.metadata <- subset(meta.data, meta.data$name == mouse.id)
      for (input.mid in unique(tmp.metadata$HTO_classification)){
        all.clones <- unique(subset(meta.data, meta.data$name == mouse.id & 
                                      meta.data$HTO_classification == input.mid)$VJcombi_CDR3_0.85)
        all.clones <- all.clones[is.na(all.clones) == FALSE]
        clone.count <- length(all.clones)
        tmp.outputdf <- data.frame(dataset = c(dataset.name),
                                   clone.case = c(clone.case),
                                   mouse = c(mouse.id),
                                   hashtag = c(input.mid),
                                   clone.count = c(clone.count)) 
        countdf <- rbind(countdf, tmp.outputdf)
      }
    }
  } else {
    for (mouse.id in unique(meta.data$name)){
      tmp.metadata <- subset(meta.data, meta.data$name == mouse.id)
        all.clones <- unique(subset(meta.data, meta.data$name == mouse.id)$VJcombi_CDR3_0.85)
        all.clones <- all.clones[is.na(all.clones) == FALSE]
        clone.count <- length(all.clones)
        tmp.outputdf <- data.frame(dataset = c(dataset.name),
                                   clone.case = c(clone.case),
                                   mouse = c(mouse.id),
                                   hashtag = c("no # information"),
                                   clone.count = c(clone.count)) 
        countdf <- rbind(countdf, tmp.outputdf)
    }
  }
}

writexl::write_xlsx(countdf, file.path(path.to.12.output, "count_clones_in_samples.single_cell.xlsx"))

all.bulk.cloneinfo <- Sys.glob(file.path(path.to.03.output, "/FASTA_output/*/VDJ_output_0.85/*/*/clonesets.csv"))

countdf <- data.frame()
for (input.file in all.bulk.cloneinfo){
  dataset.name <- str_split(input.file, "/")[[1]][[11]]
  clone.case <- str_split(input.file, "/")[[1]][[14]]
  mouse.id <- str_split(input.file, "/")[[1]][[13]]
  
  tmpdf <- read.csv(input.file)
  for (input.mid in unique(tmpdf$id)){
    clone.count <- length(unique(subset(tmpdf, tmpdf$id == input.mid)$VJcombi_CDR3_0.85))
    tmp.outputdf <- data.frame(dataset = c(dataset.name),
                               clone.case = c(clone.case),
                               mouse = c(mouse.id),
                               SampleID = c(input.mid),
                               clone.count = c(clone.count))  
    countdf <- rbind(countdf, tmp.outputdf)
  }
}
writexl::write_xlsx(countdf, file.path(path.to.12.output, "count_clones_in_samples.bulk.xlsx"))
