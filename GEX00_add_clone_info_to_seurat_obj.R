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
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

# dataset.name <- "241002_BSimons"
for (dataset.name in unique(names(path.to.all.VDJ.output)) ){
  path.to.00.output <- file.path(outdir, "GEX_output", "00_output", dataset.name)
  dir.create(path.to.00.output, showWarnings = FALSE, recursive = TRUE)
  savename <- str_replace(basename(path.to.all.s.obj[[dataset.name]]), ".rds", "")
  if (file.exists(file.path(path.to.00.output, sprintf("%s.addedCloneInfo.rds", savename))) == FALSE){
    s.obj <- readRDS(path.to.all.s.obj[[dataset.name]])
    DefaultAssay(s.obj) <- "RNA"
    
    vdj.output <- readRDS(path.to.all.VDJ.output[[dataset.name]])
    vdj.output <- vdj.output[unique(s.obj$name)]
    
    vdjdf <- data.frame()
    for (sample.id in names(vdj.output)){
      vdjdf <- rbind(vdjdf, vdj.output[[sample.id]])
    }
    
    meta.data <- s.obj@meta.data %>% rownames_to_column("barcode")
    ct.columns <- to_vec(
      for (item in colnames(vdjdf)){
        if (grepl("CT", item) == TRUE | grepl("cdr3", item) == TRUE){
          item
        }
      }
    )
    
    ct.columns <- c(ct.columns, "IGH", "IGLC")
    meta.data <- meta.data[setdiff(colnames(meta.data), ct.columns)]
    meta.data <- merge(meta.data, vdjdf, by.x = "barcode", by.y = "barcode", all.x = TRUE) %>%
      rowwise() %>%
      mutate(V.gene = ifelse(is.na(IGH) == FALSE, str_split(IGH, "[.]")[[1]][[1]], NA)) %>%
      mutate(J.gene = ifelse(is.na(IGH) == FALSE, str_split(IGH, "[.]")[[1]][[3]], NA)) %>%
      mutate(lenCDR3 = ifelse(is.na(IGH) == FALSE, nchar(cdr3_nt1), NA)) %>%
      mutate(aaSeqCDR3 = ifelse(is.na(IGH) == FALSE, cdr3_aa1, NA) ) %>%
      mutate(nSeqCDR3 = ifelse(is.na(IGH) == FALSE, cdr3_nt1, NA) ) %>%
      mutate(VJseq.combi = sprintf("%s_%s_%s_%s", V.gene, J.gene, aaSeqCDR3, nSeqCDR3)) %>%
      mutate(VJ.combi = sprintf("%s_%s", V.gene, J.gene)) %>%
      mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3))) 
    
    clonesets <- meta.data[, c("barcode", "V.gene", "J.gene", "aaSeqCDR3", "nSeqCDR3", "VJseq.combi", "VJ.combi", "VJ.len.combi")] %>%
      subset(VJ.combi != "NA_NA")
    thres.dis <- 0.15
    thres <- 0.85
    
    clonedf <- data.frame()
    for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
      tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
      seqs <- unique(tmpdf$aaSeqCDR3)
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
      clonedf <- rbind(clonedf, tmpdf)
    }
    
    meta.data[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
      meta.data$barcode, function(x){
        if (x %in% clonedf$barcode){
          return(subset(clonedf, clonedf$barcode == x)[[sprintf("VJcombi_CDR3_%s", thres)]])
        } else {
          return(NA)
        }
      }
    ))
    
    meta.data <- meta.data %>% column_to_rownames("barcode") 
    meta.data <- meta.data[row.names(s.obj@meta.data), ]
    new.cols <- c("CTgene",
                  "CTnt",
                  "CTaa",
                  "CTstrict",
                  "V.gene", 
                  "J.gene", 
                  "aaSeqCDR3", 
                  "nSeqCDR3", 
                  "VJseq.combi", 
                  "VJ.combi", 
                  "VJ.len.combi")
    for (c in new.cols){
      s.obj <- AddMetaData(object = s.obj, col.name = c, metadata = meta.data[[c]])
    }
    
    saveRDS(s.obj, file.path(path.to.00.output, sprintf("%s.addedCloneInfo.rds", savename)) )
  } else {
    print(sprintf("File exists at %s", file.path(path.to.00.output, sprintf("%s.addedCloneInfo.rds", savename))))
  }
}

