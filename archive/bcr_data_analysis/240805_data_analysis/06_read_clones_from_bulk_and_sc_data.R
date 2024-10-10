#> 06: read clones from bulk BCR data and the single cell data. 
#> The single cell dataset consists of 6 samples: M1 M2 M3 (MLN samples) and 
#> P1 P2 P3 (PP samples). 
#> See the file "240829 sample sheet.xlsx" for a metadata of the bulk dataset. 
#> We import all clones and place them into one single dataframe (13 bulk samples
#> and 6 single cell samples). For further analysis, one can easily extract 
#> the desired samples from this dataframe. 
#> 
gc()
rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.project.src <- "/media/hieunguyen/HNSD01/src/bcr_data_analysis/240805_data_analysis"
source(file.path(path.to.project.src, "config.R"))

outdir <- "/media/hieunguyen/HNSD01/outdir"
PROJECT <- "240805_BSimons"
output.version <- "20240820"
config.version <- "default"
chosen.quantile <- 0.85
integration.case <- "all_samples"
threshold <- 0.85

path.to.main.input <- file.path(outdir, PROJECT, output.version, config.version)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis", output.version, config.version)

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.04.output <- file.path(path.to.main.output, "04_output", sprintf("quantile_%s", chosen.quantile))
path.to.05.output <- file.path(path.to.main.output, "05_output", sprintf("quantile_%s", chosen.quantile), integration.case)
path.to.06.output <- file.path(path.to.main.output, "06_output", sprintf("quantile_%s", chosen.quantile), integration.case)
dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)

if (file.exists(file.path(path.to.06.output, "combine_clonedf.xlsx")) == FALSE){
  path.to.bulk.data <- file.path(outdir, "240826_BSimons", "data_analysis")
  path.to.04.bulk.output <- file.path(path.to.bulk.data, "04_output")
  path.to.02.bulk.output <- file.path(path.to.bulk.data, "02_output")
  
  if (file.exists(file.path(path.to.06.output, "sc_clones.xlsx")) == FALSE){
    all.clone.files <- Sys.glob(file.path(path.to.04.bulk.output, "CDR3_0.15", "clonesets_*.xlsx"))
    all.clone.files <- all.clone.files[grepl("split", all.clone.files) == FALSE]
    
    names(all.clone.files) <- to_vec(for (item in basename(all.clone.files)) str_replace(str_replace(item, "clonesets_", ""), ".xlsx", ""))
    
    bulkdf <- data.frame()
    for (i in names(all.clone.files)){
      file <- all.clone.files[[i]]
      tmpdf <- readxl::read_excel(file)
      tmpdf$MID <- i 
      bulkdf <- rbind(bulkdf, tmpdf)
    }
    bulk.selected.cols <- c("V.gene", "J.gene", "nSeqCDR3", "aaSeqCDR3", "MID", "name")
    bulkdf.raw <- bulkdf
    bulkdf$V.gene <- unlist(lapply(bulkdf$V.gene, function(x){
      str_split(x, "[*]")[[1]][[1]]
    }))
    bulkdf$J.gene <- unlist(lapply(bulkdf$J.gene, function(x){
      str_split(x, "[*]")[[1]][[1]]
    }))
    bulkdf$name <- bulkdf$MID
    bulkdf <- bulkdf[, bulk.selected.cols]
    path.to.sc.data <- file.path(outdir, "240805_BSimons", "20240820", "default", "VDJ_output")
    all.vdj.outputs <- Sys.glob(file.path(path.to.sc.data, "annotated_contigs_clonaltype*.csv"))
    
    scdf <- data.frame()
    for (file in all.vdj.outputs){
      sample.id <- str_replace(str_replace(basename(file), "annotated_contigs_clonaltype_", ""), ".csv", "")
      tmpdf <- read.csv(file)
      tmpdf$name <- sample.id
      scdf <- rbind(scdf, tmpdf)
    }
    scdf.raw <- scdf
    
    scdf <- subset(scdf, select = -c(X))
    selected.cols <- c("barcode", "name", "IGH", "cdr3_aa1", "cdr3_nt1")
    scdf <- scdf[, selected.cols]
    
    # remove NA clones
    scdf <- subset(scdf, is.na(scdf$IGH) == FALSE)
    
    ##### convert scdf to match bulkdf
    scdf$MID <- scdf$name
    scdf$nSeqCDR3 <- scdf$cdr3_nt1
    scdf$aaSeqCDR3 <- scdf$cdr3_aa1
    
    scdf$V.gene <- unlist(lapply(scdf$IGH, function(x){
      str_split(x, "[.]")[[1]][[1]]
    }))
    scdf$J.gene <- unlist(lapply(scdf$IGH, function(x){
      str_split(x, "[.]")[[1]][[3]]
    }))
    
    scdf <- scdf[, c("barcode", bulk.selected.cols)]
    
    bulkdf$barcode <- to_vec( for(item in seq(1, nrow(bulkdf))) sprintf("seq%s", item) )
    writexl::write_xlsx(scdf, file.path(path.to.06.output, "sc_clones.xlsx"))
    writexl::write_xlsx(bulkdf, file.path(path.to.06.output, "bulk_clones.xlsx"))
  } else {
    scdf <- readxl::read_excel(file.path(path.to.06.output, "sc_clones.xlsx"))
    bulkdf <- readxl::read_excel(file.path(path.to.06.output, "bulk_clones.xlsx"))
  }
  
  scdf$origin <- "single-cell"
  bulkdf$origin <- "bulk"
  
  clonedf <- rbind(scdf, bulkdf[, colnames(scdf)])
  
  ##### add some more information on the sequences
  clonedf <- clonedf %>% rowwise() %>% 
    mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)) )
  
  sample.metadata <- readxl::read_excel(file.path(path.to.project.src, "240829 sample sheet.xlsx"))
  convert.label <- list()
  for (mid in unique(sample.metadata$MID)){
    organ <- subset(sample.metadata, sample.metadata$MID == mid)$sample
    organ <- paste(str_split(organ, " ")[[1]][2:3], collapse = "_")
    organ <- str_replace(organ, "_NA", "")
    convert.label[[mid]] <- organ
  }
  for (mid in c("M1", "M2", "M3", "P1", "P2", "P3")){
    if (grepl("M", mid) == TRUE){
      convert.label[[mid]] <- "MLN"
    } else {
      convert.label[[mid]] <- "PP"
    }
  }
  
  clonedf$organ <- unlist(lapply(
    clonedf$name, function(x){
      return(convert.label[[x]])
    }
  ))
  clonedf$mouseid <- unlist(lapply(
    clonedf$MID, function(x){
      if (x %in% c(unique(scdf$name))){
        return(sprintf("m%s", str_replace(str_replace(x, "M", ""), "P", "")))
      } else {
        return(subset(sample.metadata, sample.metadata$MID == x)$mouse)
      }
    }
  ))
  writexl::write_xlsx(clonedf, file.path(path.to.06.output, "combine_clonedf.xlsx"))
} else {
  clonedf <- readxl::read_excel(file.path(path.to.06.output, "combine_clonedf.xlsx"))
  scdf <- readxl::read_excel(file.path(path.to.06.output, "sc_clones.xlsx"))
  bulkdf <- readxl::read_excel(file.path(path.to.06.output, "bulk_clones.xlsx"))
}

clonedf <- clonedf %>% rowwise() %>% 
  mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)) )

# Group clone by sequence similarity and V-J genes
source(file.path(path.to.project.src, "helper_functions.R"))
clonedf <- clonedf %>% rowwise() %>% 
  mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)) )

##### Assign clones to clusters
if (file.exists(file.path(path.to.06.output, "full_clonedf.xlsx")) == FALSE){
  new.clonedf <- data.frame()
  for (input.VJ.combi in unique(clonedf$VJ.len.combi)){
    tmpdf <- subset(clonedf, clonedf$VJ.len.combi == input.VJ.combi)
    seqs <- unique(tmpdf$aaSeqCDR3)
    if (length(seqs) >= 2){
      cluster.output <- assign_clusters_to_sequences(seqs, threshold = threshold)$res
      tmpdf[[sprintf("VJcombi_CDR3_%s", threshold)]] <- unlist(lapply(
        tmpdf$aaSeqCDR3, function(x){
          return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
        }
      ))    
    } else {
      tmpdf[[sprintf("VJcombi_CDR3_%s", threshold)]] <- sprintf("%s_1", input.VJ.combi)
    }
    new.clonedf <- rbind(new.clonedf, tmpdf)
  }
  writexl::write_xlsx(new.clonedf, file.path(path.to.06.output, "full_clonedf.xlsx"))
  
  cloneCluster.summary <- data.frame(clone = unique(new.clonedf[[sprintf("VJcombi_CDR3_%s", threshold)]])) %>% 
    rowwise() %>%
    mutate(num_samples = length(unique(subset(new.clonedf, new.clonedf[[sprintf("VJcombi_CDR3_%s", threshold)]] == clone)$name))) %>%
    mutate(num_total_seq = nrow(subset(new.clonedf, new.clonedf[[sprintf("VJcombi_CDR3_%s", threshold)]] == clone))) %>%
    mutate(num_unique_seq = length(unique(subset(new.clonedf, new.clonedf[[sprintf("VJcombi_CDR3_%s", threshold)]] == clone)$nSeqCDR3))) %>%
    arrange(desc(num_total_seq))
  writexl::write_xlsx(cloneCluster.summary, file.path(path.to.06.output, "cloneCluster_summary.xlsx"))
  
} else {
  new.clonedf <- readxl::read_excel(file.path(path.to.06.output, "full_clonedf.xlsx"))
}


