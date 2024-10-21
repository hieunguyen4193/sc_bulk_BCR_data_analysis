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

path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240805_BSimons"
sample.metadata <- readxl::read_excel(file.path(path.to.project.src, "240829 sample sheet.xlsx"))

outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "240805_BSimons"
PROJECT.bulk <- "240826_BSimons"

chosen.quantile <- 0.85
thres <- 0.85

all.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
all.quantiles <- c(0.99, 0.95, 0.90, 0.85)
all.integration.case <- list(
  all_samples = all.samples,
  mouse1 = c("M1", "P1"),
  mouse2 = c("M2", "P2"),
  mouse3 = c("M3", "P3")
)

all.bulk.mouses <- list(
  all_samples = c("m1", "m2", "m3"),
  mouse1 = c("m1"),
  mouse2 = c("m2"),
  mouse3 = c("m3")
)

#####----------------------------------------------------------------------#####
##### GENERATE LABEL (ORGANS) FOR ALL BULK AND SC SAMPLES
#####----------------------------------------------------------------------#####
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

#####----------------------------------------------------------------------#####
##### MAIN RUN
#####----------------------------------------------------------------------#####
for (integration.case in names(all.integration.case)){
  path.to.main.input <- file.path(outdir, PROJECT)
  path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
  
  path.to.01.output <- file.path(path.to.main.output, "01_output")
  path.to.04.output <- file.path(path.to.main.output, "04_output", sprintf("quantile_%s", chosen.quantile))
  path.to.05.output <- file.path(path.to.main.output, "05_output", sprintf("quantile_%s", chosen.quantile), integration.case)
  path.to.06.output <- file.path(path.to.main.output, "06_output", sprintf("quantile_%s", chosen.quantile), integration.case)
  dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.06.output, "combine_clonedf.xlsx")) == FALSE){
    if (file.exists(file.path(path.to.06.output, "sc_clones.xlsx")) == FALSE){
      #####----------------------------------------------------------------#####
      ##### PROCESS bulk clone data frames
      #####----------------------------------------------------------------#####
      all.clone.files <- Sys.glob(file.path(outdir, 
                                            "VDJ_output", 
                                            PROJECT.bulk, 
                                            sprintf("VDJ_output_%s", thres), 
                                            "preprocessed_files", 
                                            "clonesets_*.xlsx"))
      names(all.clone.files) <- to_vec(
          for (item in basename(all.clone.files)) 
          str_replace(str_replace(item, "clonesets_", ""), ".xlsx", "")
        )
      selected.MIDs <- subset(sample.metadata, 
                              sample.metadata$mouse %in% all.bulk.mouses[[integration.case]])$MID
      
      bulkdf <- data.frame()
      for (i in names(all.clone.files)){
        file <- all.clone.files[[i]]
        tmpdf <- readxl::read_excel(file)
        bulkdf <- rbind(bulkdf, tmpdf)
      }
      bulkdf$MID <- bulkdf$id
      
      ##### keep sequences/clones in selected MIDs only
      bulkdf <- subset(bulkdf, bulkdf$MID %in% selected.MIDs)
      print(sprintf("Selected MIDs (bulk) in this analysis: %s", paste(selected.MIDs, collapse = ", ")))
      print(sprintf("Selected mice (bulk) in this analysis: %s", paste(all.bulk.mouses[[integration.case]], collapse = ", ")))
      bulk.selected.cols <- c("V.gene", "J.gene", "nSeqCDR3", "aaSeqCDR3", "MID")
      bulkdf.raw <- bulkdf
      
      ##### remove * in V J gene name, to match single cell gene annotation
      bulkdf$V.gene <- unlist(lapply(bulkdf$V.gene, function(x){
        str_split(x, "[*]")[[1]][[1]]
      }))
      bulkdf$J.gene <- unlist(lapply(bulkdf$J.gene, function(x){
        str_split(x, "[*]")[[1]][[1]]
      }))
      bulkdf <- bulkdf[, bulk.selected.cols]
      # generate pseudo barcode for the sequence in bulkdf
      bulkdf$barcode <- to_vec( for(item in seq(1, nrow(bulkdf))) sprintf("seq%s", item) )
      
      #####----------------------------------------------------------------#####
      ##### PROCESS single cell clones data frames
      #####----------------------------------------------------------------#####
      scdf <- readxl::read_excel(file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files", sprintf("clonesets_%s.split_clones.xlsx", PROJECT)))
      scdf$MID <- scdf$id
      
      scdf <- subset(scdf, scdf$MID %in% all.integration.case[[integration.case]])
      print(sprintf("Selected MIDs (single cell data) in this analysis: %s", paste(unique(scdf$MID), collapse = ", ")))
      
      # remove NA clones
      ##### convert scdf to match bulkdf
      
      scdf <- scdf[, c("barcode", bulk.selected.cols)]
      writexl::write_xlsx(scdf, file.path(path.to.06.output, "sc_clones.xlsx"))
      writexl::write_xlsx(bulkdf, file.path(path.to.06.output, "bulk_clones.xlsx"))
    } else {
      scdf <- readxl::read_excel(file.path(path.to.06.output, "sc_clones.xlsx"))
      bulkdf <- readxl::read_excel(file.path(path.to.06.output, "bulk_clones.xlsx"))
    }
    
    #####------------------------------------------------------------------#####
    ##### merge dataframes: scdf and bulkdf
    #####------------------------------------------------------------------#####
    scdf$origin <- "single-cell"
    bulkdf$origin <- "bulk"
    clonedf <- rbind(scdf, bulkdf[, colnames(scdf)])
    
    ##### add some more information on the sequences
    clonedf <- clonedf %>% rowwise() %>% 
      mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)) )
    
    ##### add organ information on each sample
    clonedf$organ <- unlist(lapply(
      clonedf$MID, function(x){
        return(convert.label[[x]])
      }
    ))
    
    ##### unify mouse ID
    clonedf$mouseid <- unlist(lapply(
      clonedf$MID, function(x){
        if (x %in% c(unique(scdf$MID))){
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
  #####--------------------------------------------------------------------#####
  ##### ASSIGN CLONES FROM BULK AND SC DATA TO GROUPS BASED ON THEIR SIMILARITY
  ##### OF THE CDR3 SEQUENCE.
  #####--------------------------------------------------------------------#####
  clonedf <- clonedf %>% rowwise() %>% 
    mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)) )
  
  ##### Group clone based on V-J genes and CDR3 sequence length
  source(file.path(path.to.project.src, "helper_functions.R"))
  clonedf <- clonedf %>% rowwise() %>% 
    mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)) )
  
  #####--------------------------------------------------------------------#####
  ##### Assign clones to clusters
  #####--------------------------------------------------------------------#####
  if (file.exists(file.path(path.to.06.output, "full_clonedf.xlsx")) == FALSE){
    new.clonedf <- data.frame()
    for (input.VJ.combi in unique(clonedf$VJ.len.combi)){
      tmpdf <- subset(clonedf, clonedf$VJ.len.combi == input.VJ.combi)
      seqs <- unique(tmpdf$aaSeqCDR3)
      if (length(seqs) >= 2){
        cluster.output <- assign_clusters_to_sequences(seqs, thres = thres)$res
        tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
          tmpdf$aaSeqCDR3, function(x){
            return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
          }
        ))    
      } else {
        tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
      }
      new.clonedf <- rbind(new.clonedf, tmpdf)
    }
    writexl::write_xlsx(new.clonedf, file.path(path.to.06.output, "full_clonedf.xlsx"))
    
    #####--------------------------------------------------------------------#####
    #####> This is the final dataframe containing all clones from the selected
    #####> bulk and single cell mouse / mice. 
    #####> Clones are grouped based on V, J genes and the thres% similarity of 
    #####> the CDR3 sequences. 
    #####--------------------------------------------------------------------#####
    cloneCluster.summary <- data.frame(clone = unique(new.clonedf[[sprintf("VJcombi_CDR3_%s", thres)]])) %>% 
      rowwise() %>%
      mutate(num_samples = length(unique(subset(new.clonedf, new.clonedf[[sprintf("VJcombi_CDR3_%s", thres)]] == clone)$name))) %>%
      mutate(num_total_seq = nrow(subset(new.clonedf, new.clonedf[[sprintf("VJcombi_CDR3_%s", thres)]] == clone))) %>%
      mutate(num_unique_seq = length(unique(subset(new.clonedf, new.clonedf[[sprintf("VJcombi_CDR3_%s", thres)]] == clone)$nSeqCDR3))) %>%
      arrange(desc(num_total_seq))
    writexl::write_xlsx(cloneCluster.summary, file.path(path.to.06.output, "cloneCluster_summary.xlsx"))
  } else {
    new.clonedf <- readxl::read_excel(file.path(path.to.06.output, "full_clonedf.xlsx"))
  }
}



