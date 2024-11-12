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
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.04.output <- file.path(outdir, "VDJ_output", "04_output")
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

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

selected.cols <- c(
  "id",
  "VJseq.combi",
  "V.gene",
  "J.gene",
  "D.gene",
  "nSeqFR1",
  "nSeqCDR1",
  "nSeqFR2",
  "nSeqCDR2",
  "nSeqFR3",
  "nSeqCDR3",
  "nSeqFR4",
  "aaSeqFR1",
  "aaSeqCDR1",
  "aaSeqFR2",
  "aaSeqCDR2",
  "aaSeqFR3",
  "aaSeqCDR3",
  "aaSeqFR4",
  "VJ.len.combi",
  "targetSequences",
  "uniqueMoleculeCount")

rerun <- TRUE

if (file.exists(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv")) == FALSE | rerun == TRUE){
  clonedf <- data.frame()
  for (dataset.name in names(all.clone.files)){
    input.file <- all.clone.files[[dataset.name]]
    tmpdf <- readxl::read_excel(input.file)
    tmpdf$dataset <- dataset.name
    dataset.type <- dataset.origin[[dataset.name]]
    if (dataset.type == "sc"){
      tmpdf <- tmpdf[, c(setdiff(selected.cols, c("targetSequences", "uniqueMoleculeCount")), "barcode")]
      tmpdf <- tmpdf %>% rowwise() %>%
        mutate(targetSequences = paste0(c(
          nSeqFR1,
          nSeqCDR1,
          nSeqFR2,
          nSeqCDR2,
          nSeqFR3,
          nSeqCDR3,
          nSeqFR4
        ), collapse = ""))
      tmpdf$uniqueMoleculeCount <- NA
    } else {
      tmpdf <- tmpdf[, selected.cols]
      tmpdf$barcode <- to_vec( for (i in seq(1, nrow(tmpdf))) sprintf("%s_%s", dataset.name, i))
    }
    tmpdf$sampletype <- dataset.type
    tmpdf$len_aaSeqCDR3 <- unlist(lapply(tmpdf$aaSeqCDR3, function(x){
      nchar(x)
    }))
    tmpdf$len_nSeqCDR3 <- unlist(lapply(tmpdf$nSeqCDR3, function(x){
      nchar(x)
    }))
    tmpdf$dataset.name <- dataset.name
    tmpdf$V.gene <- unlist(
      lapply(tmpdf$V.gene, function(x){
        if (grepl("*", x) == TRUE){
          x <- str_split(x, "[*]")[[1]][[1]]    
        }
        return(x)
      })
    )
    tmpdf$J.gene <- unlist(
      lapply(tmpdf$J.gene, function(x){
        if (grepl("*", x) == TRUE){
          x <- str_split(x, "[*]")[[1]][[1]]    
        }
        return(x)
      })
    )
    clonedf <- rbind(clonedf, tmpdf)
  }
  
  #####---------------------------------------------------------------------------#####
  ##### GET GERMLINE SEQUENCES
  #####---------------------------------------------------------------------------#####
  ref.gene.config <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/ref_gene_config.R"
  source(ref.gene.config) # path to the configuration file for the input BCR reference genes
  
  cellranger.ref <- readDNAStringSet(ref.genes$`10x`)
  names(cellranger.ref) <- to_vec(
    for (item in names(cellranger.ref)) str_split(str_split(item, " ")[[1]][[1]], "[|]")[[1]][[2]]
  )
  s.V.genes <- readDNAStringSet(ref.genes$IMGT$V.gene)
  names(s.V.genes) <- lapply(names(s.V.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    if (grepl("*", x) == TRUE){
      x <- str_split(x, "[*]")[[1]][[1]]    
    }
    return(x)
  })
  s.J.genes <- readDNAStringSet(ref.genes$IMGT$J.gene)
  names(s.J.genes) <- lapply(names(s.J.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    if (grepl("*", x) == TRUE){
      x <- str_split(x, "[*]")[[1]][[1]]    
    }
  })
  
  ##### check mutations
  count_mutation_for_line_i <- function(i){
    V.gene <- clonedf[i, ]$V.gene
    J.gene <- clonedf[i, ]$J.gene
    if (is.na(V.gene) == TRUE | is.na(J.gene) == TRUE){
      return("NA-skip")
    }
    CDR3.seq <- clonedf[i, ]$nSeqCDR3
    CDR3.length <- nchar(CDR3.seq)
    clone.seq <- paste0(clonedf[i, ]$nSeqFR1,
                        clonedf[i, ]$nSeqCDR1,
                        clonedf[i, ]$nSeqFR2,
                        clonedf[i, ]$nSeqCDR2,
                        clonedf[i, ]$nSeqFR3, 
                        clonedf[i, ]$nSeqCDR3,
                        clonedf[i, ]$nSeqFR4, collapse = "")
    if (grepl("region_not_covered", clone.seq) == TRUE){
      return("region_not_covered-skip")
    }
    sampletype <- clonedf[i, ]$sampletype
    if (sampletype == "bulk"){
      GL.V.gene <- s.V.genes[[V.gene]] %>% as.character()
      GL.J.gene <- s.J.genes[[J.gene]] %>% as.character()  
    } else if (sampletype == "sc"){
      GL.V.gene <- cellranger.ref[[V.gene]] %>% as.character()
      GL.J.gene <- cellranger.ref[[J.gene]] %>% as.character()  
    }
    
    repN.seq <- paste(replicate(n = CDR3.length, expr = "N"), collapse = "")
    
    GL.seq.full <- sprintf("%s%s%s", GL.V.gene, repN.seq, GL.J.gene)
    input.seqs <- c(clone.seq, GL.seq.full) %>% DNAStringSet()
    msa.res <- msa(input.seqs, method = "Muscle")
    aligned.clone.seq <- toString(unmasked(msa.res)[[1]])
    aligned.GL.seq <- toString(unmasked(msa.res)[[2]])
    
    s1 <- aligned.clone.seq
    s2 <- aligned.GL.seq
    sdf <- data.frame(s1 = str_split(s1, "")[[1]],
                      s2 = str_split(s2, "")[[1]])
    
    check_mutation <- function(x1, x2){
      if (x1 == "-" | x2 == "-"){
        output <- "undefined"
      } else if (x1 == "N" | x2 == "N"){
        output <- "undefined"
      } else {
        if (x1 == x2){
          output <- "matched"
        } else {
          output <- "mutation"
        }
      }
      return(output)
    }
    
    sdf <- sdf %>% rowwise() %>% mutate(status = check_mutation(s1, s2))
    count.mutation <- table(sdf$status)
    if ("mutation" %in% names(count.mutation) == FALSE){
      count.mutation[["mutation"]] <- 0
    }
    return(count.mutation[["mutation"]])
  }
  
  clonedf$num_mutation <- unlist(lapply(
    seq(1, nrow(clonedf)), function(x){
      if (x%%1000 == 0){
        print(sprintf("Working at step %s", x))
      }
      # print(x)
      count_mutation_for_line_i(x)
    }
  ))
  write.csv(clonedf, file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"))
} else {
  print(sprintf("Clone dataframe with mutation rate exists at %s, reading in ... ", file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv")))
  clonedf <- read.csv(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"), row.names = "X")
}

