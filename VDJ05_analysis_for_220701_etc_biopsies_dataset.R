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
##### READ INPUT SUMMARY CLONE DF
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.04.output <- file.path(outdir, "VDJ_output", "04_output")
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", "220701_etc_biopsies")
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
thres <- 0.85

all.clone.files <- Sys.glob(file.path(outdir, "VDJ_output", "*", sprintf("VDJ_output_%s", thres), "preprocessed_files", "clonesets*.split_clones.xlsx" ))

clonedf <- read.csv(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"), row.names = "X")
clonedf <- subset(clonedf, clonedf$num_mutation != "region_not_covered-skip")
clonedf <- subset(clonedf, clonedf$dataset.name == "220701_etc_biopsies")

clonedf$num_mutation <- as.numeric(clonedf$num_mutation)
mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")
count.mid.in.mouse <- table(mid.metadata$population, mid.metadata$mouse) %>% colSums()

#####----------------------------------------------------------------------#####
##### PREPARE LIST OF MID W.R.T EACH CASE, ALL SAMPLES, BIOPSY ONLY OR YPF+-
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.05.output, "sample_list_based_on_YFP.rds")) == FALSE){
  yfp.mids <- list()
  for (mouse.id in names(count.mid.in.mouse[count.mid.in.mouse >= 4])){
    yfp.mids[[mouse.id]] <- list(all = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
                                 pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$X,
                                 neg = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)$X,
                                 biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$population == "biopsy")$X)   
  }
  saveRDS(yfp.mids, file.path(path.to.05.output, "sample_list_based_on_YFP.rds"))
} else {
  yfp.mids <- readRDS(file.path(path.to.05.output, "sample_list_based_on_YFP.rds"))
}

#####----------------------------------------------------------------------#####
##### REF GENES
#####----------------------------------------------------------------------#####
ref.gene <- "IMGT"
ref.gene.config <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/ref_gene_config.R"
source(ref.gene.config) # path to the configuration file for the input BCR reference genes
if (ref.gene == "10x"){
  ref.fasta <- readDNAStringSet(ref.genes$`10x`)
} else if (ref.gene == "IMGT"){
  s.V.genes <- readDNAStringSet(ref.genes$IMGT$V.gene)
  names(s.V.genes) <- lapply(names(s.V.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    x <- str_split(x, "[*]")[[1]][[1]]
    return(x)
  })
  s.J.genes <- readDNAStringSet(ref.genes$IMGT$J.gene)
  names(s.J.genes) <- lapply(names(s.J.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    x <- str_split(x, "[*]")[[1]][[1]]
    return(x)
  })
}

#####----------------------------------------------------------------------#####
##### MAIN RUN: GENERATE FASTA FOR EACH YFP CASE AND BIOPSY SAMPLES
#####----------------------------------------------------------------------#####
summary.fastadf <- data.frame()
for (mouse.id in names(yfp.mids)){
  for (yfp.case in c("all", "neg", "pos", "biopsy")){
    
    if (yfp.case == "biopsy"){
      save.folder.name <- yfp.case
    } else {
      save.folder.name <- sprintf("%s_YFP", yfp.case)
    }
    
    path.to.save.fasta <- file.path(path.to.05.output, "FASTA", mouse.id, save.folder.name)
    dir.create(path.to.save.fasta, showWarnings = FALSE, recursive = TRUE)
    
    inputdf <- subset(clonedf, clonedf$id %in% yfp.mids[[mouse.id]][[yfp.case]])
    
    new.clonesetsdf <- data.frame()
    for (input.VJ.combi in unique(inputdf$VJ.len.combi)){
      tmpdf <- subset(inputdf, inputdf$VJ.len.combi == input.VJ.combi)
      seqs <- unique(tmpdf$aaSeqCDR3)
      if (length(seqs) >= 2){
        cluster.output <- assign_clusters_to_sequences(seqs, threshold = 1 - thres)$res
        tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
          tmpdf$aaSeqCDR3, function(x){
            return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
          }
        ))    
      } else {
        tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
      }
      new.clonesetsdf <- rbind(new.clonesetsdf, tmpdf)
    }
    
    for (input.VJ.combi in unique(new.clonesetsdf[[sprintf("VJcombi_CDR3_%s", thres)]])){
      path.to.output.fasta <- file.path(path.to.save.fasta, sprintf("%s_%s_%s.fasta", mouse.id, save.folder.name, input.VJ.combi))
      fasta.df <- subset(new.clonesetsdf, new.clonesetsdf[[sprintf("VJcombi_CDR3_%s", thres)]] == input.VJ.combi)
      
      V.gene <- str_split(input.VJ.combi, "_")[[1]][[1]]
      J.gene <- str_split(input.VJ.combi, "_")[[1]][[2]]
      CDR3.length <- as.numeric(str_split(input.VJ.combi, "_")[[1]][[3]])
      
      ##### get germline sequences and merge with the real data sequences
      GL.V.gene <- s.V.genes[[V.gene]] %>% as.character()
      GL.J.gene <- s.J.genes[[J.gene]] %>% as.character()
      repN.seq <- paste(replicate(n = CDR3.length, expr = "N"), collapse = "")
      GL.seq <- sprintf("%s%s%s", GL.V.gene, repN.seq, GL.J.gene)
      
      fasta.output <- fasta.df[, c("targetSequences",
                                   "uniqueMoleculeCount", 
                                   "V.gene", 
                                   "D.gene", 
                                   "J.gene", 
                                   "id", 
                                   "aaSeqCDR3", 
                                   "nSeqCDR3")]
      
      colnames(fasta.output) <- c("seq", "abundance", "V.gene", "D.gene", "J.gene", "SampleID", "CDR3aa", "CDR3nt")
      
      all.seqs <- c(fasta.output %>% pull(`seq`), GL.seq)
      
      #####----------------------------------------------------------------#####
      ##### start writing fasta file
      #####----------------------------------------------------------------#####
      if (nrow(fasta.output) > 1){
        ##### multiple alignment sequences, package MSA. 
        MiXCRtreeVDJ <- all.seqs %>% DNAStringSet()
        msaMiXCRtreeVDJ <- msa(inputSeqs = MiXCRtreeVDJ, verbose = TRUE)
        sink(path.to.output.fasta)
        for (i in seq(1, length(all.seqs))){
          if (i == length(all.seqs)){
            output.info <- ">GL"
          } else {
            sample.id <- fasta.output[i, ]$SampleID
            cdr3aa <- fasta.output[i, ]$CDR3aa
            cdr3nt <- fasta.output[i, ]$CDR3nt
            abundance <- fasta.output[i, ]$abundance
            output.info <- sprintf(">Sample:%s|Mouse:%s|CDR3aa:%s|CDR3nt:%s|Index:%s|Abundance:%s", 
                                   sample.id, 
                                   mouse.id,
                                   cdr3aa,
                                   cdr3nt,
                                   i,
                                   abundance)            
          }
          
          output.seq <- toString(unmasked(msaMiXCRtreeVDJ)[[i]])  
          # print(nchar(output.seq))
          cat(output.info)
          cat("\n")
          cat(output.seq)
          cat("\n")
        }
        sink()
        #####----------------------------------------------------------------#####
        # finished writing to fasta file
        #####----------------------------------------------------------------#####
        tmp.sumdf <- data.frame(
          MouseID = c(mouse.id),
          Case = c(yfp.case),
          MID = c(paste( yfp.mids[[mouse.id]][[yfp.case]], collapse = ",")),
          cloneID = c(input.VJ.combi),
          fastaName = c(basename(path.to.output.fasta)),
          numSeq = c(length(all.seqs))
        )
        summary.fastadf <- rbind(summary.fastadf, tmp.sumdf)
      }
    }
  }
}
write.csv(summary.fastadf, file.path(path.to.05.output, "summary_all_FASTA.csv"))


