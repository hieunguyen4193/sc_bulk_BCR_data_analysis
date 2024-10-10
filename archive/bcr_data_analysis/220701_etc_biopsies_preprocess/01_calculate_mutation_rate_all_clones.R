gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/bcr_data_analysis/220701_etc_biopsies_preprocess"
source(file.path(path.to.project.src, "import_libraries.R"))
source(file.path(path.to.project.src, "helper_functions.R"))

new.pkgs <- c("msa")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() ==  FALSE){
    BiocManager::install(pkg, update = FALSE)
  }
}
library(tidyr)
library(readr)
library(stringr)
library(stringdist)
library(ggplot2)
library(msa)
library(magrittr)
library(ggbeeswarm)
library(tibble)
library(dplyr)
library(plotly)
library(ggtree)
library(ggrepel)
library(data.table)
library(ggtree)
library(ape)
library(igraph)

outdir <- "/media/hieunguyen/HNSD01/outdir"
PROJECT <- "220701_etc_biopsies"
thres <- 0.85

path.to.main.input <- file.path(outdir, PROJECT, "mixcr_pipeline_output")
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.01.output <- file.path(path.to.main.output, "01_output", sprintf("CDR3_%s", thres))
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.mid.output <- file.path(path.to.main.input, "mid_based_output")

##### Get germline sequences for V, D, J genes
source(file.path(path.to.project.src, "get_germline_sequences.R"))

##### METADATA
path.to.mid.metadata <- file.path(path.to.project.src, "mid_labels.csv")
mid.metadata <- read.csv(path.to.mid.metadata, sep = ";")

mid.metadata$timepoint <- mid.metadata$population

#####----------------------------------------------------------------------#####
##### Get clone information dataframes
#####----------------------------------------------------------------------#####
all.clones.tables <- Sys.glob(file.path(path.to.mid.output, "*", "*.reassigned_IGH.tsv"))
names(all.clones.tables) <- unlist(lapply(
  basename(all.clones.tables), function(x){
    return(str_replace(x, ".reassigned_IGH.tsv", ""))
  }
))

for (mouse.id in unique(mid.metadata$mouse)){
  dir.create(file.path(path.to.01.output, mouse.id), showWarnings = FALSE, recursive = TRUE)
  print(sprintf("generating clone dataframe for mouse %s", mouse.id))
  all.mids <- subset(mid.metadata, mid.metadata$mouse == mouse.id)$X
  
  cloneset_files <- all.clones.tables[all.mids]
  
  clonesets <- read_tsv(cloneset_files, id="fileName") %>% 
    mutate(id=str_remove_all(fileName,".*\\/|.reassigned_IGH.tsv"),
           isotype=str_remove(allCHitsWithScore, "\\*.*")) 
  
  clonesets <- clonesets %>% rowwise() %>%
    mutate(V.gene = str_split(str_split(allVHitsWithScore, ",")[[1]][[1]], "[()]")[[1]][[1]]) %>%
    # mutate(V.gene = str_split(V.gene, "[*]")[[1]][[1]]) %>% # do not remove * <<<<<
    mutate(V.gene = paste(str_split(V.gene, "-")[[1]][1:2], collapse = "-")) %>% 
    mutate(J.gene = str_split(str_split(allJHitsWithScore, ",")[[1]][[1]], "[()]")[[1]][[1]]) %>%
    # mutate(J.gene = str_split(J.gene, "[*]")[[1]][[1]]) %>% # do not remove * <<<<<
    mutate(J.gene = str_split(J.gene, "-")[[1]][[1]]) %>%
    mutate(D.gene = str_split(str_split(allDHitsWithScore, ",")[[1]][[1]], "[()]")[[1]][[1]]) %>%
    mutate(D.gene = paste(str_split(D.gene, "-")[[1]][1:2], collapse = "-")) %>% 
    mutate(VJseq.combi = sprintf("%s_%s_%s_%s", V.gene, J.gene, aaSeqCDR3, nSeqCDR3)) %>%
    mutate(VJ.combi = sprintf("%s_%s", V.gene, J.gene)) %>%
    mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)))
  writexl::write_xlsx(clonesets, file.path(path.to.01.output, sprintf("clonesets_%s.xlsx", mouse.id)))
  
  #####---------------------------------------------------------------------#####
  ##### Group AA sequences to group based on 85% similarity
  #####---------------------------------------------------------------------#####    
  clonesetsdf <- readxl::read_excel(file.path(path.to.01.output, sprintf("clonesets_%s.xlsx", mouse.id)))
  new.clonesetsdf <- data.frame()
  
  for (input.VJ.combi in unique(clonesetsdf$VJ.len.combi)){
    tmpdf <- subset(clonesetsdf, clonesetsdf$VJ.len.combi == input.VJ.combi)
    seqs <- unique(tmpdf$aaSeqCDR3)
    if (length(seqs) >= 2){
      cluster.output <- assign_clusters_to_sequences(seqs, threshold = 1- thres)$res
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
  writexl::write_xlsx(new.clonesetsdf, file.path(path.to.01.output, sprintf("clonesets_%s.split_clones_%s.xlsx", mouse.id, thres)))
  
  # replace the clonesets by new.clonesets for downstream analysis
  clonesets <- new.clonesetsdf
  #####---------------------------------------------------------------------#####
  ##### Generate summary table for clones, not for the fasta files
  #####---------------------------------------------------------------------#####
  if (file.exists(file.path(path.to.01.output, sprintf("clones_%s.xlsx", mouse.id))) == FALSE){
    clonedf <- data.frame(VJseq.combi.tmp = unique(clonesets$VJseq.combi)) %>%
      rowwise() %>%
      mutate(CDR3aa = str_split(VJseq.combi.tmp, "_")[[1]][[3]]) %>%
      mutate(V.gene = str_split(VJseq.combi.tmp, "_")[[1]][[1]]) %>%
      mutate(J.gene = str_split(VJseq.combi.tmp, "_")[[1]][[2]]) %>%
      mutate(CDR3nt = str_split(VJseq.combi.tmp, "_")[[1]][[4]]) %>%
      mutate(cloneSize = nrow(subset(clonesets, clonesets$VJseq.combi == VJseq.combi.tmp))) %>% 
      mutate(CDR3aa.length = nchar(CDR3aa)) %>% 
      mutate(CDR3nt.length = nchar(CDR3nt)) %>% 
      mutate(samples = paste(unique(subset(clonesets, clonesets$VJseq.combi == VJseq.combi.tmp)$id), collapse = ",")) %>%
      mutate(uniqueMoleculeCount = subset(clonesets, clonesets$VJseq.combi == VJseq.combi.tmp)$uniqueMoleculeCount %>% sum()) %>%
      arrange(desc(cloneSize))
    
    clonedf <- subset(clonedf, select = -c(VJseq.combi.tmp))
    
    writexl::write_xlsx(clonedf, file.path(path.to.01.output, sprintf("clones_%s.xlsx", mouse.id)))
  } else {
    clonedf <-  readxl::read_excel(file.path(path.to.01.output, sprintf("clones_%s.xlsx", mouse.id)))
  }
}