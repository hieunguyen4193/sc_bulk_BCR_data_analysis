new.pkgs <- c("svglite", "ggpubr")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }   
}
if (packageVersion("tidyr") != "1.3.1"){
  install.packages("tidyr")  
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
library(ggpubr)
library(circlize)
library(ggalluvial)

`%ni%` = Negate(`%in%`)


#####----------------------------------------------------------------------#####
##### function to create an interactive data table
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

#####----------------------------------------------------------------------#####
##### Compute distance matrix for each pair of AA sequences
#####----------------------------------------------------------------------#####
compute_distance_matrix <- function(seqs) {
  n <- length(seqs)
  dist_matrix <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist <- stringdist(seqs[i], seqs[j], method = "lv")
      dist_matrix[i, j] <- dist
      dist_matrix[j, i] <- dist
    }
  }
  return(dist_matrix)
}

assign_clusters_to_sequences <- function(seqs, threshold = 0.15){
  distance_matrix <- compute_distance_matrix(seqs)
  max_length <- max(nchar(seqs))
  normalized_matrix <- distance_matrix / max_length
  graph <- graph_from_adjacency_matrix(normalized_matrix < threshold, mode = "undirected", diag = FALSE)
  cl <- cluster_optimal(graph)
  clus.res <- data.frame(seq = seqs, cluster = cl$membership)
  # function to plot the graph and its node community
  # plot(cl, graph, layout=layout_with_fr(graph))
  return(list(
    graph = graph,
    cl = cl, 
    res = clus.res
  ))
}

#####----------------------------------------------------------------------#####
##### FUNCTIONS to preprocess the bulk VDJ data
#####----------------------------------------------------------------------#####
run_preprocessing_bulk_VDJ_data <- function(outdir, PROJECT, thres, ref.gene){
  path.to.02.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
  dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)
  path.to.fasta <- file.path(path.to.main.src, "FASTA", ref.gene, "V_J_genes")
  path.to.storage <- "/media/hieunguyen/GSHD_HN01/storage/all_BSimons_datasets"
  
  # input files from mixcr pipeline output, stored in storage.
  path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output", "mid_based_output")
  
  all.clones.tables <- Sys.glob(file.path(path.to.mid.output, "*", "*.reassigned_IGH.tsv"))
  names(all.clones.tables) <- unlist(lapply(
    basename(all.clones.tables), function(x){
      return(str_replace(x, ".reassigned_IGH.tsv", ""))
    }
  ))
  
  #####----------------------------------------------------------------------#####
  ##### read REFERENCE GENES
  #####----------------------------------------------------------------------#####
  s.V.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHV.fasta")) 
  names(s.V.genes) <- lapply(names(s.V.genes), function(x){
    str_split(x, "[|]")[[1]][[2]]
  })
  
  s.J.genes <- readDNAStringSet(file.path(path.to.fasta, "IGHJ.fasta")) 
  names(s.J.genes) <- lapply(names(s.J.genes), function(x){
    str_split(x, "[|]")[[1]][[2]]
  })
  
  #####----------------------------------------------------------------------#####
  ##### READ DATA FROM ALL CLONE TABLES IN THE INPUT DIR
  #####----------------------------------------------------------------------#####
  cloneset_files <- all.clones.tables
  
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
  
  #####----------------------------------------------------------------------#####
  ##### Group sequences + Gene usages to clones
  #####----------------------------------------------------------------------#####
  new.clonesets <- data.frame()
  for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
    tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
    seqs <- unique(tmpdf$aaSeqCDR3)
    if (length(seqs) >= 2){
      cluster.output <- assign_clusters_to_sequences(seqs, threshold = thres)$res
      tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
        tmpdf$aaSeqCDR3, function(x){
          return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
        }
      ))    
    } else {
      tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
    }
    
    new.clonesets <- rbind(new.clonesets, tmpdf)
  }
  writexl::write_xlsx(new.clonesets, file.path(path.to.02.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT)))
  #####----------------------------------------------------------------------#####
  ##### Create the final clone dataframe
  #####----------------------------------------------------------------------#####
  clonedf <- data.frame(VJseq.combi.tmp = unique(new.clonesets$VJseq.combi)) %>%
    rowwise() %>%
    mutate(CDR3aa = str_split(VJseq.combi.tmp, "_")[[1]][[3]]) %>%
    mutate(V.gene = str_split(VJseq.combi.tmp, "_")[[1]][[1]]) %>%
    mutate(J.gene = str_split(VJseq.combi.tmp, "_")[[1]][[2]]) %>%
    mutate(CDR3nt = str_split(VJseq.combi.tmp, "_")[[1]][[4]]) %>%
    mutate(cloneSize = nrow(subset(new.clonesets, new.clonesets$VJseq.combi == VJseq.combi.tmp))) %>% 
    mutate(CDR3aa.length = nchar(CDR3aa)) %>% 
    mutate(CDR3nt.length = nchar(CDR3nt)) %>% 
    mutate(samples = paste(unique(subset(new.clonesets, new.clonesets$VJseq.combi == VJseq.combi.tmp)$id), collapse = ",")) %>%
    mutate(uniqueMoleculeCount = subset(new.clonesets, new.clonesets$VJseq.combi == VJseq.combi.tmp)$uniqueMoleculeCount %>% sum()) %>%
    arrange(desc(cloneSize))
  clonedf <- subset(clonedf, select = -c(VJseq.combi.tmp))
  
  writexl::write_xlsx(clonedf, file.path(path.to.02.output, sprintf("%s.clonedf.xlsx", PROJECT)))
}