new.pkgs <- c("svglite", "Matrix", "car")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }   
}

# if (packageVersion("Matrix") != "1.5.4.1" | packageVersion("ggplot2") != "3.4.4"){
#   install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", type = "source", repos = NULL)
#   install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", type = "source", repos = NULL)  
# }

# if (packageVersion("igraph") != "2.1.1"){
#   install.packages("https://cran.r-project.org/src/contrib/igraph_2.1.1.tar.gz", type = "source", repos = NULL)
# }
if (packageVersion("ggplot2") != "3.4.4"){
  install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", type = "source", repos = NULL)
}

if ("msa" %in% installed.packages() == FALSE){
  BiocManager::install("msa", update = FALSE)
}

if (packageVersion("tidyr") != "1.3.1"){
  install.packages("tidyr")  
  install.packages("tidyselect")
  library(tidyr)
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
formatBulkVDJtable <- function(cloneset_files, 
                               PROJECT,
                               thres,
                               thres.dis,
                               savefile, 
                               path.to.save.output,
                               define.clone.clusters,
                               rerun){
  if (file.exists(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT))) == FALSE |
      file.exists(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT))) == FALSE | rerun == TRUE){
      
    print("Generate File clonedf.xlsx and flie split_clones.xlsx ...")
    clonesets <- read_tsv(cloneset_files, id="fileName") %>% 
      mutate(id=str_remove_all(fileName,".*\\/|.reassigned_IGH.tsv"),
             isotype=str_remove(allCHitsWithScore, "\\*.*")) 
    
    clonesets <- clonesets %>% rowwise() %>%
      mutate(V.gene = str_split(str_split(allVHitsWithScore, ",")[[1]][[1]], "[()]")[[1]][[1]]) %>%
      mutate(V.gene = paste(str_split(V.gene, "-")[[1]][1:2], collapse = "-")) %>% 
      # mutate(V.gene = ifelse(grepl("*", V.gene) == TRUE, str_split(V.gene, "[*]")[[1]][[1]] , V.gene)) %>%
      mutate(J.gene = str_split(str_split(allJHitsWithScore, ",")[[1]][[1]], "[()]")[[1]][[1]]) %>%
      mutate(J.gene = str_split(J.gene, "-")[[1]][[1]]) %>%
      # mutate(J.gene = ifelse(grepl("*", J.gene) == TRUE, str_split(J.gene, "[*]")[[1]][[1]] , J.gene)) %>%
      mutate(D.gene = str_split(str_split(allDHitsWithScore, ",")[[1]][[1]], "[()]")[[1]][[1]]) %>%
      mutate(D.gene = paste(str_split(D.gene, "-")[[1]][1:2], collapse = "-")) %>% 
      mutate(VJseq.combi = sprintf("%s_%s_%s_%s", V.gene, J.gene, aaSeqCDR3, nSeqCDR3)) %>%
      mutate(VJ.combi = sprintf("%s_%s", V.gene, J.gene)) %>%
      mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)))
    
    if (define.clone.clusters == TRUE){
      ##### Group sequences + Gene usages to clones
      ##### Definition: A clone = V + J gene, length of the CDR3 sequence.
      ##### use similarity distance to define sequences in a clone. 
      new.clonesets <- data.frame()
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
        new.clonesets <- rbind(new.clonesets, tmpdf)
      }
    } else {
      new.clonesets <- clonesets
    }
    
    ##### Create the final clone dataframe
    ##### definition: A clone = V + J gene + Sequence. 
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
    if (savefile == TRUE){
      writexl::write_xlsx(clonedf, file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT)))
      writexl::write_xlsx(new.clonesets, file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT)))    
    }
  } else {
    print("File clonedf.xlsx and flie split_clones.xlsx exist, reading in ...")
    clonedf <- readxl::read_excel(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT)))
    new.clonesets <- readxl::read_excel(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT))) 
  }
  
  output = list(
    clonedf = clonedf, 
    clonesets = new.clonesets
  )
  return(output)
}

#####----------------------------------------------------------------------#####
##### get all *.reassigned_IGH.tsv from preprocessed_files folder, format tables
#####----------------------------------------------------------------------#####
run_preprocessing_all_bulk_VDJ_data <- function(path.to.mid.output,
                                                path.to.save.output,
                                                PROJECT,
                                                thres, 
                                                thres.dis,
                                                savefile,
                                                verbose,
                                                rerun,
                                                define.clone.clusters){
  if (rerun == TRUE){
    ##### CLEAN UP OLD RESULTS
    if (file.exists(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT))) == TRUE){
      file.remove(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT)))
      print(sprintf("remove old data at %s", file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT))))
    }
    if (file.exists(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT))) == TRUE){
      file.remove(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT)))
      print(sprintf("remove old data at %s", file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT))))
    }
  }
  if (verbose == TRUE){
    print(sprintf("Output will be saved at %s", path.to.save.output))
  }
  if (verbose == TRUE){
    print(sprintf("fetching files *.reassigned_IGH.tsv from %s", path.to.mid.output)) 
  }
  all.clones.tables <- Sys.glob(file.path(path.to.mid.output, "*", "*.reassigned_IGH.tsv"))
  if (verbose == TRUE){
    print("reading and processing all clone tables ...") 
    print(sprintf("Number of files in the input path %s, %s", path.to.mid.output, length(all.clones.tables)))
  }
  names(all.clones.tables) <- unlist(lapply(
    basename(all.clones.tables), function(x){
      return(str_replace(x, ".reassigned_IGH.tsv", ""))
    }
  ))
  if (file.exists(file.path(path.to.save.output, sprintf("%s.rds", PROJECT))) == FALSE | rerun == TRUE){
    if (file.exists(file.path(path.to.save.output, sprintf("%s.rds", PROJECT))) == TRUE){
      file.remove(file.path(path.to.save.output, sprintf("%s.rds", PROJECT)))
    }
    #> by running this function, we read all the clone table "*.reassigned_IGH.tsv"
    #> from the input path.to.storage, path.to.mid.output and preprocess them all.
    #> Generate one final big VDJ data table containing all clones from all samples / all MIDs
    #> and all mice
    if (verbose == TRUE){
      print("format the table to the unified VDJ table format...") 
    }
    output <- formatBulkVDJtable(cloneset_files = all.clones.tables, 
                                 PROJECT = PROJECT,
                                 thres = thres,
                                 thres.dis = thres.dis,
                                 savefile = savefile,
                                 path.to.save.output = path.to.save.output,
                                 rerun = rerun,
                                 define.clone.clusters = define.clone.clusters)
    saveRDS(output, file.path(path.to.save.output, sprintf("%s.rds", PROJECT)))
  } else {
    output <- readRDS(file.path(path.to.save.output, sprintf("%s.rds", PROJECT)))
  }
  if (verbose == TRUE){
   print("Finished!") 
  }
  return(output)
}

#####----------------------------------------------------------------------#####
##### FUNCTIONS to preprocess the single cell VDJ data
#####----------------------------------------------------------------------#####
formatScVDJtable <- function(all.VDJ.files,
                             PROJECT, 
                             thres,
                             thres.dis,
                             savefile, 
                             path.to.save.output,
                             rerun,
                             define.clone.clusters){
  if (file.exists(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT))) == FALSE |
      file.exists(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT))) == FALSE |
      rerun == TRUE){
    print("Generate File clonedf.xlsx and flie split_clones.xlsx ...")
    clonesets <- data.frame()
    for (file in all.VDJ.files){
      print(file)
      sampleid <- basename(file) %>% str_replace(".xlsx", "")
      tmpdf <- readxl::read_xlsx(file)  
      tmpdf$id <- sampleid
      clonesets <- rbind(clonesets, tmpdf)
    }
    
    # Keep IGH chain only, to compare with the BULK data
    clonesets <- subset(clonesets, clonesets$chain == "IGH") %>%
      rowwise() %>%
      mutate(VJseq.combi = sprintf("%s_%s_%s_%s", v_gene, j_gene, cdr3, cdr3_nt)) %>%
      mutate(V.gene = v_gene) %>%
      # mutate(V.gene = ifelse(grepl("*", V.gene) == TRUE, str_split(V.gene, "[*]")[[1]][[1]] , V.gene)) %>%
      mutate(J.gene = j_gene) %>%
      # mutate(J.gene = ifelse(grepl("*", J.gene) == TRUE, str_split(J.gene, "[*]")[[1]][[1]] , J.gene)) %>%
      mutate(D.gene = d_gene) %>% 
      mutate(nSeqCDR1 = cdr1_nt) %>%
      mutate(nSeqCDR2 = cdr2_nt) %>%
      mutate(nSeqCDR3 = cdr3_nt) %>%
      mutate(aaSeqCDR1 = cdr1) %>%
      mutate(aaSeqCDR2 = cdr2) %>%
      mutate(aaSeqCDR3 = cdr3) %>%
      mutate(nSeqFR1 = fwr1_nt) %>%
      mutate(nSeqFR2 = fwr2_nt) %>%
      mutate(nSeqFR3 = fwr3_nt) %>%
      mutate(nSeqFR4 = fwr4_nt) %>%
      mutate(aaSeqFR1 = fwr1) %>%
      mutate(aaSeqFR2 = fwr2) %>%
      mutate(aaSeqFR3 = fwr3) %>%
      mutate(aaSeqFR4 = fwr4) %>%
      mutate(VJseq.combi = sprintf("%s_%s_%s_%s", V.gene, J.gene, aaSeqCDR3, nSeqCDR3)) %>%
      mutate(VJ.combi = sprintf("%s_%s", V.gene, J.gene)) %>%
      mutate(VJ.len.combi = sprintf("%s_%s_%s", V.gene, J.gene,  nchar(nSeqCDR3)))
    clonesets$uniqueMoleculeCount <- NA
    clonesets <- clonesets %>% rowwise() %>%
      mutate(targetSequences = paste0(c(
        nSeqFR1,
        nSeqCDR1,
        nSeqFR2,
        nSeqCDR2,
        nSeqFR3,
        nSeqCDR3,
        nSeqFR4
      ), collapse = ""))
    if (define.clone.clusters == TRUE){
      ##### Group sequences + Gene usages to clones
      new.clonesets <- data.frame()
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
        new.clonesets <- rbind(new.clonesets, tmpdf)
      }
    } else {
      new.clonesets <- clonesets
    }
    
    clonedf <- data.frame(
      VJseq.combi.tmp = unique(new.clonesets$VJseq.combi)) %>%
      rowwise() %>%
      mutate(CDR3aa = str_split(VJseq.combi.tmp, "_")[[1]][[3]]) %>%
      mutate(V.gene = str_split(VJseq.combi.tmp, "_")[[1]][[1]]) %>%
      mutate(J.gene = str_split(VJseq.combi.tmp, "_")[[1]][[2]]) %>%
      mutate(CDR3nt = str_split(VJseq.combi.tmp, "_")[[1]][[4]]) %>%
      mutate(cloneSize = nrow(subset(new.clonesets, new.clonesets$VJseq.combi == VJseq.combi.tmp))) %>% 
      mutate(CDR3aa.length = nchar(CDR3aa)) %>% 
      mutate(CDR3nt.length = nchar(CDR3nt)) %>% 
      mutate(samples = paste(unique(subset(new.clonesets, new.clonesets$VJseq.combi == VJseq.combi.tmp)$id), collapse = ",")
      )
    
    clonedf <- subset(clonedf, select = -c(VJseq.combi.tmp))
    if (savefile == TRUE){
      writexl::write_xlsx(clonedf, file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT)))
      writexl::write_xlsx(new.clonesets, file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT)))
    }    
  } else {
    print("File clonedf.xlsx and flie split_clones.xlsx exist, reading in ...")
    clonedf <- readxl::read_excel(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT)))
    new.clonesets <- readxl::read_excel(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT)))
  }
  output = list(
    clonedf = clonedf, 
    clonesets = new.clonesets
  )
  return(output)
}

#####----------------------------------------------------------------------#####
##### Function to get all VDJ input files and format them all into VDJ tables
#####----------------------------------------------------------------------#####
run_preprocessing_all_sc_data <- function(path.to.VDJ.output, 
                                          path.to.save.output, 
                                          PROJECT,
                                          thres, 
                                          thres.dis,
                                          savefile,
                                          rerun,
                                          define.clone.clusters){
  if (rerun == TRUE){
    ##### CLEAN UP OLD RESULTS
    if (file.exists(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT))) == TRUE){
      file.remove(file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT)))
      print(sprintf("Remove old data at %s", file.path(path.to.save.output, sprintf("%s.clonedf.xlsx", PROJECT))))
    }
    if (file.exists(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT))) == TRUE){
      file.remove(file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT)))
      print(sprintf("Remove old data at %s", file.path(path.to.save.output, sprintf("clonesets_%s.split_clones.xlsx", PROJECT))))
    }
  }
  dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
  all.VDJ.files <- Sys.glob(file.path(path.to.VDJ.output, "*.xlsx"))
  names(all.VDJ.files) <- unlist(lapply(all.VDJ.files, function(x){
    str_replace(basename(x), ".xlsx", "")
  }))
  # temp remove 
  all.VDJ.files <- all.VDJ.files[names(all.VDJ.files) != "GF_7w_14w"]
  if (file.exists(file.path(path.to.save.output, sprintf("%s.rds", PROJECT))) == FALSE | rerun == TRUE){
    if (file.exists(file.path(path.to.save.output, sprintf("%s.rds", PROJECT))) == TRUE){
      file.remove(file.path(path.to.save.output, sprintf("%s.rds", PROJECT)))
    }
    print(sprintf("Generating file %s", file.path(path.to.save.output, sprintf("%s.rds", PROJECT))))
    output <- formatScVDJtable(all.VDJ.files = all.VDJ.files,
                               PROJECT = PROJECT, 
                               thres = thres,
                               thres.dis = thres.dis,
                               savefile = savefile, 
                               path.to.save.output = path.to.save.output,
                               rerun = rerun,
                               define.clone.clusters = define.clone.clusters)
    saveRDS(output, file.path(path.to.save.output, sprintf("%s.rds", PROJECT)))
  } else { 
    print(sprintf("File %s exists", file.path(path.to.save.output, sprintf("%s.rds", PROJECT))))
    output <- readRDS(file.path(path.to.save.output, sprintf("%s.rds", PROJECT)))
  }
  return(output)
}

#####----------------------------------------------------------------------#####
##### Generate FASTA files from clonesets
#####----------------------------------------------------------------------#####
generate_fasta <- function(clonesets, 
                           mouse.id,
                           path.to.save.output, 
                           ref.gene, 
                           ref.gene.config,
                           PROJECT,
                           thres = 0.85,
                           thres.dis = 0.15,
                           save_fasta = TRUE,
                           re_define_clone_cluster = FALSE){
    dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
    if (re_define_clone_cluster == TRUE){
      print("RE-DEFINE THE CLONE CLUSTERS BASED ON CDR3 SEQUENCE SIMILARITY AND V-J GENE USAGES ON SELECTED SAMPLES/MIDS ONLY")
      print(sprintf("Rename the column %s to %s", 
                    sprintf("VJcombi_CDR3_%s", thres),
                    sprintf("VJcombi_CDR3_%s_prev", thres)))
      print(sprintf("Remove the column %s", sprintf("VJcombi_CDR3_%s", thres)))
      clonesets[[sprintf("VJcombi_CDR3_%s_prev", thres)]] <- clonesets[[sprintf("VJcombi_CDR3_%s", thres)]]
      clonesets[[sprintf("VJcombi_CDR3_%s", thres)]] <- NULL
      
      ##### Group sequences + Gene usages to clones
      new.clonesets <- data.frame()
      for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
        tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
        seqs <- unique(tmpdf$aaSeqCDR3)
        print(sprintf("VJ.len.combi: %s, num seqs: %s", input.VJ.combi, length(seqs)))
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
        new.clonesets <- rbind(new.clonesets, tmpdf)
      }
      clonesets <- new.clonesets
      print("IMPORTANT NOTE: THE CLONE CLUSTERS ARE GENERATED BASED ON SAMPLES/MIDS IN THIS SELECTED MOUSE/SAMPLE ONLY; NOT THE AGGREGATED TABLE.")
      print(sprintf("Saving new cloneset files after assign clones to clusters to %s", file.path(path.to.save.output, "clonesets.csv")))
      write.csv(clonesets, file.path(path.to.save.output, "clonesets.csv"))
    }
    source(ref.gene.config) # path to the configuration file for the input BCR reference genes
    if (ref.gene == "10x"){
      ref.fasta <- readDNAStringSet(ref.genes$`10x`)
    } else if (ref.gene == "IMGT"){
      s.V.genes <- readDNAStringSet(ref.genes$IMGT$V.gene)
      names(s.V.genes) <- lapply(names(s.V.genes), function(x){
        x <- str_split(x, "[|]")[[1]][[2]]
        # x <- str_split(x, "[*]")[[1]][[1]]
        return(x)
      })
      s.J.genes <- readDNAStringSet(ref.genes$IMGT$J.gene)
      names(s.J.genes) <- lapply(names(s.J.genes), function(x){
        x <- str_split(x, "[|]")[[1]][[2]]
        # x <- str_split(x, "[*]")[[1]][[1]]
        return(x)
      })
    }
    for (input.VJ.combi in unique(clonesets[[sprintf("VJcombi_CDR3_%s", thres)]])){
      V.gene <- str_split(input.VJ.combi, "_")[[1]][[1]]
      J.gene <- str_split(input.VJ.combi, "_")[[1]][[2]]
      CDR3.length <- as.numeric(str_split(input.VJ.combi , "_")[[1]][[3]])
      # remove the * sign in the file name
      path.to.fasta.file <- file.path(path.to.save.output, 
                                      sprintf("%s.fasta", str_replace_all(input.VJ.combi, "[*]", "-")))
      
      if (file.exists(path.to.fasta.file) == FALSE){
        fasta.output <- subset(clonesets, clonesets[[sprintf("VJcombi_CDR3_%s", thres)]] == input.VJ.combi)[, c("targetSequences", 
                                                                                                                "uniqueMoleculeCount", 
                                                                                                                "V.gene", 
                                                                                                                "D.gene", 
                                                                                                                "J.gene", 
                                                                                                                "id", 
                                                                                                                "aaSeqCDR3", 
                                                                                                                "nSeqCDR3")]
        colnames(fasta.output) <- c("seq", 
                                    "abundance", 
                                    "V.gene", 
                                    "D.gene", 
                                    "J.gene", 
                                    "SampleID",
                                    "CDR3aa", 
                                    "CDR3nt")
        
        ##### get germline sequences and merge with the real data sequences
        if (ref.gene == "IMGT"){
          GL.V.gene <- s.V.genes[[V.gene]] %>% as.character()
          GL.J.gene <- s.J.genes[[J.gene]] %>% as.character()
        } else if (ref.gene == "10x"){
          GL.V.gene <- c[[V.gene]] %>% as.character()
          GL.J.gene <- ref.fasta[[J.gene]] %>% as.character()
        }
        
        repN.seq <- paste(replicate(n = CDR3.length, expr = "N"), collapse = "")
        GL.seq <- sprintf("%s%s%s", GL.V.gene, repN.seq, GL.J.gene)
        # merge the clone sequences with the reference sequence. 
        all.seqs <- c(fasta.output %>% pull(`seq`), GL.seq)
        
        if (nrow(fasta.output) > 1){
          ##### multiple alignment sequences, package MSA. 
          MiXCRtreeVDJ <- all.seqs %>% DNAStringSet()
          msaMiXCRtreeVDJ <- msa(inputSeqs = MiXCRtreeVDJ, verbose = TRUE)
          if (save_fasta == TRUE){
            print(sprintf("File fasta: %s", basename(path.to.fasta.file)))
            sink(path.to.fasta.file)
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
          }
        } 
      }
    }
} 

##### thresholds
# low mutated 0-3 
# intermediate 4-15
# highly mutated >15 
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
