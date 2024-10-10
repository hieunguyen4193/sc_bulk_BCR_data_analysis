gc()
rm(list = ls())

path.to.project.src <- "/home/hieu/src/BCRTree_release/gctree/data_analysis"
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
library(RColorBrewer)
library(scales)
library(ggpubr)
library(circlize)
library(ggalluvial)
outdir <- "/home/hieu/outdir"
PROJECT <- "mixcr_pipeline_output"
thres <- 0.15

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output", sprintf("CDR3_%s", thres))
path.to.04.output <- file.path(path.to.main.output, "04_output", sprintf("CDR3_%s", thres))
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

path.to.mid.metadata <- file.path(outdir, "FT_output", "mid_labels.csv")
mid.metadata <- read.csv(path.to.mid.metadata, sep = ";")
mid.metadata$MID <- mid.metadata$X
mid.metadata <- mid.metadata %>% rowwise() %>%
  mutate(X = sprintf("%s_%s", mouse, str_replace(str_replace(population, "Ly6c[+]", ""), "Ly6c[-]", "")))
count.mid.in.mouse <- table(mid.metadata$mouse)

for (mouse.id in unique(mid.metadata$mouse)){
  # mouse.id <- "m28"
  print(sprintf("Working on mouse %s", mouse.id))
  all.mids <- subset(mid.metadata, mid.metadata$mouse == mouse.id)$X %>% unique()
  clonedf <- readxl::read_excel(file.path(path.to.02.output, sprintf("clones_%s.xlsx", mouse.id))) %>%
    rowwise() %>%
    mutate(merge.samples = paste(to_vec( for(item in str_split(samples, ",")[[1]]) subset(mid.metadata, mid.metadata$MID == item)$X), collapse = ","))
  clonesetsdf <- readxl::read_excel(file.path(path.to.02.output, sprintf("clonesets_%s.split_clones_%s.xlsx", mouse.id, thres)))
  clonesetsdf <- clonesetsdf %>% rowwise() %>%
    mutate(merge.id = subset(mid.metadata, mid.metadata$MID == id)$X)
  ##### NEED FIX!
  #> 27.06.2024: currently using the combination of V and J gene and length of CDR3 
  #> AA as CloneID, this hasn't taken into account the similarity/difference between
  #> CDR3 sequences. ---> DONE
  clonesetsdf$CloneID <- clonesetsdf[[sprintf("VJcombi_CDR3_%s", thres)]]
  all.clones <- unique(clonesetsdf$CloneID)
  #####----------------------------------------------------------------------#####
  ##### bar plot: number of occurrences for each gene usage in each MID
  #####----------------------------------------------------------------------#####
  path.to.save.plot <- file.path(path.to.04.output, "barplot", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  # input.mid <- "MID4"
  for (input.mid in all.mids){
    if (file.exists(file.path(path.to.save.plot, sprintf("%s.barplot.svg", input.mid))) == FALSE){
      tmpdf <- subset(clonedf, grepl(input.mid, clonedf$merge.samples) == TRUE)
      count.V.genes <- table(tmpdf$V.gene) %>% as.data.frame() %>% arrange(desc(Freq))
      count.V.genes <- count.V.genes %>% rowwise() %>%
        mutate(pct = Freq/sum(count.V.genes$Freq))
      # Count
      p.count <- count.V.genes %>% ggplot(aes(x = reorder(Var1, Freq), y = Freq)) + 
        geom_bar(stat = "identity") +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 90)) +
        xlab("V gene") + ylab("Count") + 
        ggtitle(sprintf("Mouse %s, %s", mouse.id, input.mid))
      ggsave(plot = p.count, filename = sprintf("%s.barplot.svg", input.mid), 
             path = path.to.save.plot,
             device = "svg",
             width = 14, 
             height = 10)     
      # Percentage
      p.pct <- count.V.genes %>% ggplot(aes(x = reorder(Var1, pct), y = pct)) + 
        geom_bar(stat = "identity") +
        theme_pubr() +
        theme(axis.text.x = element_text(angle = 90)) +
        xlab("V gene") + ylab("Count") + 
        ggtitle(sprintf("Mouse %s, %s", mouse.id, input.mid))
      ggsave(plot = p.pct, filename = sprintf("%s.pct.barplot.svg", input.mid), 
             path = path.to.save.plot,
             device = "svg",
             width = 14, 
             height = 10)   
    }
  }
  
  #####----------------------------------------------------------------------#####
  ##### chord diagram
  #####----------------------------------------------------------------------#####
  path.to.save.plot <- file.path(path.to.04.output, "chord_diagram", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  if (length(all.mids) > 1){
    count.share.clones <- data.frame(MID = all.mids)
    for (input.mid in all.mids){
      clones1 <- subset(clonesetsdf, clonesetsdf$merge.id == input.mid)$CloneID %>% unique()
      count.share.clones[[input.mid]] <- unlist(lapply(
        count.share.clones$MID, function(x){
          clones <- subset(clonesetsdf, clonesetsdf$merge.id == x)$CloneID
          return(length(intersect(clones, clones1)))
        }
      ))
    }
    
    count.share.clones <- count.share.clones %>% column_to_rownames("MID") %>% as.matrix()
    diag(count.share.clones) <- 0
    svg(file.path(path.to.save.plot, sprintf("%s.svg", mouse.id)), width = 10, height = 10) 
    chordDiagram(count.share.clones, transparency = 0.8) 
    dev.off()
    circos.clear()
  }
  
  
  #####----------------------------------------------------------------------#####
  ##### Morisita-Horn index to measure the similarity between two groups of clones
  #####----------------------------------------------------------------------#####
  path.to.save.plot <- file.path(path.to.04.output, "MHI_heatmap", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  if (length(all.mids) > 1){
    mhidf <- data.frame(MID = all.mids)
    for (input.mid1 in all.mids){
      mhidf[[input.mid1]] <- unlist(lapply(
        mhidf$MID, function(input.mid2){
          S <- subset(clonesetsdf, clonesetsdf$merge.id %in% c(input.mid1, input.mid2))$CloneID %>% unique() %>% length()
          X <- subset(clonesetsdf, clonesetsdf$merge.id == input.mid1) %>% nrow()
          Y <- subset(clonesetsdf, clonesetsdf$merge.id == input.mid2) %>% nrow()
          
          x <- unlist(lapply(
            subset(clonesetsdf, clonesetsdf$merge.id %in% c(input.mid1, input.mid2))$CloneID %>% unique(),
            function(x){
              return(subset(clonesetsdf, clonesetsdf$merge.id == input.mid1 & clonesetsdf$CloneID == x) %>% nrow())
            }
          ))
          y <- unlist(lapply(
            subset(clonesetsdf, clonesetsdf$merge.id %in% c(input.mid1, input.mid2))$CloneID %>% unique(),
            function(x){
              return(subset(clonesetsdf, clonesetsdf$merge.id == input.mid2 & clonesetsdf$CloneID == x) %>% nrow())
            }
          ))
          nom <- 2 * sum(x * y)
          det <- (sum(x^2)/X^2) + (sum(y^2)/Y^2)
          mhi <- nom/(X*Y * det)
          return(mhi)
        }
      ))
    }
    
    svg(file.path(path.to.save.plot, sprintf("%s.svg", mouse.id)), width = 10, height = 10) 
    mhidf %>% column_to_rownames("MID") %>% pheatmap::pheatmap()
    dev.off()
  }
  
  #####----------------------------------------------------------------------#####
  path.to.save.plot <- file.path(path.to.04.output, "alluvial_tracking", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  if (length(all.mids) > 1){
    for (input.mid in all.mids){
      if (file.exists(file.path(path.to.save.plot, sprintf("%s.tracking.svg", input.mid))) == FALSE){
        mid.clones <- subset(clonesetsdf, clonesetsdf$merge.id == input.mid)
        mid.clones <- table(mid.clones$CloneID) %>% sort() %>% tail(10) %>% names()
        clone.in.mid <- data.frame(MID = all.mids)
        for (cloneid in mid.clones){
          clone.in.mid[[cloneid]] <- unlist(lapply(
            all.mids, function(x){
              (subset(clonesetsdf, clonesetsdf$merge.id == x & clonesetsdf$CloneID == cloneid) %>% nrow())/nrow(clonesetsdf)
            }
          ))
        }
        
        clone.ind.mid.pivot <- clone.in.mid %>% pivot_longer(!MID, names_to = "Clone", values_to = "Count")
        
        is_lodes_form(clone.ind.mid.pivot, key = MID, value = Count, id = Clone, silent = TRUE)
        # clone.ind.mid.pivot$Clone <- factor(clone.ind.mid.pivot$Clone, levels = unique(clone.ind.mid.pivot$Clone))
        p <- ggplot(clone.ind.mid.pivot,
                    aes(x = MID, stratum = Clone, alluvium = Clone, y = Count, fill = Clone)) +
          geom_flow(stat = "alluvium", lode.guidance = "frontback",
                    color = "darkgray") +
          geom_stratum() +
          theme(legend.position = "bottom") 
        
        ggsave(plot = p, filename = sprintf("%s.tracking.svg", input.mid), 
               path = path.to.save.plot,
               device = "svg",
               width = 14, 
               height = 10) 
      }
    }
  }
}

