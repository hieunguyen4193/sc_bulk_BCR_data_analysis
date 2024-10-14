gc()
rm(list = ls())

source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/00_helper_functions.R")
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis"
PROJECT <- "240826_BSimons"
thres <- 0.85

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.VDJ.input <- file.path(path.to.main.output, "01_output", sprintf("CDR3_%s", thres))
path.to.01.output <- file.path(path.to.main.output, "02_output", sprintf("CDR3_%s", thres))
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.mid.metadata <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240826_BSimons/240829 sample sheet.xlsx"

if (grepl(".xlsx", path.to.mid.metadata) == TRUE){
  mid.metadata <- readxl::read_excel(path.to.mid.metadata)
} else {
  mid.metadata <- read.csv(path.to.mid.metadata)
}

count.mid.in.mouse <- table(mid.metadata$mouse)
for (mouse.id in unique(mid.metadata$mouse)){
  # mouse.id <- "m28"
  print(sprintf("Working on mouse %s", mouse.id))
  all.mids <- subset(mid.metadata, mid.metadata$mouse == mouse.id)$MID %>% unique()
  clonedf <- readxl::read_excel(file.path(path.to.VDJ.input, sprintf("clones_%s.xlsx", mouse.id)))
  clonesetsdf <- readxl::read_excel(file.path(path.to.VDJ.input, sprintf("clonesets_%s.split_clones_%s.xlsx", mouse.id, thres)))
  
  ##### NEED FIX!
  #> 27.06.2024: currently using the combination of V and J gene and length of CDR3 
  #> AA as CloneID, this hasn't taken into account the similarity/difference between
  #> CDR3 sequences. ---> DONE
  clonesetsdf$CloneID <- clonesetsdf[[sprintf("VJcombi_CDR3_%s", thres)]]
  all.clones <- unique(clonesetsdf$CloneID)
  #####----------------------------------------------------------------------#####
  ##### bar plot: number of occurrences for each gene usage in each MID
  #####----------------------------------------------------------------------#####
  path.to.save.plot <- file.path(path.to.01.output, "barplot", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  # input.mid <- "MID4"
  for (input.mid in all.mids){
    if (file.exists(file.path(path.to.save.plot, sprintf("%s.barplot.svg", input.mid))) == FALSE){
      tmpdf <- subset(clonedf, grepl(input.mid, clonedf$samples) == TRUE)
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
  path.to.save.plot <- file.path(path.to.01.output, "chord_diagram", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  if (length(all.mids) > 1){
    count.share.clones <- data.frame(MID = all.mids)
    for (input.mid in all.mids){
      clones1 <- subset(clonesetsdf, clonesetsdf$id == input.mid)$CloneID %>% unique()
      count.share.clones[[input.mid]] <- unlist(lapply(
        count.share.clones$MID, function(x){
          clones <- subset(clonesetsdf, clonesetsdf$id == x)$CloneID
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
  path.to.save.plot <- file.path(path.to.01.output, "MHI_heatmap", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  if (length(all.mids) > 1){
    mhidf <- data.frame(MID = all.mids)
    for (input.mid1 in all.mids){
      mhidf[[input.mid1]] <- unlist(lapply(
        mhidf$MID, function(input.mid2){
          S <- subset(clonesetsdf, clonesetsdf$id %in% c(input.mid1, input.mid2))$CloneID %>% unique() %>% length()
          X <- subset(clonesetsdf, clonesetsdf$id == input.mid1) %>% nrow()
          Y <- subset(clonesetsdf, clonesetsdf$id == input.mid2) %>% nrow()
          
          x <- unlist(lapply(
            subset(clonesetsdf, clonesetsdf$id %in% c(input.mid1, input.mid2))$CloneID %>% unique(),
            function(x){
              return(subset(clonesetsdf, clonesetsdf$id == input.mid1 & clonesetsdf$CloneID == x) %>% nrow())
            }
          ))
          y <- unlist(lapply(
            subset(clonesetsdf, clonesetsdf$id %in% c(input.mid1, input.mid2))$CloneID %>% unique(),
            function(x){
              return(subset(clonesetsdf, clonesetsdf$id == input.mid2 & clonesetsdf$CloneID == x) %>% nrow())
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
  path.to.save.plot <- file.path(path.to.01.output, "alluvial_tracking", mouse.id)
  dir.create(path.to.save.plot, showWarnings = FALSE, recursive = TRUE)
  if (length(all.mids) > 1){
    for (input.mid in all.mids){
      if (file.exists(file.path(path.to.save.plot, sprintf("%s.tracking.svg", input.mid))) == FALSE){
        mid.clones <- subset(clonesetsdf, clonesetsdf$id == input.mid)
        mid.clones <- table(mid.clones$CloneID) %>% sort() %>% tail(10) %>% names()
        clone.in.mid <- data.frame(MID = all.mids)
        for (cloneid in mid.clones){
          clone.in.mid[[cloneid]] <- unlist(lapply(
            all.mids, function(x){
              (subset(clonesetsdf, clonesetsdf$id == x & clonesetsdf$CloneID == cloneid) %>% nrow())/nrow(clonesetsdf)
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

