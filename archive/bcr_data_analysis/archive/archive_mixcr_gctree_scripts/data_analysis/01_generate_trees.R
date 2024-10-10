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


outdir <- "/home/hieu/outdir"
PROJECT <- "mixcr_pipeline_output"

path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

path.to.mouse.output <- file.path(path.to.main.input, "mouse_based_output")
path.to.mid.output <- file.path(path.to.main.input, "mid_based_output")
##### METADATA
path.to.mid.metadata <- file.path(outdir, "FT_output", "mid_labels.csv")
mid.metadata <- read.csv(path.to.mid.metadata, sep = ";")

convert.timepoint <- to_vec(for (item in seq(1, length(unique(mid.metadata$population)))) sprintf("T%s", item))
names(convert.timepoint) <- unique(mid.metadata$population)
mid.metadata$timepoint <- unlist(lapply(mid.metadata$population, function(x){
  return(convert.timepoint[[x]])
}))

mid.metadata$label1 <- unlist(lapply(mid.metadata$population, function(x){
  return(str_replace(str_replace(x, "Ly6c[+]", ""), "Ly6c[-]", ""))
}))

##### Get tree information dataframes
all.tree.tsv.files <- Sys.glob(file.path(path.to.mouse.output, "*", "*trees.tsv"))
names(all.tree.tsv.files) <- unlist(lapply(
  all.tree.tsv.files, function(x){
    str_replace(str_replace(basename(x), "MID_list_", ""), "_trees.tsv", "")
  }
))

##### Get clone information dataframes
all.clones.tables <- Sys.glob(file.path(path.to.mid.output, "*", "*.reassigned_IGH.tsv"))
names(all.clones.tables) <- unlist(lapply(
  basename(all.clones.tables), function(x){
    return(str_replace(x, ".reassigned_IGH.tsv", ""))
  }
))

##### perform analysis for each mouse, all MID
# mouse.id <- "m11"

excluded.trees <- list(
  m32 = c(452),
  m53 = c(129)
)

my_color_palette <- hue_pal()(5)
names(my_color_palette) <- unique(mid.metadata$timepoint)

my_color_palette1 <- hue_pal()(3)
names(my_color_palette1) <- unique(mid.metadata$label1)

mid_color_pal <- hue_pal()(length(unique(mid.metadata$X)))
names(mid_color_pal) <- unique(mid.metadata$X)

mouse_color_pal <- hue_pal()(length(unique(mid.metadata$mouse)))
names(mouse_color_pal) <- unique(mid.metadata$mouse)

write.csv(data.frame(mid_color_pal), file.path(path.to.01.output, "mid_color_pal.csv"))
write.csv(data.frame(mouse_color_pal), file.path(path.to.01.output, "mouse_color_pal.csv"))

for (mouse.id in names(all.tree.tsv.files)){
# for (mouse.id in c("m11")){
  if (file.exists(file.path(path.to.01.output, sprintf("finished_saving_mouse_ID_%s.20240702.csv", mouse.id))) == FALSE){
    path.to.mouse <- file.path(path.to.mouse.output, sprintf("MID_list_%s", mouse.id))
    all.mids <- subset(mid.metadata, mid.metadata$mouse == mouse.id)$X
    
    ##### Get clone set files
    print("reading all clone set files ...")
    cloneset_files <- all.clones.tables[all.mids]
    
    clonesets <- read_tsv(cloneset_files, id="fileName") %>% 
      mutate(id=str_remove_all(fileName,".*\\/|.reassigned_IGH.tsv"),
             isotype=str_remove(allCHitsWithScore, "\\*.*"))
    
    print("generating statistics ...")
    stats_clonesets<-clonesets %>%
      group_by(id) %>%
      summarise(totalCloneInSample=n(),
                totalUMIsInSample=sum(uniqueMoleculeCount))
    
    ##### read tree information dataframe
    lng_db <- read_tsv(file.path(all.tree.tsv.files[[mouse.id]])) %>% 
      mutate(id=str_remove_all(fileName,".*\\/|.reassigned.clns"),
             isotype=str_remove(bestCHit,"\\*00")) %>% 
      left_join(stats_clonesets,by ="id")
    
    print("filtering tree data ...")
    lng_db <- subset(lng_db, is.na(lng_db$id) == FALSE)
    lng_db <- merge(lng_db, mid.metadata, by.x = "id", by.y = "X")
    lng_db <- merge(lng_db, subset(clonesets, select = c(cloneId, uniqueMoleculeCount)), by.x = "cloneId", by.y = "cloneId")
    
    ##### Pivot to wide format
    print("pivot tree dataframe into wide format ...")
    lngs_trees_wide<-lng_db %>%
      filter(!fileName=="") %>% 
      group_by(treeId,timepoint) %>% 
      summarise(nClones=n(),
                nClonesNormalized=n()/dplyr::first(totalCloneInSample),
                nUMIsNormalized=sum(uniqueMoleculeCount)/dplyr::first(totalUMIsInSample),
                nUMIs=sum(uniqueMoleculeCount),
                isotypes=list(unique(isotype)),
                totalCloneInSample=dplyr::first(totalCloneInSample),
                totalUMIsInSample=dplyr::first(totalUMIsInSample)
      ) %>% 
      ungroup() %>%
      group_by(treeId) %>% 
      mutate(timepoint=factor(timepoint, levels= unique(mid.metadata$timepoint)),
             timepoints=timepoint %>%unique() %>% sort() %>%  paste(collapse=","),
             isotypes=c(isotypes) %>% unlist() %>% unique() %>%sort() %>%  paste(collapse=","),
             totalNClones=sum(nClones),
             totalCount=sum(nUMIs)) %>% 
      ungroup() %>% 
      arrange(timepoint) %>% 
      pivot_wider(id_cols = c(treeId,isotypes,timepoints,totalNClones,totalCount), 
                  names_from = c(timepoint),
                  values_from = c(nClones,nUMIs,nClonesNormalized,nUMIsNormalized),
                  values_fill = 0.00) %>% 
      mutate(across(starts_with("nClonesN"),.fns = ~ ifelse(.x==0, min(.x[.x>0])/2, .x )),
             across(starts_with("nUmisN"),.fns = ~ ifelse(.x==0, min(.x[.x>0])/2, .x ))) %>% 
      mutate(switchedIsotypePresent=ifelse(str_detect(isotypes, "IGHA|IGHG|IGHE"), TRUE, FALSE))
    
    ##### Generate BCR MIXCR trees
    dir.create(file.path(path.to.01.output, "color_by_AASeqCDR3", mouse.id), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(path.to.01.output, "color_by_ID", mouse.id), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(path.to.01.output, "color_by_YFP", mouse.id), showWarnings = FALSE, recursive = TRUE)
    dir.create(file.path(path.to.01.output, "color_by_YFP_noLabel", mouse.id), showWarnings = FALSE, recursive = TRUE)
    all.trees <- setdiff(unique(lng_db$treeId), excluded.trees[[mouse.id]])
    save.plot.list <- list()
    for (treeID in all.trees){
      print(sprintf("plotting tree iD %s", treeID))
      path.to.tree.file <- file.path(path.to.mouse, "IM_newick", sprintf("%s.tree", treeID))
      tree <- read.tree(path.to.tree.file)
      tree_meta <- lng_db %>%
        filter(treeId==treeID) %>%
        select(nodeId,
               isotype,
               timepoint, 
               label1,
               aaSeqCDR3,
               uniqueMoleculeCount,
               cloneId) %>% 
        dplyr::rename(ID=nodeId)
      
      ## draw initial tree 
      tryCatch(
        {
          if (file.exists(file.path(path.to.01.output, "color_by_ID", mouse.id, sprintf("tree_%s.svg", treeID))) == FALSE){
            g.timepoint <- ggtree(tree) %<+% tree_meta +
              geom_text(aes(label=timepoint), hjust=-.5,size=4)+
              geom_tippoint(aes(color = timepoint , size=uniqueMoleculeCount))+
              theme(legend.position = 'bottom')+
              scale_size_continuous(range = c(2, 7)) + 
              ggplot2:::manual_scale(
                'color', 
                values = setNames(my_color_palette, names(my_color_palette)))
            ggsave(plot = g.timepoint, filename = sprintf("tree_%s.svg", treeID), path = file.path(path.to.01.output, "color_by_ID", mouse.id), device = "svg", width = 14, height = 10, dpi = 300)
          }
          
          if (file.exists(file.path(path.to.01.output, "color_by_YFP", mouse.id, sprintf("tree_%s.svg", treeID))) == FALSE){
            g.label1 <- ggtree(tree) %<+% tree_meta +
              geom_text(aes(label=label1), hjust=-.5,size=4)+
              geom_tippoint(aes(color = label1 , size=uniqueMoleculeCount))+
              theme(legend.position = 'bottom')+
              scale_size_continuous(range = c(2, 7)) + 
              ggplot2:::manual_scale(
                'color', 
                values = setNames(my_color_palette1, names(my_color_palette1)))
            ggsave(plot = g.label1, filename = sprintf("tree_%s.svg", treeID), path = file.path(path.to.01.output, "color_by_YFP", mouse.id), device = "svg", width = 14, height = 10, dpi = 300)
          }
          
          if (file.exists(file.path(path.to.01.output, "color_by_YFP_noLabel", mouse.id, sprintf("tree_%s.svg", treeID))) == FALSE){
            g.label1 <- ggtree(tree) %<+% tree_meta +
              geom_tippoint(aes(color = label1 , size=uniqueMoleculeCount))+
              theme(legend.position = 'bottom')+
              scale_size_continuous(range = c(2, 7)) + 
              ggplot2:::manual_scale(
                'color', 
                values = setNames(my_color_palette1, names(my_color_palette1)))
            ggsave(plot = g.label1, filename = sprintf("tree_%s.svg", treeID), path = file.path(path.to.01.output, "color_by_YFP_noLabel", mouse.id), device = "svg", width = 14, height = 10, dpi = 300)
          }
          
          if (file.exists(file.path(path.to.01.output, "color_by_AASeqCDR3", mouse.id, sprintf("tree_%s.svg", treeID))) == FALSE){
            g.AAseq <- ggtree(tree) %<+% tree_meta +
              geom_text(aes(label=timepoint), hjust=-.5,size=4)+
              geom_tippoint(aes(color = timepoint , size=uniqueMoleculeCount))+
              theme(legend.position = 'bottom')+
              scale_size_continuous(range = c(2, 7))
            ggsave(plot = g.AAseq, filename = sprintf("tree_%s.svg", treeID), path = file.path(path.to.01.output, "color_by_AASeqCDR3", mouse.id), device = "svg", width = 14, height = 10, dpi = 300)
          }
        },
        error = function(cond) {
        }
      )
    }
    write.csv(data.frame(status = c(sprintf("finished_ID_%s", mouse.id))), file.path(path.to.01.output, sprintf("finished_saving_mouse_ID_%s.20240702.csv", mouse.id)))
  }
}
