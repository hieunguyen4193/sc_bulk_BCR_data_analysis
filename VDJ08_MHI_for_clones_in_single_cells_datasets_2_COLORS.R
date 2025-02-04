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
library(viridis)
#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

# for (input.case in c(
#   "_filterHT_cutoff_40",
#   "_filterHT_cutoff_40_selected_cells_group1",
#   "_filterHT_cutoff_40_selected_cells_group2",
#   "_filterHT_cutoff_40_selected_cells_group3",
#   "_filterHT_remove_last_ht",
#   "_filterHT_remove_last_ht_selected_cells_group1",
#   "_filterHT_remove_last_ht_selected_cells_group2",
#   "_filterHT_remove_last_ht_selected_cells_group3")){

# to.run.projects <- c("240805_BSimons_240826_BSimons",
#                      "241002_BSimons_241104_BSimons_241031_BSimons")

bulk.metadata <- list(
  `240826_BSimons` = readxl::read_excel(file.path(path.to.main.src, 
                                                  "preprocessing", 
                                                  "240826_BSimons", 
                                                  "240829 sample sheet.xlsx")),
  `240805_BSimons_filterHT_cluster_renamed_240826_BSimons` = readxl::read_excel(file.path(path.to.main.src, 
                                                                                          "preprocessing", 
                                                                                          "240826_BSimons", 
                                                                                          "240829 sample sheet.xlsx")),
  `241002_241104_BSimons_241031_BSimons` = readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241031_BSimons/241031_sample_sheet.xlsx")
)

tmp.metadata <- bulk.metadata[["241002_241104_BSimons_241031_BSimons"]]
colnames(tmp.metadata) <- c("MID", "mouse", "organ", "YFP", "sample")
tmp.metadata <- tmp.metadata %>% rowwise() %>%
  mutate(MID = sprintf("MID%s", MID))
bulk.metadata[["241002_241104_BSimons_241031_BSimons"]] <- tmp.metadata

to.run.projects <- c("241002_241104_BSimons_241031_BSimons",
                     "240805_BSimons_filterHT_cluster_renamed_240826_BSimons")

input.case <- ""

for (input.PROJECT in to.run.projects){
  print(sprintf("Working on project %s", input.PROJECT))
  circos.group.type <- "VJcombi_CDR3_0.85"
  print(sprintf("Working on circos group type: %s", circos.group.type))
  outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
  path.to.05.output <- file.path(outdir, "VDJ_output", 
                                 sprintf("05_output%s", input.case))
  path.to.08.output <- file.path(outdir, 
                                 "VDJ_output", 
                                 sprintf("08_output_2COLOR_%s", input.case),
                                 input.PROJECT,
                                 circos.group.type,
                                 "filter_clone_FALSE",
                                 "MHI_results")
  dir.create(path.to.08.output, showWarnings = FALSE, recursive = TRUE)
  path.to.input <- file.path(path.to.05.output, 
                             input.PROJECT, 
                             circos.group.type, 
                             "circos_plot", 
                             "filter_clone_FALSE")
  all.files <- Sys.glob(file.path(path.to.input, "*.csv"))
  
  for (input.file in all.files){
    mouse.id <- str_split(basename(input.file), "_")[[1]][[1]]
    print(sprintf("Working on sample %s", input.file))
    tmpdf <- read.csv(input.file) %>% subset(select = -c(X))
    colnames(tmpdf) <- c("CloneID", colnames(tmpdf)[2:length(colnames(tmpdf))])
    filename <- basename(input.file) %>% str_replace(".csv", "")
    if (file.exists(file.path(path.to.08.output, sprintf("%s.MHI.csv", filename))) == FALSE){
      mhidf <- data.frame(MID = to_vec(
        for (item in colnames(tmpdf)){
          if (grepl("_Freq", item) == TRUE){
            str_replace(item, "_Freq", "")
          }
        }
      ))
      
      #####----------------------------------------------------------------#####
      # update 16.01.2025: change the order of M samples, P samples and 
      # MID of bulk samples. 
      #####----------------------------------------------------------------#####
      mid.metadata <- bulk.metadata[[input.PROJECT]] %>%
        subset(mouse == mouse.id)
      if (input.PROJECT == "240805_BSimons_filterHT_cluster_renamed_240826_BSimons"){
        MID.samples <- to_vec(
          for (item in c("SI prox", "SI mid", "SI dist", "colon")){
            subset(mid.metadata, mid.metadata$sample == sprintf("%s %s", mouse.id, item))$MID
          }
        )
        all.plot.samples <- mhidf$MID %>% unique()
        p.sample <- all.plot.samples[grepl("P", all.plot.samples)]
        p.sample <- sort(p.sample, decreasing = TRUE)
        
        m.sample <- all.plot.samples[grepl("M", all.plot.samples) == TRUE & 
                                       grepl("MID", all.plot.samples) == FALSE]
        m.sample <- sort(m.sample, decreasing = FALSE)
        
        all.plot.samples <- c(p.sample, m.sample, MID.samples)
      } else if (input.PROJECT == "241002_241104_BSimons_241031_BSimons"){
        if (grepl("MERGE_YFP", basename(input.file)) == FALSE){
          MID.samples <- mid.metadata$MID
          all.plot.samples <- mhidf$MID %>% unique()
          if (mouse.id == "m3"){
            if (grepl("hashtags", basename(input.file)) == TRUE){
              p.sample <- c("PP3_HT1", "PP3_HT2", "PP3_HT3")
            } else{
              p.sample <- c("PP3")
            }
          } else if (mouse.id == "m7"){
            if (grepl("hashtags", basename(input.file)) == TRUE){
              p.sample <- c("PP7_HT3", "PP7_HT1", "PP7_HT2")
            } else{
              p.sample <- c("PP7")
            }
          }
          all.plot.samples <- c(p.sample, MID.samples)
        } else {
          all.plot.samples <- mhidf$MID
        }
      }
      
      mhidf <- data.frame(MID = factor(all.plot.samples, levels = all.plot.samples))
      mhidf.count <- data.frame(MID = factor(all.plot.samples, levels = all.plot.samples))
      #####----------------------------------------------------------------#####
      ##### calculate MHI
      #####----------------------------------------------------------------#####
      for (input.mid1 in mhidf$MID){
        print(sprintf("working on %s", input.mid1))
        mhidf[[input.mid1]] <- unlist(
          lapply(mhidf$MID, function(input.mid2){
            X <- tmpdf[[sprintf("%s_Count", input.mid1)]] %>% sum()
            Y <- tmpdf[[sprintf("%s_Count", input.mid2)]] %>% sum()
            
            x <- unlist(lapply(
              tmpdf$CloneID %>% unique(),
              function(x){
                return(subset(tmpdf, tmpdf$CloneID == x)[[sprintf("%s_Count", input.mid1)]])
              }
            ))
            y <- unlist(lapply(
              tmpdf$CloneID %>% unique(),
              function(x){
                return(subset(tmpdf, tmpdf$CloneID == x)[[sprintf("%s_Count", input.mid2)]])
              }
            ))
            
            nom <- 2 * sum(x * y)
            det <- (sum(x^2)/X^2) + (sum(y^2)/Y^2)
            mhi <- nom/(as.double(X)*as.double(Y) * det)
            return(mhi)
          })
        )
      }
      #####----------------------------------------------------------------#####
      ##### calculate num share clones
      #####----------------------------------------------------------------#####
      for (input.mid1 in mhidf.count$MID){
        print(sprintf("working on %s", input.mid1))
        mhidf.count[[input.mid1]] <- unlist(
          lapply(mhidf.count$MID, function(input.mid2){
            checkdf <- tmpdf[, c(sprintf("%s_Count", input.mid1), sprintf("%s_Count", input.mid2))]
            colnames(checkdf) <- c("Sample1", "Sample2")
            return(nrow(subset(checkdf, Sample1 != 0 & Sample2 != 0)))
          }
        )
      )
      }
      write.csv(mhidf.count, file.path(path.to.08.output, sprintf("%s.MHI_countShare.csv", filename)))
      write.csv(mhidf, file.path(path.to.08.output, sprintf("%s.MHI.csv", filename)))
    } else {
      print("File exists, reading in....")
      mhidf.count <- read.csv(file.path(path.to.08.output, sprintf("%s.MHI_countShare.csv", filename)))
      mhidf <- read.csv(file.path(path.to.08.output, sprintf("%s.MHI.csv", filename)))
    }

    #####----------------------------------------------------------------#####
    ##### plot MHI heatmap
    #####----------------------------------------------------------------#####
    if ("X" %in% colnames(mhidf)){
      mhidf <- subset(mhidf, select = -c(X))
      all.plot.samples <- mhidf$MID
    }
    if ("X" %in% colnames(mhidf.count)){
      mhidf.count <- subset(mhidf.count, select = -c(X))
      all.plot.samples <- mhidf.count$MID
    }
    mhidf.pivot <- mhidf %>% pivot_longer(!MID, names_to = "SampleID", values_to = "MHI") %>% 
      rowwise() %>%
      mutate(MHI.round = round(MHI, 3))
    mhidf.pivot$MID <- factor(mhidf.pivot$MID, levels = all.plot.samples)
    mhidf.pivot$SampleID <- factor(mhidf.pivot$SampleID, levels = all.plot.samples)
    print(paste(all.plot.samples, collapse = ", "))
    
    mhi.plot <- mhidf.pivot  %>%
      ggplot(aes(x = MID, y = SampleID, fill = MHI)) + 
      geom_tile(color = "white") + 
      theme(axis.text.x = element_text(angle = 90)) + 
      scale_fill_gradient(high = "red", low = "#f7fafd") + 
      geom_text(aes(label = MHI.round), color = "black", size = 4) 
    ggsave(plot = mhi.plot, filename = sprintf("%s.MHI.svg", filename),
           path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
    
    
    palette_red_blue <- colorRampPalette(colors = c("#f7fafd", "lightblue", "red"))
    
    mhidf.pivot$z_categorised <- cut(mhidf.pivot$MHI, c(seq(-0.01, 0.1, length.out = 100), 
                                                        seq(0.101, 1, length.out = 100)) )
    mhi.plot <- mhidf.pivot  %>%
      ggplot(aes(x = MID, y = SampleID, fill = z_categorised)) + 
      geom_tile(color = "white") + 
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "none") + 
      scale_fill_manual(values = palette_red_blue(length(unique(mhidf.pivot$z_categorised)))) + 
      geom_text(aes(label = MHI.round), color = "black", size = 4) 
    ggsave(plot = mhi.plot, filename = sprintf("%s.MHI_2COLOR.svg", filename),
           path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
    
    #####----------------------------------------------------------------#####
    ##### plot SHARE CLONE COUNT HEATMMAP
    #####----------------------------------------------------------------#####
    if ("X" %in% colnames(mhidf.count)){
      mhidf.count <- subset(mhidf.count, select = -c(X))
      all.plot.samples <- mhidf.count$MID
    }
    mhidf.count.pivot <- mhidf.count %>% pivot_longer(!MID, names_to = "SampleID", values_to = "MHI") %>% 
      rowwise() %>%
      mutate(MHI.round = round(MHI, 3))
    mhidf.count.pivot$MID <- factor(mhidf.count.pivot$MID, levels = all.plot.samples)
    mhidf.count.pivot$SampleID <- factor(mhidf.count.pivot$SampleID, levels = all.plot.samples)
    print(paste(all.plot.samples, collapse = ", "))
    mhi.plot <- mhidf.count.pivot  %>%
      ggplot(aes(x = MID, y = SampleID, fill = MHI)) + 
      geom_tile(color = "white") + 
      theme(axis.text.x = element_text(angle = 90)) + 
      scale_fill_gradient(high = "red", low = "#f7fafd") + 
      geom_text(aes(label = MHI.round), color = "black", size = 4) 
    ggsave(plot = mhi.plot, filename = sprintf("%s.MHI_countShare.svg", filename),
           path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
    
    palette_red_blue <- colorRampPalette(colors = c("#f7fafd", "lightblue", "red"))
    mhidf.count.pivot$z_categorised <- cut(mhidf.count.pivot$MHI, c(seq(-0.01, 20, length.out = 100), 
                                                        seq(20.01, max(mhidf.count.pivot$MHI), 
                                                            length.out = max(mhidf.count.pivot$MHI) - 20)) )
    mhi.plot <- mhidf.count.pivot  %>%
      ggplot(aes(x = MID, y = SampleID, fill = z_categorised)) + 
      geom_tile(color = "white") + 
      theme(axis.text.x = element_text(angle = 90),
            legend.position = "none") + 
      scale_fill_manual(values = palette_red_blue(length(unique(mhidf.count.pivot$z_categorised)))) + 
      geom_text(aes(label = MHI.round), color = "black", size = 4) 
    ggsave(plot = mhi.plot, filename = sprintf("%s.MHI_2COLOR_countShare.svg", filename),
           path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
  }
}


