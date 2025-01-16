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
                                                                                          "240829 sample sheet.xlsx"))
)

to.run.projects <- c("240805_BSimons_filterHT_cluster_renamed_240826_BSimons")

for (input.case in c("")){
  for (input.PROJECT in to.run.projects){
    print(sprintf("Working on project %s", input.PROJECT))
    # circos.group.type <- "VJcombi_CDR3_0.85"
    for (circos.group.type in c("VJaa", "VJnt", "VJcombi_CDR3_0.85")){
      print(sprintf("Working on circos group type: %s", circos.group.type))
      outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
      path.to.05.output <- file.path(outdir, "VDJ_output", 
                                     sprintf("05_output%s", input.case))
      path.to.08.output <- file.path(outdir, 
                                     "VDJ_output", 
                                     sprintf("08_output%s", input.case),
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
          
          # update 16.01.2025: change the order of M samples, P samples and 
          # MID of bulk samples. 
          mid.metadata <- bulk.metadata[[input.PROJECT]] %>%
            subset(mouse == mouse.id)
          
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
          mhidf <- data.frame(MID = factor(all.plot.samples, levels = all.plot.samples))
          
          #####
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
          write.csv(mhidf, file.path(path.to.08.output, sprintf("%s.MHI.csv", filename)))
        }else {
          print("File exists, reading in....")
          mhidf <- read.csv(file.path(path.to.08.output, sprintf("%s.MHI.csv", filename)))
        }
        library(viridis)
        if ("X" %in% colnames(mhidf)){
          mhidf <- subset(mhidf, select = -c(X))
        }
        mhidf.pivot <- mhidf %>% pivot_longer(!MID, names_to = "SampleID", values_to = "MHI") 
        mhidf.pivot$MID <- factor(mhidf.pivot$MID, levels = all.plot.samples)
        mhidf.pivot$SampleID <- factor(mhidf.pivot$SampleID, levels = all.plot.samples)
        print(paste(all.plot.samples, collapse = ", "))
        mhi.plot <- mhidf.pivot %>% 
          rowwise() %>%
          mutate(MHI.round = round(MHI, 3)) %>%
          ggplot(aes(x = MID, y = SampleID, fill = MHI)) + 
          geom_tile(color = "white") + 
          theme(axis.text.x = element_text(angle = 90)) + 
          scale_fill_gradient(high = "red", low = "gray") + 
          geom_text(aes(label = MHI.round), color = "white", size = 4) 
        ggsave(plot = mhi.plot, filename = sprintf("%s.MHI.svg", filename),
               path = path.to.08.output, device = "svg", dpi = 300, width = 14, height = 10)
      }
    }
  }
}





