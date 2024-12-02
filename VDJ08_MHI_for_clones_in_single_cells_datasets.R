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

input.PROJECT <- "240805_BSimons_240826_BSimons"
circos.group.type <- "VJaa"

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output")
path.to.08.output <- file.path(outdir, 
                               "VDJ_output", 
                               "08_output",
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
    
    for (input.mid1 in mhidf$MID){
      print(sprintf("working on %s", input.mid1))
      mhidf[[input.mid1]] <- unlist(
        lapply(mhidf$MID, function(input.mid2){
          X <- tmpdf[tmpdf[[sprintf("%s_Freq", input.mid1)]] != 0,] %>% nrow()
          Y <- tmpdf[tmpdf[[sprintf("%s_Freq", input.mid2)]] != 0,] %>% nrow()
          
          x <- unlist(lapply(
            tmpdf$CloneID %>% unique(),
            function(x){
              return(subset(tmpdf, 
                            tmpdf[[sprintf("%s_Count", input.mid1)]]!= 0 & 
                              tmpdf$CloneID == x) %>% nrow())
            }
          ))
          y <- unlist(lapply(
            tmpdf$CloneID %>% unique(),
            function(x){
              return(subset(tmpdf, 
                            tmpdf[[sprintf("%s_Count", input.mid2)]]!= 0 &
                              tmpdf$CloneID == x) %>% nrow())
            }
          ))
          
          nom <- 2 * sum(x * y)
          det <- (sum(x^2)/X^2) + (sum(y^2)/Y^2)
          mhi <- nom/(X*Y * det)
          return(mhi)
        })
      )
    }
    write.csv(mhidf, file.path(path.to.08.output, sprintf("%s.MHI.csv", filename)))
  }else {
    mhidf <- read.csv(file.path(path.to.08.output, sprintf("%s.MHI.csv", filename)))
  }
  
}



