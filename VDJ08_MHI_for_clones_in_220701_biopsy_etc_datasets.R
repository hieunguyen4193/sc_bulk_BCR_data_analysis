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
library(circlize)
#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

mid.metadata <- read.csv("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/220701_etc_biopsies/metadata.csv", sep =";")
yfp.mids <- list()
sample.list <- list()
for (mouse.id in unique(mid.metadata$mouse)){
  sample.list[[mouse.id]] <- subset(mid.metadata, mid.metadata$mouse == mouse.id)$X
  yfp.mids[[mouse.id]] <- list(
    all_w_biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id)$X,
    all_yfp = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP", mid.metadata$population) == TRUE)$X,
    pos = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[+]", mid.metadata$population) == TRUE)$X,
    neg = subset(mid.metadata, mid.metadata$mouse == mouse.id & grepl("YFP[-]", mid.metadata$population) == TRUE)$X,
    biopsy = subset(mid.metadata, mid.metadata$mouse == mouse.id & mid.metadata$population == "biopsy")$X)   
}

PROJECT <- "220701_etc_biopsies"
path.to.VDJ.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> NEW DATA
#####----------------------------------------------------------------------#####
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output_NEW_20250205", PROJECT)
path.to.08.output <- file.path(outdir, "VDJ_output", "08_output", sprintf("%s_new_MHI_20250212_biopsy", PROJECT))
dir.create(path.to.08.output, showWarnings = FALSE, recursive = TRUE)

mid.metadata <- mid.metadata %>% rowwise() %>%
  mutate(sample_type = str_replace(str_replace(population, "Ly6c[+]", ""), "Ly6c[-]", "")) %>%
  mutate(sample_type = str_replace(sample_type, "YFP[+]", "YFP_positive")) %>%
  mutate(sample_type = str_replace(sample_type, "YFP[-]", "YFP_negative"))

count.mice <- table(mid.metadata$mouse)
plot.mice <- count.mice[count.mice >= 2] %>% names()

for (mouse.id in plot.mice){
  all.input.files <- Sys.glob(file.path(path.to.05.output, "circos_plot", "filter_clone_FALSE", 
                                        sprintf("%s*.csv", mouse.id) ))
  path.to.file <- all.input.files[[1]]
  filename <- basename(path.to.file) %>% str_replace(".csv", "")
  tmpdf <- read.csv(path.to.file) %>%
    subset(select = -c(X))
  
  all.plot.samples <- c("biopsy", "YFP_negative", "YFP_positive")
  
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
          tmpdf$id %>% unique(),
          function(x){
            return(subset(tmpdf, tmpdf$id == x)[[sprintf("%s_Count", input.mid1)]])
          }
        ))
        y <- unlist(lapply(
          tmpdf$id %>% unique(),
          function(x){
            return(subset(tmpdf, tmpdf$id == x)[[sprintf("%s_Count", input.mid2)]])
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
}
