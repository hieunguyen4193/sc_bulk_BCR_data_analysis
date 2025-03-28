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
library(circlize)

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

circos.group.type <- "VJaa"

source(file.path(path.to.main.src, "convert_sampleID_to_locationName.R"))

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
bulk.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240826_BSimons/240829 sample sheet.xlsx")
sc.projects <- c("240805_BSimons_filterHT_cluster_renamed")
bulk.projects <- c("240826_BSimons")
list.of.PROJECT <- c(sc.projects, bulk.projects)

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
path.to.15.output <- file.path(outdir, "VDJ_output", "15_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.15.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
#### read the Seurat object of single cell dataset
######----------------------------------------------------------------------#####
s.obj <- list()
sc.meta.data <- list()
list.of.ht <- list()
for (PROJECT in sc.projects){
  s.obj[[PROJECT]] <- readRDS(path.to.all.s.obj[[PROJECT]])
  sc.meta.data[[PROJECT]] <- s.obj[[PROJECT]]@meta.data %>% rownames_to_column("barcode")
  list.of.ht[[PROJECT]] <- list()
  for (sample.id in unique(sc.meta.data[[PROJECT]]$name)){
    list.of.ht[[PROJECT]][[sample.id]] <- unique(subset(sc.meta.data[[PROJECT]], sc.meta.data[[PROJECT]]$name == sample.id)$HTO_classification)
  }
}

#####----------------------------------------------------------------------#####
##### metadata  preprocessing
#####----------------------------------------------------------------------#####
excluded.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
meta.data <- readxl::read_excel(file.path(path.to.05.output, "all_data_metadata.xlsx"))
meta.data.splitted.or.not <- list(
  with_hashtags = subset(meta.data, meta.data$SampleID %in% excluded.samples == FALSE),
  without_hashtags = subset(meta.data, grepl("_", meta.data$SampleID) == FALSE)
)

tmp.metadata <- meta.data.splitted.or.not$with_hashtags

mouse.sample.list <- list()
for (mouse.id in unique(tmp.metadata$mouse)){
  mouse.sample.list[[mouse.id]] <- subset(tmp.metadata, tmp.metadata$mouse == mouse.id)$SampleID %>% sort()
}

#####----------------------------------------------------------------------#####
##### PREPROCESS FILE: ASSIGN CLONES TO CLUSTERS 0.85 SIMILARITY
#####----------------------------------------------------------------------#####
all.input.files <- Sys.glob(file.path(path.to.05.output,
                                      circos.group.type,
                                      "*.simplified.csv"))

input.metadata <- data.frame(
  path = all.input.files,
  SampleID = to_vec(for (item in all.input.files){
    str_replace(basename(item), ".simplified.csv", "") 
  }))

all.input.files <- input.metadata$path
names(all.input.files) <- input.metadata$SampleID

##### MAIN RUN
all.samples <- names(all.input.files)
keep.samples <- setdiff(all.samples, excluded.samples)

input.files <- all.input.files[keep.samples] # run for all samples, not for each mouse only.

thres.dis <- 0.15
thres <- 0.85

if (file.exists(file.path(path.to.15.output, "new.clonesets.csv")) == FALSE){
  clonesets <- read_tsv(input.files, id = "fileName") %>%
    rowwise() %>%
    mutate(fileName = basename(fileName) %>% str_replace(".simplified.csv", "")) %>%
    mutate(bestVHit = str_split(bestVHit, "[*]")[[1]][[1]]) %>%
    mutate(bestJHit = str_split(bestJHit, "[*]")[[1]][[1]]) %>%
    mutate(len = nchar(nSeqCDR3)) %>%
    mutate(VJ.len.combi = sprintf("%s_%s_%s", bestVHit, bestJHit, len ))
  
  colnames(clonesets) <- c("fileName", "id", "cloneCount", "bestVHit", "bestJHit", "seq", "len", "VJ.len.combi")
  
  ##### Group sequences + Gene usages to clones
  new.clonesets <- data.frame()
  for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
    tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
    seqs <- unique(tmpdf$seq)
    print(sprintf("VJ.len.combi: %s, num seqs: %s", input.VJ.combi, length(seqs)))
    if (length(seqs) >= 2){
      cluster.output <- assign_clusters_to_sequences(seqs = seqs, threshold = thres.dis)$res
      tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
        tmpdf$seq, function(x){
          return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
        }
      ))    
    } else {
      tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
    }
    new.clonesets <- rbind(new.clonesets, tmpdf)
  }
  write.csv(new.clonesets, file.path(path.to.15.output, "new.clonesets.csv"))
  
} else {
  print("reading in existing clonesets data")
  new.clonesets <- read.csv(file.path(path.to.15.output, "new.clonesets.csv"))
}

#####----------------------------------------------------------------------#####
##### summary of all clones and their counts
#####----------------------------------------------------------------------#####
all.plot.samples <- c(mouse.sample.list$m1, mouse.sample.list$m2, mouse.sample.list$m3)

if (file.exists(file.path(path.to.15.output, "clonedf_allSamples_allMice.csv")) == FALSE){
  clonedf <- data.frame(CloneID = unique(new.clonesets$VJcombi_CDR3_0.85))
  for (sample.id in all.plot.samples){
    print(sprintf("wokring on sample %s", sample.id))
    clonedf[[sprintf("%s_Count", sample.id)]] <- unlist(lapply(
      clonedf$CloneID, function(x){
        tmpdf <- subset(new.clonesets, new.clonesets$VJcombi_CDR3_0.85 == x & new.clonesets$fileName == sample.id)
        if (nrow(tmpdf) == 0){
          output <- 0
        } else {
          output <- tmpdf$cloneCount %>% sum()
        }
        if (length(output) != 1){
          print(x)
        }
        return(output)
      }
    ))
  }
  write.csv(clonedf, file.path(path.to.15.output, "clonedf_allSamples_allMice.csv"))
} else {
  print("reading in existing data...")
  clonedf <- read.csv(file.path(path.to.15.output, "clonedf_allSamples_allMice.csv"))
}

#####----------------------------------------------------------------------#####
##### generate MHI values for all pair of samples in all mice
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.15.output, "mhidf.csv")) == FALSE){
  mhidf <- data.frame(MID = factor(all.plot.samples, levels = all.plot.samples))
  
  for (input.mid1 in mhidf$MID){
    print(sprintf("working on %s", input.mid1))
    mhidf[[input.mid1]] <- unlist(
      lapply(mhidf$MID, function(input.mid2){
        X <- clonedf[[sprintf("%s_Count", input.mid1)]] %>% sum()
        Y <- clonedf[[sprintf("%s_Count", input.mid2)]] %>% sum()
        
        x <- unlist(lapply(
          clonedf$CloneID %>% unique(),
          function(x){
            return(subset(clonedf, clonedf$CloneID == x)[[sprintf("%s_Count", input.mid1)]])
          }
        ))
        y <- unlist(lapply(
          clonedf$CloneID %>% unique(),
          function(x){
            return(subset(clonedf, clonedf$CloneID == x)[[sprintf("%s_Count", input.mid2)]])
          }
        ))
        
        nom <- 2 * sum(x * y)
        det <- (sum(x^2)/X^2) + (sum(y^2)/Y^2)
        mhi <- nom/(as.double(X)*as.double(Y) * det)
        return(mhi)
      })
    )
  }
  write.csv(mhidf, file.path(path.to.15.output, "mhidf.csv"))  
} else {
  mhidf <- read.csv(file.path(path.to.15.output, "mhidf.csv")) %>%
    subset(select = -c(X))
}

colnames(mhidf) <- to_vec(for (item in colnames(mhidf)){
  if (item == "MID"){
    "SampleID"
  } else {
    item
  }
})

#####----------------------------------------------------------------------#####
##### HEAT MAP FOR MHI
#####----------------------------------------------------------------------#####
mhidf$SampleID <- factor(mhidf$SampleID, 
                         levels = all.plot.samples)
pivot.mhidf <- mhidf %>% pivot_longer(!SampleID, names_to = "SampleID2", values_to = "mhi")
pivot.mhidf$SampleID <- factor(pivot.mhidf$SampleID, 
                               levels = all.plot.samples)
pivot.mhidf$SampleID2 <- factor(pivot.mhidf$SampleID2, 
                                levels = all.plot.samples)

heatmap.mhi <- pivot.mhidf %>% ggplot(aes(x = SampleID, y = SampleID2, fill = mhi)) + 
  geom_tile(color = "white") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_fill_gradient(high = "red", low = "#f7fafd") 
#####----------------------------------------------------------------------#####
##### Generate count dataframe, number of share clones between all samples in all mice
#####----------------------------------------------------------------------#####
countdf <- data.frame(SampleID = unique(new.clonesets$fileName) %>% sort())

for (sample.id in countdf$SampleID){
  countdf[[sample.id]] <- unlist(lapply(
    countdf$SampleID, function(x){
      x.clones <- subset(new.clonesets, new.clonesets$fileName == x)$VJcombi_CDR3_0.85
      sample.id.clones <- subset(new.clonesets, new.clonesets$fileName == sample.id)$VJcombi_CDR3_0.85
      return(length(intersect(x.clones, sample.id.clones)))
    }
  ))
}

countdf$SampleID <- factor(countdf$SampleID, 
                           levels = all.plot.samples)
pivot.countdf <- countdf %>% pivot_longer(!SampleID, names_to = "SampleID2", values_to = "count")
pivot.countdf$SampleID <- factor(pivot.countdf$SampleID, 
                                 levels = all.plot.samples)
pivot.countdf$SampleID2 <- factor(pivot.countdf$SampleID2, 
                                 levels = all.plot.samples)

heatmap.count <- pivot.countdf %>% ggplot(aes(x = SampleID, y = SampleID2, fill = count)) + 
  geom_tile(color = "white") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_fill_gradient(high = "red", low = "#f7fafd") +
  geom_text(aes(label = count), color = "black", size = 4) 
