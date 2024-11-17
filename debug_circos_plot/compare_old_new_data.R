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
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

# outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
outdir <- "/media/hieunguyen/HNHD01/outdir/sc_bulk_BCR_data_analysis_v0.1"
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
path.to.VDJ.output <- file.path( outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> NEW DATA
#####----------------------------------------------------------------------#####
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", PROJECT)
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, "circos_plot"), showWarnings = FALSE, recursive = TRUE)

path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

clone.obj <- run_preprocessing_all_bulk_VDJ_data(path.to.mid.output = path.to.mid.output,
                                                 path.to.save.output = path.to.save.output,
                                                 PROJECT = PROJECT,
                                                 thres = thres, 
                                                 thres.dis = thres.dis,
                                                 savefile = savefile,
                                                 verbose = verbose,
                                                 rerun = rerun,
                                                 define.clone.clusters = define.clone.clusters)   


full.clonedf <- clone.obj$clonesets
full.clonedf <- subset(full.clonedf, full.clonedf$aaSeqCDR3 != "region_not_covered")

new.mids <- list()
for (mid in unique(full.clonedf$id)){
 tmpdf <- subset(full.clonedf, full.clonedf$id == mid) 
 tmpdf <- subset(tmpdf, 
                 select = c(
                   nSeqCDR3,
                   aaSeqCDR3,
                   V.gene,
                   J.gene,
                   VJ.len.combi
                 ))
 tmpdf <- tmpdf %>% rowwise() %>%
   mutate(VJnt = sprintf("%s_%s_%s", 
                         str_split(V.gene, "[*]")[[1]][[1]], 
                         str_split(J.gene, "[*]")[[1]][[1]], 
                         nSeqCDR3))
 tmp.countdf <- data.frame(table(tmpdf$VJnt))
 colnames(tmp.countdf) <- c("id", "cloneCount")
 tmp.countdf <- tmp.countdf %>% rowwise() %>%
   mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
   mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
   mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]]) %>%
   arrange(desc(cloneCount))
 new.mids[[mid]] <- tmp.countdf
}


#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> OLD DATA
#####----------------------------------------------------------------------#####

old.mids <- list()
path.to.old.data <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets/test_circos_inputs"
all.old.files <- Sys.glob(file.path(path.to.old.data, "mid*_clones.csv"))
for (i in all.old.files){
  mid <- basename(i)
  mid <- str_split(mid, "_")[[1]][[1]] %>% toupper()
  if (mid %in% c("MID4", "MID55", "MID5")){
    sep <- ";"
  } else {
    sep = "\t"
  }
  tmp.olddf <- read.csv(i, sep = sep) %>%
    subset(select = c(nSeqCDR3, bestVHit, bestJHit, cloneCount)) %>%
    rowwise() %>%
    mutate(id = sprintf("%s_%s_%s", bestVHit, bestJHit, nSeqCDR3))
  tmp.olddf <- tmp.olddf[, c(colnames(new.mids$MID10))]
 
  old.mids[[mid]] <- tmp.olddf
}

data <- list(
  old = old.mids,
  new = new.mids
  )

print(sprintf("Number of samples in old MID: %s", length(names(old.mids))))
print(sprintf("Number of samples in new MID: %s", length(names(new.mids))))
print(sprintf("Diff samples: %s", paste(setdiff(names(new.mids), names(old.mids)), collapse = ", ")))

checkdf <- data.frame(MID = intersect(names(old.mids), names(new.mids)))

for (input.case in c("old", "new")){
  input.data <- data[[input.case]]
  ##### num clone
  checkdf[[sprintf("%s_num_clone", input.case)]] <- unlist(
    lapply(
      checkdf$MID, function(x){
        nrow(input.data[[x]])
      }
    )
  )
}

mid <- "MID4"
olddf <- old.mids[[mid]] %>%
  rowwise() %>%
  mutate(checkACGT = ifelse(
    length(unique(str_split(nSeqCDR3, "")[[1]])) == 4, "yes", "no"
  ))
newdf <- new.mids[[mid]]

intersect.clones <- intersect(newdf$id, olddf$id)

olddf.share <- subset(olddf, olddf$id %in% intersect.clones)
olddf.notshare <- subset(olddf, olddf$id %in% intersect.clones == FALSE)
olddf.check.fail <- subset(olddf, olddf$checkACGT == "no")

newdf.share <- subset(newdf, newdf$id %in% intersect.clones)
newdf.notshare <- subset(newdf, newdf$id %in% intersect.clones == FALSE)

raw.newdf <- subset(clone.obj$clonesets, clone.obj$clonesets$id == mid & clone.obj$clonesets$nSeqCDR3 != "region_not_covered")
raw.newdf <- raw.newdf %>% rowwise() %>%
  mutate(VJnt = sprintf("%s_%s_%s", 
                        str_split(V.gene, "[*]")[[1]][[1]], 
                        str_split(J.gene, "[*]")[[1]][[1]], 
                        nSeqCDR3)) %>%
  mutate(VJ = sprintf("%s_%s", 
                        str_split(V.gene, "[*]")[[1]][[1]], 
                        str_split(J.gene, "[*]")[[1]][[1]]))
raw.newdf.in <-  subset(raw.newdf, raw.newdf$VJnt %in% intersect.clones)
raw.newdf.notin <-  subset(raw.newdf, raw.newdf$VJnt %in% intersect.clones == FALSE)
tmpdf <- subset(raw.newdf, select = c(VJ, VJnt, uniqueMoleculeCount))

# "TGTGCAAGAGAGGGAATTACTACGGTACCATTTGTTTACTGG"
# "TGTGCAAGAGAGGGRATTACTWCGGTACCATTTGTTTACTGG"

olddf.notshare
