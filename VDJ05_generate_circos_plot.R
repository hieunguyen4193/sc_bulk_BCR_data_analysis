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
path.to.VDJ.output <- file.path( outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> NEW DATA
#####----------------------------------------------------------------------#####
path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", PROJECT)
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, "circos_plot", "Hieu_version"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, "circos_plot", "Fabio_version"), showWarnings = FALSE, recursive = TRUE)

path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
path.to.save.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")

dir.create(file.path(path.to.save.output, "single_MID_clone_df"), showWarnings = FALSE, recursive = TRUE)

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

for (mid in unique(mid.metadata$X)){
  if (file.exists(file.path(path.to.save.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid))) == FALSE){
    print(sprintf("Working on sample MID %s", mid))
    clonedf <- subset(full.clonedf, full.clonedf$id == mid)
    input.circos <- subset(clonedf, select = c(V.gene, 
                                               J.gene, 
                                               aaSeqCDR3, 
                                               nSeqCDR3, 
                                               uniqueMoleculeCount)) %>%
      rowwise() %>%
      mutate(VJnt = sprintf("%s_%s_%s", V.gene, J.gene, nSeqCDR3)) %>%
      group_by(VJnt) %>%
      summarise(cloneCount = sum(uniqueMoleculeCount)) %>% 
      arrange(desc(cloneCount))
    colnames(input.circos) <- c("id", "cloneCount")
    input.circos <- input.circos %>% rowwise() %>%
      mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
      mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
      mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]])
    write.table(input.circos, file.path(path.to.save.output, "single_MID_clone_df", sprintf("%s.simplified.csv", mid)), quote = FALSE, sep = "\t", row.names = FALSE)  
  }
}

all.mid.files <- Sys.glob(file.path(path.to.save.output, "single_MID_clone_df", "*.simplified.csv"))
names(all.mid.files) <- to_vec(
  for (item in all.mid.files){
    str_replace(basename(item), ".simplified.csv", "")
  }
)

#####----------------------------------------------------------------------#####
##### GENERATE CIRCOS PLOT WITH FABIO SCRIPT
#####----------------------------------------------------------------------#####
source(file.path(path.to.main.src, "circos.R"))
count.mice <- table(mid.metadata$mouse)

# for (mouse.id in count.mice[count.mice >= 2] %>% names()){
#   selected.mids <- yfp.mids[[mouse.id]]$all_w_biopsy
#   fileAliases <- to_vec( for (item in selected.mids){
#     subset(mid.metadata, mid.metadata$X == item)$population
#   })
#   cutoff <- 1.0
#   files <- all.mid.files[selected.mids]
#   input.sort <- FALSE
#   countColors <- c("#FFFFFFFF", "#0000FFFF")
#   linkColors <- rep("#FF000080", ifelse(length(files) == 1, 1, length(combn(length(files), 2))))
#   dir.create(file.path(path.to.save.output, "circos_plots_FT_script", mouse.id), showWarnings = FALSE, recursive = TRUE)
#   dir.create(file.path(path.to.save.output, "circos_plots_Hieu_script", mouse.id), showWarnings = FALSE, recursive = TRUE)
#   
#   circos(files = files,
#          fileAliases = fileAliases,
#          saveFolder = file.path(path.to.save.output, "circos_plots_FT_script", sprintf("%s/", mouse.id)),
#          cutoff = cutoff,
#          sort = input.sort,
#          countColors = countColors,
#          linkColors = linkColors,
#          showLinks = rep(TRUE, length(combn(length(files), 2))))
#   
#   
# }

#####----------------------------------------------------------------------#####
##### GENERATE CIRCOS PLOT WITH HIEU SCRIPT
#####----------------------------------------------------------------------#####
mouse.id <- "m11"
input.case <- "all_w_biopsy"
selected.mids <- yfp.mids[[mouse.id]][[input.case]]
filter.clone <- FALSE
filter.clone.cutoff <- 10
saveFileName <- sprintf("%s_%s_circos.svg", mouse.id, paste(selected.mids, collapse = "_"))

if (filter.clone == TRUE){
  path.to.save.svg <- file.path(path.to.05.output, 
                                "circos_plot", 
                                "Hieu_version", 
                                sprintf("filter_clone_%s",filter.clone), 
                                sprintf("filter_cutoff_%s", filter.clone.cutoff),
                                saveFileName)
  dir.create(file.path(path.to.05.output, 
                       "circos_plot", 
                       "Hieu_version", 
                       sprintf("filter_clone_%s",filter.clone),
                       sprintf("filter_cutoff_%s", filter.clone.cutoff),),
             showWarnings = FALSE, 
             recursive = TRUE)
} else {
  path.to.save.svg <- file.path(path.to.05.output, 
                                "circos_plot", 
                                "Hieu_version", 
                                sprintf("filter_clone_%s",filter.clone), 
                                saveFileName)
  dir.create(file.path(path.to.05.output, 
                       "circos_plot", 
                       "Hieu_version", 
                       sprintf("filter_clone_%s",filter.clone)),
             showWarnings = FALSE, 
             recursive = TRUE)
}


cloneCountdf <- data.frame()
for (input.mid in selected.mids){
  tmp.clonedf <- read_tsv(all.mid.files[input.mid])
  if (filter.clone == TRUE){
    tmp.clonedf <- subset(tmp.clonedf, tmp.clonedf$cloneCount > filter.clone.cutoff)
  }
  totalCloneCount <- tmp.clonedf$cloneCount %>% sum()
  tmp.clonedf <- tmp.clonedf %>%
    rowwise() %>%
    mutate(Freq = cloneCount/totalCloneCount) %>%
    arrange(Freq)
  accum.sum <- 0
  all.accum.sum <- list()
  for (i in seq(1, nrow(tmp.clonedf))){
    accum.sum <- tmp.clonedf[i, ][["Freq"]] + accum.sum
    all.accum.sum <- c(all.accum.sum, c(accum.sum))
  }
  tmp.clonedf$accum.Freq <- unlist(all.accum.sum)
  tmp.clonedf <- subset(tmp.clonedf, select = c(id, Freq, accum.Freq))
  tmp.clonedf$SampleID <- input.mid
  cloneCountdf <- rbind(cloneCountdf, tmp.clonedf)
}
cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = unique(cloneCountdf$SampleID))

plotdf <- subset(cloneCountdf, cloneCountdf$SampleID == levels(cloneCountdf$SampleID)[[1]]) %>%
  subset(select = c(id, Freq, accum.Freq))
colnames(plotdf) <- c("id", 
                      sprintf("%s_Freq", levels(cloneCountdf$SampleID)[[1]]), 
                      sprintf("%s_accumFreq", levels(cloneCountdf$SampleID)[[1]]))

for (sample.id in levels(cloneCountdf$SampleID)[2: length(levels(cloneCountdf$SampleID))]){
  tmpdf <- subset(cloneCountdf, cloneCountdf$SampleID == sample.id) %>%
    subset(select = c(id, Freq, accum.Freq))
  colnames(tmpdf) <- c("id", 
                        sprintf("%s_Freq", sample.id), 
                        sprintf("%s_accumFreq", sample.id))
  plotdf <- merge(plotdf, tmpdf, by.x = "id", by.y = "id", all.x = TRUE, all.y = TRUE)
}
plotdf[is.na(plotdf)] <- 0

fileAliases <- to_vec(
  for (item in levels(cloneCountdf$SampleID)){
    subset(mid.metadata, mid.metadata$X == item)$population
  }
)

##### Define COUNT colors
countColors <- c("#FFFFFFFF", "#0000FFFF")
maxCount <- cloneCountdf$Freq %>% max()
countRamp <- function(x) {
  ramp <- colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
  color <- ramp(x / maxCount)
  rgb(color, alpha = color[4], maxColorValue = 255)
}

##### Define LINK colors
linkColors <- c("#FF000080", "#FF000080")
linkRampAB <- function(x, maxLinkSize) {
  ramp <- colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
  color <- ramp(x / maxLinkSize)
  rgb(color, alpha = color[4], maxColorValue = 255)
}
  
##### initialize the structure of CIRCOS plot. 
svg(path.to.save.svg)
circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5, track.height = 0.1, start.degree = -2.5, points.overflow.warning = FALSE)
circos.initialize(factors = cloneCountdf$SampleID, x = cloneCountdf$accum.Freq)

# Add ring with labels, ticks and colored boxes
circos.track(
  factors = cloneCountdf$SampleID, 
  x = cloneCountdf$Freq, 
  y = cloneCountdf$accum.Freq, 
  bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, 
                CELL_META$ylim[2] + mm_y(4), 
                fileAliases[CELL_META$sector.numeric.index], 
                niceFacing = TRUE)
    circos.lines(c(CELL_META$xlim[1], 
                   CELL_META$xlim[2]), 
                 c((CELL_META$ylim[1] + CELL_META$ycenter) * 0.5, 
                   (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
    for (index in 1:length(x)) {
      circos.lines(c(y[index], 
                     y[index]), 
                   c(CELL_META$ylim[1], (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
    }
    n <- length(x) - 1
    circos.rect(y[2:(n + 1)] - x[2:(n + 1)], 
                rep(CELL_META$ycenter, n), 
                y[2:(n + 1)], 
                rep(CELL_META$ylim[2], n), 
                col = countRamp(x[2:(n + 1)]), 
                border = NA)
  }
)

all.combi <- combn(levels(cloneCountdf$SampleID), 2)
for (j in seq(1, ncol(all.combi))){
  sample1 <- all.combi[, j][[1]]
  sample2 <- all.combi[, j][[2]]
  
  tmp.plotdf <- plotdf[, c( 
    sprintf("%s_Freq", sample1), sprintf("%s_accumFreq", sample1),
    sprintf("%s_Freq", sample2), sprintf("%s_accumFreq", sample2)
  )]
  colnames(tmp.plotdf) <- c("sample1_Freq", "sample1_accumFreq",
                            "sample2_Freq", "sample2_accumFreq")
  tmp.plotdf <- subset(tmp.plotdf, tmp.plotdf$sample1_Freq > 0 & tmp.plotdf$sample2_Freq > 0)
  maxLinkSize <- max(tmp.plotdf$sample1_Freq + tmp.plotdf$sample2_Freq)
  
  for (index in seq(1, nrow(tmp.plotdf))){
    circos.link(sample1, 
                c(tmp.plotdf$sample1_accumFreq[index] - tmp.plotdf$sample1_Freq[index], 
                  tmp.plotdf$sample1_accumFreq[index]),
                sample2, 
                c(tmp.plotdf$sample2_accumFreq[index] - tmp.plotdf$sample2_Freq[index], 
                  tmp.plotdf$sample2_accumFreq[index]), 
                col = linkRampAB(tmp.plotdf$sample1_Freq[index] + tmp.plotdf$sample2_Freq[index], maxLinkSize),
                border = NA)
  }  
}
# Close diagramm and output file
circos.clear()
dev.off()

  


