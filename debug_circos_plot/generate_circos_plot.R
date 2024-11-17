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
##### READ CLONE DATA
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


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
path.to.old.data <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets/test_circos_inputs"

mouse.id <- "m11"
selected.mids <- sample.list[[mouse.id]]

clonedf <- full.clonedf %>% subset(id %in% selected.mids)

simplified.clonedf <- clonedf %>% subset(select = c(id, V.gene, J.gene, aaSeqCDR3, nSeqCDR3))
simplified.clonedf$mouse.id <- mouse.id
simplified.clonedf <- simplified.clonedf %>% rowwise() %>%
  mutate(V.gene = str_split(V.gene, "[*]")[[1]][[1]]) %>%
  mutate(J.gene = str_split(J.gene, "[*]")[[1]][[1]]) %>%
  mutate(clone = sprintf("%s_%s_%s", V.gene, J.gene, nSeqCDR3)) %>%
  mutate(group = subset(mid.metadata, mid.metadata$X == id)$population)

# count.clonedf <- table(simplified.clonedf$clone, simplified.clonedf$group) %>% 
#   as.data.frame() %>%
#   pivot_wider(names_from = "Var2", values_from = "Freq")


count.clonedf <- table(simplified.clonedf$clone, simplified.clonedf$group) %>%
  as.data.frame() %>%
  rowwise() %>% 
  mutate(Freq2 = Freq/nrow(subset(simplified.clonedf, simplified.clonedf$group == Var2))) 

maindf <- data.frame(clone = unique(count.clonedf$Var1))
newdf <- data.frame()
available.groups <- unique(count.clonedf$Var2)

for (g in unique(count.clonedf$Var2)){
  tmpdf <- subset(count.clonedf, count.clonedf$Var2 == g) %>% arrange(Freq)
  print(dim(tmpdf))
  all.accum.sum <- c()
  accum.sum <- 0
  for (i in seq(1, nrow(tmpdf))){
    c <- tmpdf[i, ][["Freq2"]]
    accum.sum <- accum.sum + c
    all.accum.sum <- c(all.accum.sum, accum.sum)
  }
  tmpdf[["accumRelCount"]] <- all.accum.sum 
  newdf <- rbind(newdf, subset(tmpdf, select = c(Var1, Var2, Freq2, accumRelCount)))
  tmpdf <- subset(tmpdf, select = -c(Freq, Var2))
  colnames(tmpdf) <- c("clone", sprintf("%s_Freq", g), sprintf("%s_accumRelCount", g))
  maindf <- merge(maindf, tmpdf, by.x = "clone", by.y = "clone")
}

for (g in setdiff(available.groups, c("biopsy"))){
  num.share <- subset(maindf, maindf[["biopsy_Freq"]] !=0 & maindf[[sprintf("%s_Freq", g)]] !=0) %>% nrow()
  print(sprintf("number of share clones: %s %s = %s", "biopsy", g, num.share))
}

colnames(newdf) <- c("clone", "SampleID", "Freq", "accumRelCount")
newdf$SampleID <- factor(newdf$SampleID, levels = available.groups)


countColors <- c("#FFFFFFFF", "#0000FFFF")
maxCount <- max(newdf$Freq)
countRamp <- function(x) {
  ramp <- colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
  color <- ramp(x / maxCount)
  rgb(color, alpha = color[4], maxColorValue = 255)
}
linkColors <- c("#FF000080", "#FF000080")

circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5, track.height = 0.1, start.degree = -2.5, points.overflow.warning = FALSE)
circos.initialize(factors = newdf$SampleID, x = newdf$accumRelCount)
circos.track(
  factors = newdf$SampleID, x = newdf$Freq, y = newdf$accumRelCount, bg.border = NA,
  panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(4), levels(newdf$SampleID)[CELL_META$sector.numeric.index], niceFacing = TRUE)
    circos.lines(c(CELL_META$xlim[1], CELL_META$xlim[2]), c((CELL_META$ylim[1] + CELL_META$ycenter) * 0.5, (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
    for (index in 1:length(x)) {
      circos.lines(c(y[index], y[index]), c(CELL_META$ylim[1], (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
    }
    n <- length(x) - 1
    circos.rect(y[2:(n + 1)] - x[2:(n + 1)], rep(CELL_META$ycenter, n), y[2:(n + 1)], rep(CELL_META$ylim[2], n), col = countRamp(x[2:(n + 1)]), border = NA)
  }
)

all.combi <- combn(levels(newdf$SampleID), 2)

for (i in seq(1, ncol(all.combi))){
  col <- "#FF000080"
  sample1 <- all.combi[1, i]
  sample2 <- all.combi[2, i]
  print(sprintf("%s vs %s", sample1, sample2))
  plotdf <- maindf[, c("clone",
                       sprintf("%s_Freq", sample1), sprintf("%s_accumRelCount", sample1),
                       sprintf("%s_Freq", sample2), sprintf("%s_accumRelCount", sample2))]
  colnames(plotdf) <- c("Clone", "freq1", "accumFreq1", "freq2", "accumFreq2")
  plotdf <- subset(plotdf, plotdf$freq1 > 0 & plotdf$freq2 > 0)
  if (nrow(plotdf) != 0){
    for (j in seq(1, nrow(plotdf))){
      sample1.freq <- plotdf[j, ]$freq1
      sample2.freq <- plotdf[j, ]$freq2
      sample1.accumFreq <- plotdf[j, ]$accumFreq1
      sample2.accumFreq <- plotdf[j, ]$accumFreq2
      circos.link(sample1, 
                  c(sample1.accumFreq - sample1.freq, 
                    sample1.accumFreq), 
                  sample2, 
                  c(sample2.accumFreq - sample2.freq, 
                    sample2.accumFreq),
                  col = col)  
    }
  }
}
circos.clear()