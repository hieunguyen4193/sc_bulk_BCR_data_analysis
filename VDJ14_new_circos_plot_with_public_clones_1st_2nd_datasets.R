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

# sc.projects <- c("240805_BSimons_filterHT_cluster_renamed")
sc.projects <- c("1st_dataset", "2nd_dataset")
# sc.projects <- c("2nd_dataset")
# sc.projects <- c("241002_241104_BSimons")

bulk.projects <- c("240826_BSimons")
list.of.PROJECT <- c(sc.projects, bulk.projects)

thres <- 0.85 

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
if (length(sc.projects) > 1){
  path.to.14.output <- file.path(outdir, "VDJ_output", "14_output_20250323", paste(sc.projects, collapse = "_"))  
}

dir.create(path.to.14.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### public clones
#####----------------------------------------------------------------------#####
raw.public.clonedf <- read.csv(file.path(path.to.main.src, "public_clones.csv")) %>%
  rowwise() %>%
  mutate(id = sprintf("%s_%s", V_gene, J_gene))

public.clonedf <- data.frame(
  id = raw.public.clonedf$id
)
public.clonedf$cloneCount <- 1
public.clonedf$Freq <- 1/nrow(public.clonedf)
accum.sum <- 0
all.accum.sum <- list()
for (i in seq(1, nrow(public.clonedf))){
  accum.sum <- public.clonedf[i, ][["Freq"]] + accum.sum
  all.accum.sum <- c(all.accum.sum, c(accum.sum))
}
public.clonedf$accum.Freq <- unlist(all.accum.sum)
public.clonedf$SampleID <- "publicClone"
public.clonedf <- rbind(data.frame(cloneCount = 0, 
                                   id = "START", 
                                   SampleID = "publicClone", 
                                   Freq = 0, 
                                   accum.Freq = 0), 
                        public.clonedf)

#####----------------------------------------------------------------------#####
##### get the clone information data
#####----------------------------------------------------------------------#####
all.files <- list()
for (input.sc.project in sc.projects){
  path.to.tree.05.output <- file.path(outdir, 
                                      "GEX_output", 
                                      "05_output", 
                                      input.sc.project, 
                                      "input_circos_public_clones")
  
  tmp.all.files <- Sys.glob(file.path(path.to.tree.05.output, "*.csv"))
  names(tmp.all.files) <- to_vec(
    for (item in tmp.all.files){
      str_replace(basename(item), "_publicClone.csv", "")
    }
  )
  all.files <- c(all.files, tmp.all.files)
}

names(all.files) <- unlist(
  lapply(names(all.files), function(x){
    x.split <- str_split(x, "_")[[1]]
    half.len <- length(x.split)/2
    output <- paste0(x.split[1:half.len], collapse = "_")
    return(output)
  })
)
#####----------------------------------------------------------------------#####
##### main analysis: generate clone count df
#####----------------------------------------------------------------------#####

cloneCountdf <- data.frame()
match.public.clones <- list()
all.cloneCountdf <- list()

for (sample.name in names(all.files)){
  print(sample.name)
  tmp.clonedf <- read.csv(all.files[[sample.name]])
  tmp.clonedf.raw <- tmp.clonedf
  colnames(tmp.clonedf) <- c("id", "cloneCount", "match_public_clone")
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
  tmp.clonedf <- subset(tmp.clonedf, select = c(id, cloneCount, Freq, accum.Freq))
  tmp.clonedf$SampleID <- sample.name
  tmp.clonedf <- rbind(data.frame(cloneCount = 0, id = "START", SampleID = sample.name, Freq = 0, accum.Freq = 0), tmp.clonedf)
  cloneCountdf <- rbind(cloneCountdf, tmp.clonedf)
  match.public.clones[[sample.name]] <- subset(tmp.clonedf.raw, tmp.clonedf.raw$match_public_clone == "yes")$cloneID
  all.cloneCountdf[[sample.name]] <- tmp.clonedf
}

cloneCountdf <- rbind(cloneCountdf, public.clonedf)

##### generate file alias
fileAliases <- names(all.files)
names(fileAliases) <- names(all.files)
fileAliases <- c(fileAliases, list(publicClone = "publicClone"))

#####----------------------------------------------------------------------#####
##### DEFINE COLORS
#####----------------------------------------------------------------------#####
linkColor1 = "#FF000080"
linkColor2 = "lightgray"
countColors <- c("#FFFFFFFF", "#0000FFFF")
maxCount <- cloneCountdf$Freq %>% max()
countRamp <- function(x) {
  ramp <- colorRamp(c(countColors[1], countColors[2]), alpha = TRUE)
  color <- ramp(x / maxCount)
  rgb(color, alpha = color[4], maxColorValue = 255)
}

##### Define LINK colors
linkColors <- c(linkColor1, linkColor1)
linkRampAB <- function(x, maxLinkSize) {
  ramp <- colorRamp(c(linkColors[1], linkColors[2]), alpha = TRUE)
  color <- ramp(x / maxLinkSize)
  rgb(color, alpha = color[4], maxColorValue = 255)
}

##### Define BACKGROUND LINK colors
linkColors2 <- c(linkColor2, linkColor2)
linkRampAB2 <- function(x, maxLinkSize) {
  ramp <- colorRamp(c(linkColors2[1], linkColors2[2]), alpha = TRUE)
  color <- ramp(x / maxLinkSize)
  rgb(color, alpha = color[4], maxColorValue = 255)
}

#####----------------------------------------------------------------------#####
path.to.save.svg <- file.path(path.to.14.output, sprintf("%s.svg", paste0(sc.projects, collapse = "_") ))

##### prepare the CIRCOS plot

cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = names(fileAliases))

svg(path.to.save.svg)
circos.par(cell.padding = c(0, 0, 0, 0), 
           gap.degree = 5, 
           track.height = 0.1, 
           start.degree = -2.5, 
           points.overflow.warning = FALSE)

# the public clone (10th component) gets half of the circle, 
# while the other components get the other half.
# custom.sector.width <- c(rep(0.5/9, 9), 0.5)

# or 1/3.
custom.sector.width <- c(rep(0.7/9, 9), 0.3)

circos.initialize(factors = cloneCountdf$SampleID, x = cloneCountdf$accum.Freq, sector.width = custom.sector.width)

##### generating base structure of the circos plot ...
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

for (sample.name in names(all.files)){
  if(length(match.public.clones[[sample.name]]) != 0){
    tmp.plotdf <- all.cloneCountdf[[sample.name]] %>% 
      subset(id %in% match.public.clones[[sample.name]]) %>%
      rowwise() %>%
      mutate(VJ = paste0(str_split(id, "_")[[1]][1:2], collapse = "_"))
    
    for (index in seq(1, nrow(tmp.plotdf))){
      tmp.cloneID <- tmp.plotdf[index,]$VJ
      tmp.public.clonedf <- subset(public.clonedf, public.clonedf$id == tmp.cloneID)
      maxLinkSize <- max(tmp.plotdf$Freq + public.clonedf$Freq)
      
      circos.link(sample.name, 
                  c(tmp.plotdf$accum.Freq[index] - tmp.plotdf$Freq[index], 
                    tmp.plotdf$accum.Freq[index]),
                  "publicClone", 
                  c(tmp.public.clonedf$accum.Freq - tmp.public.clonedf$Freq, 
                    tmp.public.clonedf$accum.Freq), 
                  col = linkRampAB( tmp.plotdf$Freq[index] + tmp.public.clonedf$Freq, maxLinkSize),
                  border = NA)
    }
  }
}
circos.clear()
dev.off()
