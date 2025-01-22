generate_circos <- function(
    input.files,
    fileAliases,
    saveFileName,
    outputdir,
    filter.clone = FALSE,
    filter.clone.cutoff = NULL,
    linkColor1 = "#FF000080",
    linkColor2 = "lightgray",
    group.to.highlight1 = NULL,
    group.to.highlight2 = NULL,
    ordered.samples = NULL
){
  if (filter.clone == TRUE){
    path.to.save.svg <- file.path(outputdir, 
                                  sprintf("filter_clone_%s",filter.clone), 
                                  sprintf("filter_cutoff_%s", filter.clone.cutoff),
                                  saveFileName)
    dir.create(file.path(outputdir, 
                         sprintf("filter_clone_%s",filter.clone),
                         sprintf("filter_cutoff_%s", filter.clone.cutoff)),
               showWarnings = FALSE, 
               recursive = TRUE)
  } else {
    path.to.save.svg <- file.path(outputdir, 
                                  sprintf("filter_clone_%s",filter.clone), 
                                  saveFileName)
    dir.create(file.path(outputdir, 
                         sprintf("filter_clone_%s", filter.clone)),
               showWarnings = FALSE, 
               recursive = TRUE)
  }
  
  
  cloneCountdf <- data.frame()
  for (input.mid in names(input.files)){
    print(sprintf("reading in clones from sample %s", input.mid))
    tmp.clonedf <- read_tsv(input.files[[input.mid]])
    print(dim(tmp.clonedf))
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
    tmp.clonedf <- subset(tmp.clonedf, select = c(id, cloneCount, Freq, accum.Freq))
    tmp.clonedf$SampleID <- input.mid
    cloneCountdf <- rbind(cloneCountdf, tmp.clonedf)
  }
  
  count.clone.in.samples <- table(cloneCountdf$SampleID)
  exclude.samples <- count.clone.in.samples[count.clone.in.samples <= 1] %>% names()
  cloneCountdf <- subset(cloneCountdf, cloneCountdf$SampleID %in% exclude.samples == FALSE)
  keep.samples <- setdiff(names(input.files), exclude.samples)
  if (is.null(ordered.samples) == TRUE){
    cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = keep.samples)
  } else {
    print(sprintf("Samples are ordered based on input order %s", paste(ordered.samples, collapse = ", ")))
    cloneCountdf$SampleID <- factor(cloneCountdf$SampleID, levels = ordered.samples)    
  }

  
  new.fileAliases <- to_vec(
    for (item in keep.samples){
      fileAliases[[item]]
    }
  )
  fileAliases <- new.fileAliases
  plotdf <- subset(cloneCountdf, cloneCountdf$SampleID == levels(cloneCountdf$SampleID)[[1]]) %>%
    subset(select = c(id, cloneCount, Freq, accum.Freq))
  colnames(plotdf) <- c("id", 
                        sprintf("%s_Count", levels(cloneCountdf$SampleID)[[1]]),
                        sprintf("%s_Freq", levels(cloneCountdf$SampleID)[[1]]), 
                        sprintf("%s_accumFreq", levels(cloneCountdf$SampleID)[[1]]))
  
  for (sample.id in levels(cloneCountdf$SampleID)[2: length(levels(cloneCountdf$SampleID))]){
    tmpdf <- subset(cloneCountdf, cloneCountdf$SampleID == sample.id) %>%
      subset(select = c(id, cloneCount, Freq, accum.Freq))
    colnames(tmpdf) <- c("id", 
                         sprintf("%s_Count", sample.id),
                         sprintf("%s_Freq", sample.id), 
                         sprintf("%s_accumFreq", sample.id))
    plotdf <- merge(plotdf, tmpdf, by.x = "id", by.y = "id", all.x = TRUE, all.y = TRUE)
  }
  plotdf[is.na(plotdf)] <- 0
  
  write.csv(plotdf, str_replace(path.to.save.svg, ".svg", ".csv"))
  ##### Define COUNT colors
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
  
  ##### initialize the structure of CIRCOS plot. 
  svg(path.to.save.svg)
  circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5, track.height = 0.1, start.degree = -2.5, points.overflow.warning = FALSE)
  circos.initialize(factors = cloneCountdf$SampleID, x = cloneCountdf$accum.Freq)
  
  print("generating base structure of the circos plot ...")
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
  all.combidf <- data.frame(all.combi %>% t)
  colnames(all.combidf) <- c("Sample1", "Sample2")
  all.combidf <- all.combidf %>% rowwise() %>%
    mutate(check = ifelse(Sample1 %in% group.to.highlight1 & Sample2 %in% group.to.highlight2, "yes", "no")) %>%
    arrange(check)
  for (j in seq(1, nrow(all.combidf))){
    sample1 <- all.combidf[j, ][["Sample1"]]
    sample2 <- all.combidf[j, ][["Sample2"]]
    print(sprintf("Generating links between %s and %s", sample1, sample2))
    tmp.plotdf <- plotdf[, c( 
      sprintf("%s_Freq", sample1), sprintf("%s_accumFreq", sample1),
      sprintf("%s_Freq", sample2), sprintf("%s_accumFreq", sample2)
    )]
    colnames(tmp.plotdf) <- c("sample1_Freq", "sample1_accumFreq",
                              "sample2_Freq", "sample2_accumFreq")
    tmp.plotdf <- subset(tmp.plotdf, tmp.plotdf$sample1_Freq > 0 & tmp.plotdf$sample2_Freq > 0)
    if (nrow(tmp.plotdf) != 0){
      maxLinkSize <- max(tmp.plotdf$sample1_Freq + tmp.plotdf$sample2_Freq)
      
      for (index in seq(1, nrow(tmp.plotdf))){
        if (sample1 %in% group.to.highlight1 & sample2 %in% group.to.highlight2){
          circos.link(sample1, 
                      c(tmp.plotdf$sample1_accumFreq[index] - tmp.plotdf$sample1_Freq[index], 
                        tmp.plotdf$sample1_accumFreq[index]),
                      sample2, 
                      c(tmp.plotdf$sample2_accumFreq[index] - tmp.plotdf$sample2_Freq[index], 
                        tmp.plotdf$sample2_accumFreq[index]), 
                      col = linkRampAB(tmp.plotdf$sample1_Freq[index] + tmp.plotdf$sample2_Freq[index], maxLinkSize),
                      border = NA)
        } else {
          circos.link(sample1, 
                      c(tmp.plotdf$sample1_accumFreq[index] - tmp.plotdf$sample1_Freq[index], 
                        tmp.plotdf$sample1_accumFreq[index]),
                      sample2, 
                      c(tmp.plotdf$sample2_accumFreq[index] - tmp.plotdf$sample2_Freq[index], 
                        tmp.plotdf$sample2_accumFreq[index]), 
                      col = linkRampAB2(tmp.plotdf$sample1_Freq[index] + tmp.plotdf$sample2_Freq[index], maxLinkSize),
                      border = NA)
        }

      }  
    }
  }
  # Close diagramm and output file
  circos.clear()
  dev.off()
}
