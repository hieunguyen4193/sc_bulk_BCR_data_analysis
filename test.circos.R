library(circlize)
df <- data.frame(table(clonesets.filtered$id_hashtag, clonesets.filtered$VJ.len.combi))
colnames(df) <- c("SampleID", "Clone", "Count") 

alldf <- list()
newdf <- data.frame()
for (sample.id in unique(df$SampleID)){
  tmpdf <- subset(df, df$SampleID == sample.id)
  sum_COUNT <- sum(tmpdf$Count)
  tmpdf <- tmpdf %>%
    rowwise() %>%
    mutate(Freq = Count/sum_COUNT)
  cur <- 0
  tmpdf <- tmpdf %>% arrange(Freq)
  tmpdf$accumRelCount <- NA
  for (index in 1:nrow(tmpdf)) {
    cur <- cur + tmpdf$Freq[index]
    tmpdf$accumRelCount[index] <- cur
  }
  tmpdf <- subset(tmpdf, tmpdf$Count > 0)
  newdf <- rbind(newdf, tmpdf)
  
  colnames(tmpdf) <- unlist(lapply(
    colnames(tmpdf), function(x){
      if (x != "Clone"){
        return(sprintf("%s_%s", sample.id, x))
      } else {
        return(x)
      }
    }
  ))
  alldf[[sample.id]] <- tmpdf
}

maindf <- alldf[[names(alldf)[[1]]]]
for (sample.id in names(alldf)[2:length(alldf)]){
  maindf <- merge(maindf, alldf[[sample.id]], all.x = TRUE, by.x = "Clone", by.y = "Clone", all.y = TRUE)
}
keep.cols <- to_vec(
  for (item in colnames(maindf)){
    if(grepl("Freq", item) == TRUE | grepl("accumRelCount", item) == TRUE){
      item
    }
  }
)
maindf <- maindf[, c("Clone", keep.cols)]
maindf[is.na(maindf)] <- 0

circos.par(cell.padding = c(0, 0, 0, 0), gap.degree = 5, track.height = 0.1, start.degree = -2.5, points.overflow.warning = FALSE)
circos.initialize(factors = newdf$SampleID, x = newdf$accumRelCount)
circos.track(
  factors = newdf$SampleID, x = newdf$Freq, y = newdf$accumRelCount, bg.border = NA,
  panel.fun = function(x, y) {
    circos.lines(c(CELL_META$xlim[1], CELL_META$xlim[2]), c((CELL_META$ylim[1] + CELL_META$ycenter) * 0.5, (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
    for (index in 1:length(x)) {
      circos.lines(c(y[index], y[index]), c(CELL_META$ylim[1], (CELL_META$ylim[1] + CELL_META$ycenter) * 0.5))
    }
    n <- length(x) - 1
  }
)

all.combi <- combn(names(alldf), 2)
for (i in seq(1, ncol(all.combi))){
  sample1 <- all.combi[1, i]
  sample2 <- all.combi[2, i]
  
  plotdf <- maindf[, c("Clone",
                       sprintf("%s_Freq", sample1), sprintf("%s_accumRelCount", sample1),
                       sprintf("%s_Freq", sample2), sprintf("%s_accumRelCount", sample2))]
  colnames(plotdf) <- c("Clone", "freq1", "accumFreq1", "freq2", "accumFreq2")
  plotdf <- subset(plotdf, plotdf$freq1 > 0 & plotdf$freq2 > 0)
  for (j in seq(1, nrow(plotdf))){
    sample1.freq <- plotdf[j, ]$freq1
    sample2.freq <- plotdf[j, ]$freq2
    sample1.accumFreq <- plotdf[j, ]$accumFreq1
    sample2.accumFreq <- plotdf[j, ]$accumFreq2
    circos.link(sample1, c(sample1.accumFreq - sample1.freq, sample1.accumFreq), sample2, c(sample2.accumFreq - sample2.freq, sample2.accumFreq))  
  }
}


circos.clear()
