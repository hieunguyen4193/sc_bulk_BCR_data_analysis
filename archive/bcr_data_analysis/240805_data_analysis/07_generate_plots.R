gc()
rm(list = ls())

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

path.to.project.src <- "/media/hieunguyen/HNSD01/src/bcr_data_analysis/240805_data_analysis"
source(file.path(path.to.project.src, "config.R"))

outdir <- "/media/hieunguyen/HNSD01/outdir"
PROJECT <- "240805_BSimons"
output.version <- "20240820"
config.version <- "default"
chosen.quantile <- 0.85
integration.case <- "all_samples"

path.to.main.input <- file.path(outdir, PROJECT, output.version, config.version)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis", output.version, config.version)

path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.04.output <- file.path(path.to.main.output, "04_output", sprintf("quantile_%s", chosen.quantile))
path.to.05.output <- file.path(path.to.main.output, "05_output", sprintf("quantile_%s", chosen.quantile), integration.case)
path.to.06.output <- file.path(path.to.main.output, "06_output", sprintf("quantile_%s", chosen.quantile), integration.case)
path.to.07.output <- file.path(path.to.main.output, "07_output", sprintf("quantile_%s", chosen.quantile), integration.case)
dir.create(path.to.07.output, showWarnings = FALSE, recursive = TRUE)

clonedf <- readxl::read_excel(file.path(path.to.06.output, "full_clonedf.xlsx"))

all.cases <- list(
  "m1_only" = c("m1"),
  "m2_only" = c("m2"),
  "m3_only" = c("m3"),
  "all" = c("m1", "m2", "m3")
)

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

# for (plot.case in names(all.cases)){
for (plot.case in c("m1_only")){
  dir.create(file.path(path.to.07.output, plot.case), showWarnings = FALSE, recursive = TRUE)
  mouse.list <- all.cases[[plot.case]]
  
  plot.clonedf <- subset(clonedf, clonedf$mouseid %in% mouse.list)
  clone.sampledf <- data.frame(SampleID = unique(plot.clonedf$MID))
  
  for (sample2 in unique(plot.clonedf$MID)){
    clone.sampledf[[sample2]] <- unlist(
      lapply(
        clone.sampledf$SampleID, function(sample1){
          return(length(unique(intersect(subset(plot.clonedf, plot.clonedf$MID == sample1)$VJcombi_CDR3_0.85, subset(plot.clonedf, plot.clonedf$MID == sample2)$VJcombi_CDR3_0.85))))
        }
      )
    )
  }
  
  clone.sampledf.pivot <- clone.sampledf %>% pivot_longer(!SampleID, names_to = "SampleID2", values_to = "count")
  writexl::write_xlsx(clone.sampledf, file.path(path.to.07.output, plot.case, "count_clone_shared_between_samples.xlsx"))
  
  writexl::write_xlsx(plot.clonedf, file.path(path.to.07.output, plot.case, sprintf("all_clone_information_in_%s.xlsx", plot.case)))
  
  heatmap.plotdf <- clone.sampledf %>% column_to_rownames("SampleID") %>% as.matrix()
  heatmap.plotdf.scaled <- (heatmap.plotdf - rowMeans(heatmap.plotdf))/rowSds(heatmap.plotdf)
  heatmap.plotdf.scaled.pivot <- data.frame(heatmap.plotdf.scaled) %>% rownames_to_column("SampleID") %>%
    pivot_longer(!SampleID, names_to = "Sample2", values_to = "count")
  
  heatmap.scaled <- heatmap.plotdf.scaled.pivot %>% ggplot(aes(x = SampleID, y = Sample2, fill = count)) + geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    coord_fixed()
  ggsave(plot = heatmap.scaled, filename = sprintf("heatmap_scaled_shared_clones_%s.svg", plot.case),
         path = file.path(path.to.07.output, plot.case), device = "svg", height = 10, width = 10, dpi = 300)
  
  heatmap.count <- clone.sampledf.pivot %>% ggplot(aes(x = SampleID, y = SampleID2, fill = count)) + geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    coord_fixed() + 
    geom_text(aes(label = count), color = "black", size = 4)
  ggsave(plot = heatmap.scaled, filename = sprintf("heatmap_raw_shared_clones_%s.svg", plot.case),
         path = file.path(path.to.07.output, plot.case), device = "svg", height = 10, width = 10, dpi = 300)
  
  svg(file.path(path.to.07.output, plot.case, sprintf("chordDiagram_%s.svg", plot.case)))
  chordDiagram(clone.sampledf.pivot, self.link = FALSE)
  dev.off()
}

