gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/bcr_data_analysis/240805_data_analysis"

path.to.rmd1 <- file.path(path.to.project.src, "01_hashtag_antibodies.Rmd")
path.to.rmd2 <- file.path(path.to.project.src, "01_comparing_different_quantile.Rmd")

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "240805_BSimons"
output.version <- "20240820"
config.version <- "default"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis", output.version, config.version)
path.to.save.html <- file.path(path.to.main.output, "html_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
all.quantiles <- c(0.99, 0.95, 0.90, 0.85)
for (sample.id in all.samples){
  for (chosen.quantile in all.quantiles){
    print(sprintf("Working on sample %s, quantile %s", sample.id, chosen.quantile))
    if (file.exists(file.path(path.to.save.html, sprintf("01_hashtag_antibodies_analysis_%s.html", sample.id))) == FALSE){
      rmarkdown::render(input = path.to.rmd1,
                        params = list(
                          outdir = outdir, 
                          PROJECT = PROJECT,
                          output.version = output.version,
                          config.version = config.version, 
                          sample.id = sample.id,
                          chosen.quantile = chosen.quantile
                        ),
                        output_file = sprintf("01_hashtag_antibodies_analysis_%s_quantile_%s.html", sample.id, chosen.quantile),
                        output_dir = path.to.save.html)
    } else {
      print(sprintf("html output for sample %s existed", sample.id))
    }
  }
  
  if (file.exists(file.path(path.to.save.html, sprintf("01_comparing_different_quantile_%s.html", sample.id))) == FALSE){
    rmarkdown::render(input = path.to.rmd2,
                      params = list(
                        outdir = outdir, 
                        PROJECT = PROJECT,
                        output.version = output.version,
                        config.version = config.version, 
                        sample.id = sample.id
                      ),
                      output_file = sprintf("01_comparing_different_quantile_%s.html", sample.id),
                      output_dir = path.to.save.html)
  } else {
    print(sprintf("html output for sample %s existed", sample.id))
  }
}