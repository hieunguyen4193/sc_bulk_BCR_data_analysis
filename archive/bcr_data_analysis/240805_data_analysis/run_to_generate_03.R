gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/bcr_data_analysis/240805_data_analysis"

path.to.rmd1 <- file.path(path.to.project.src, "03_downstream_analysis_singlet_only.Rmd")

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "240805_BSimons"
output.version <- "20240820"
config.version <- "default"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis", output.version, config.version)
path.to.save.html <- file.path(path.to.main.output, "html_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
all.quantiles <- c(0.99, 0.95, 0.90, 0.85)



for (chosen.quantile in all.quantiles){
  for (sample.id in all.samples){
    print("--------------------------------------------------------------------")
    print(sprintf("%s --  %s", chosen.quantile, sample.id))
    print("--------------------------------------------------------------------")
    
    if (file.exists(file.path(path.to.save.html, sprintf("03_downstream_analysis_singlet_only_quantile_%s_sample_%s.html", chosen.quantile, sample.id))) == FALSE){
      rmarkdown::render(input = path.to.rmd1,
                        params = list(
                          outdir = outdir,
                          PROJECT = PROJECT,
                          output.version = output.version,
                          config.version = config.version,
                          sample.id = sample.id,
                          chosen.quantile = chosen.quantile
                        ),
                        output_file = sprintf("03_downstream_analysis_singlet_only_quantile_%s_sample_%s.html", chosen.quantile, sample.id),
                        output_dir = path.to.save.html)
    } else {
      print(sprintf("Html output file for sample %s, quantile %s existed", sample.id, chosen.quantile))
    }

  }
}
