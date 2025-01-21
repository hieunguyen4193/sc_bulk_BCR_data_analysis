gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"

path.to.rmd <- file.path(path.to.main.src, "preprocessing", "241002_241104_BSimons", "02_integrated_data_analysis.Rmd")

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "241002_241104_BSimons"
path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.01.output <- file.path(path.to.main.output, "01_output")
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

path.to.save.html <- file.path(path.to.main.output, "html_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

rmarkdown::render(input = path.to.rmd, 
                  output_file = sprintf("02_integrated_data_analysis_%s.Rmd", PROJECT),
                  output_dir = path.to.save.html)