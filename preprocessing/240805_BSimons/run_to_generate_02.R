gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240805_BSimons"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "240805_BSimons"

path.to.rmd1 <- file.path(path.to.project.src, "02_downstream_analysis_all_cells.Rmd")

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.save.html <- file.path(path.to.main.output, "html_output", "02_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
for (sample.id in all.samples){
  if (file.exists(file.path(path.to.save.html, sprintf("02_downstream_analysis_%s.html", sample.id))) == FALSE){
    print(sprintf("Working on sample %s", sample.id))
    rmarkdown::render(input = path.to.rmd1,
                      params = list(
                        outdir = outdir, 
                        PROJECT = PROJECT,
                        sample.id = sample.id
                      ),
                      output_file = sprintf("02_downstream_analysis_%s.html", sample.id),
                      output_dir = path.to.save.html)
  } else {
    print(sprintf("html output for sample %s existed", sample.id))
  }
}