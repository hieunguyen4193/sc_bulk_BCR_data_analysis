gc()
rm(list = ls())

path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/240805_BSimons"
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "240805_BSimons"

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.save.html <- file.path(path.to.main.output, "html_output", "05_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.samples <- c("M1", "M2", "M3", "P1", "P2", "P3")
all.quantiles <- c(0.99, 0.95, 0.90, 0.85)
all.integration.case <- list(
  all_samples = all.samples,
  mouse1 = c("M1", "P1"),
  mouse2 = c("M2", "P2"),
  mouse3 = c("M3", "P3")
)

path.to.rmd <- file.path(path.to.project.src, "05_integrated_data_analysis.Rmd")

for (chosen.quantile in all.quantiles){
  for (integration.case in names(all.integration.case)){
    print(sprintf("Working on integration case %s and quantile %s", integration.case, chosen.quantile))
    rmarkdown::render(input = path.to.rmd,
                      params = list(
                        outdir = outdir, 
                        PROJECT = PROJECT,
                        integration.case = integration.case,
                        chosen.quantile = chosen.quantile
                      ),
                      output_file = sprintf("05_integrated_data_analysis_quantile_%s_integration_%s.html", chosen.quantile, integration.case),
                      output_dir = path.to.save.html)
  }
}

