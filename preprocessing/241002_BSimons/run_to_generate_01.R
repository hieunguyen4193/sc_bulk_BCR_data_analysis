gc()
rm(list = ls())

if (packageVersion("Matrix") != "1.5.4.1"){
  install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", type = "source", repos = NULL)
  install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.3.tar.gz", type = "source", repos = NULL)
  devtools::install_github("ncborcherding/scRepertoire@v1", upgrade = FALSE)
}

path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241002_BSimons"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis"
PROJECT <- "241002_BSimons"

path.to.rmd1 <- file.path(path.to.project.src, "01_hashtag_antibodies.Rmd")
path.to.rmd2 <- file.path(path.to.project.src, "01_comparing_different_quantile.Rmd")

path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.save.html <- file.path(path.to.main.output, "html_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.samples <- c("PP3")
all.quantiles <- c(0.99, 0.95, 0.90, 0.85)
for (sample.id in all.samples){
  for (chosen.quantile in all.quantiles){
    print(sprintf("Working on sample %s, quantile %s", sample.id, chosen.quantile))
    if (file.exists(file.path(path.to.save.html, sprintf("01_hashtag_antibodies_analysis_%s.html", sample.id))) == FALSE){
      rmarkdown::render(input = path.to.rmd1,
                        params = list(
                          outdir = outdir, 
                          PROJECT = PROJECT,
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
                        sample.id = sample.id
                      ),
                      output_file = sprintf("01_comparing_different_quantile_%s.html", sample.id),
                      output_dir = path.to.save.html)
  } else {
    print(sprintf("html output for sample %s existed", sample.id))
  }
}