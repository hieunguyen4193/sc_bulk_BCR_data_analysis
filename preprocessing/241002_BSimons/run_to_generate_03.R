gc()
rm(list = ls())

if (packageVersion("Matrix") != "1.5.4.1"){
  install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", type = "source", repos = NULL)
  install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.3.tar.gz", type = "source", repos = NULL)
}
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241002_BSimons"
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "241002_BSimons"

path.to.rmd1 <- file.path(path.to.project.src, "03_downstream_analysis_singlet_only.Rmd")


path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.save.html <- file.path(path.to.main.output, "html_output")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

all.samples <- c("PP3")
all.quantiles <- c(0.85, 0.99, 0.95, 0.90)

rerun <- TRUE
for (chosen.quantile in all.quantiles){
  for (sample.id in all.samples){
    print("--------------------------------------------------------------------")
    print(sprintf("%s --  %s", chosen.quantile, sample.id))
    print("--------------------------------------------------------------------")
    
    if (file.exists(file.path(path.to.save.html, 
                              sprintf("03_downstream_analysis_singlet_only_quantile_%s_sample_%s.html", 
                                      chosen.quantile, 
                                      sample.id))) == FALSE | rerun == TRUE){
      rmarkdown::render(input = path.to.rmd1,
                        params = list(
                          outdir = outdir,
                          PROJECT = PROJECT,
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
