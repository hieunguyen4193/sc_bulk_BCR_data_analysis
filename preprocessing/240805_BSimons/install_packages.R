gc()
rm(list = ls())

install.packages("https://cran.r-project.org/src/contrib/Archive/reticulate/reticulate_1.32.0.tar.gz", type = "source", repos = NULL)
library(reticulate)

reticulate::install_python(version = '3.9')

py_config() # run this, yes to create miniconda env and restart R session.

install.packages("keras")
use_virtualenv(virtualenv = "/home/rstudio/.virtualenvs/r-reticulate/", required = TRUE)
library(keras)
install_keras()

# #### install scRepertoire 2.0
# if (packageVersion("scRepertoire") != "2.0.0"){
#   install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-5.tar.gz", type = "source", repos = NULL)
#   install.packages(c("ggdendro", "iNEXT", "quantreg"))
#   install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/scRepertoire_2.0.0.tar.gz", type = "source", repos = NULL)
# }

##### install the legacy version v1 of scRepertoire package
devtools::install_github("ncborcherding/scRepertoire@v1", upgrade = FALSE)

devtools::install_github("ncborcherding/Ibex", upgrade = FALSE)
devtools::install_github("ncborcherding/Trex", upgrade = FALSE)

library("Trex")
library("Ibex")

# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", repos = NULL, type = "source")
# install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", repos = NULL, type = "source")
# devtools::install_github("ncborcherding/scRepertoire@v1")
