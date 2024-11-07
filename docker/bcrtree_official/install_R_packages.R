################################################################################
# This script is used to clean the envrionment and import all necessary packages
################################################################################
# Specify the list of packages that need to be imported ########################
list.of.packages <- c("BiocManager",
                      "optparse", 
                      "tidyverse", 
                      "ggplot2", 
                      "SoupX",
                      "comprehenr",
                      "DoubletFinder",
                      "vroom",
                      "hash",
                      "DT",
                      "janitor",
                      "knitr",
                      "circlize",
                      "formattable",
                      "htmlwidgets",
                      "plotly",
                      "stringr",
                      "rcompanion",
                      "argparse",
                      "scatterpie", 
                      "scales",
                      "rstatix",
                      "remotes",
                      "PPCI",
                      "writexl",
                      "devtools",
                      "svglite",
                      "lme4",
                      "ggpubr", 
                      "testit", 
                      "heatmaply",
                      "devtools",
                      "igraph"
)

bioc.packages <- c("Seurat",
                   "SingleCellExperiment", 
                   "celda", 
                   "BiocSingular", 
                   "PCAtools", 
                   "SingleCellExperiment",
                   "scRepertoire", 
                   "sctransform", 
                   "progeny",
                   "powerTCR",
                   "scRepertoire",
                   "org.Hs.eg.db",
                   "org.Mm.eg.db",
                   "DESeq2",
                   "diptest",
                   "CATALYST",
                   "HDCytoData",
                   "diffcyt",
                   "FlowSOM",
                   "flowCore", 
                   "scater",
                   "enrichplot",
                   "msa")

# Check if packages are installed ##############################################

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages #########################################################
if(length(new.packages)) install.packages(new.packages)
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages, update = FALSE, ask = TRUE)

# Import all packages ##########################################################
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)

##### install ClusterProfiler
BiocManager::install("msigdbr", update = FALSE)
BiocManager::install("org.Hs.eg.db", update = FALSE)
for (pkg in c("clusterProfiler", "DOSE", "GOSemSim")){
    if (pkg %in% installed.packages() == TRUE){
        remove.packages(pkg)
    }
}
devtools::install_github("YuLab-SMU/yulab.utils", upgrade = "never")
if ("clusterProfiler" %in% installed.packages() == TRUE){
  if (packageVersion("clusterProfiler") != "4.13.4"){
    remove.packages("DOSE")
    remove.packages("GOSemSim")
    remove.packages("yulab.utils")
    remove.packages("clusterProfiler")
    remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
    remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")
  }
} else {
  remotes::install_github("GuangchuangYu/GOSemSim", upgrade = "never")
  remotes::install_github("GuangchuangYu/clusterProfiler", upgrade = "never")
}

new.pkgs <- c("svglite", "Matrix", "car")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() == FALSE){
    install.packages(pkg)
  }   
}

if (packageVersion("igraph") != "2.1.1"){
  install.packages("https://cran.r-project.org/src/contrib/igraph_2.1.1.tar.gz", type = "source", repos = NULL)
}
if (packageVersion("ggplot2") != "3.4.4"){
  install.packages("https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.4.4.tar.gz", type = "source", repos = NULL)
}

if ("msa" %in% installed.packages() == FALSE){
  BiocManager::install("msa", update = FALSE)
}

if (packageVersion("tidyr") != "1.3.1"){
  install.packages("tidyr")  
  install.packages("tidyselect")
  library(tidyr)
}
# EOF ##########################################################################
