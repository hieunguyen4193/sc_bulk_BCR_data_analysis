#####----------------------------------------------------------------------#####
##### INSTALL MONOCLE2
#####----------------------------------------------------------------------#####
remove.packages("monocle3")
devtools::install_github("cysouw/qlcMatrix")
install.packages("DDRTree")
install.packages("densityClust")
BiocManager::install("monocle.objSingleCell", update = FALSE)
install.packages("fastICA")
BiocManager::install("biocViews", update = FALSE)
# remove.packages("BiocGenerics")
BiocManager::install("HSMMSingleCell", update = FALSE)
# install.packages("/media/hieunguyen/HD0/storage/offline_pkgs/monocle_2.30.0.tar.gz", type = "source", repos = NULL)
install.packages("/media/hieunguyen/HD0/storage/offline_pkgs/monocle_2.30.0/monocle", type = "source", repos = NULL)
BiocManager::install("tradeSeq", update = FALSE)