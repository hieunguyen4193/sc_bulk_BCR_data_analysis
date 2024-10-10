gc()
rm(list = ls())
maindir <- "/home/hieunguyen/CRC1382/outdir/molmed_server"
data.version <- "Fabio_data"
data.type <- "mouse based trees"
path.to.main.input <- file.path(maindir, data.version, data.type)
all.mouses <- Sys.glob(file.path(path.to.main.input, "m*_full"))

mouse.id <- "m11"
path.to.mouse.tree.data <- file.path(path.to.main.input, sprintf("%s_full", mouse.id))

path.to.mixcr.output <- file.path(maindir, data.version, "analysis fabio without trees/full")

all.mid.tables <- Sys.glob(file.path(path.to.mixcr.output, "*")) 

tmpdf <- read.table(all.mid.tables[[1]], sep = "\t")
tmpdf
