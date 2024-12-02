gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.04.output <- file.path(outdir, "VDJ_output", "04_output")
path.to.09.output <- file.path(outdir, "VDJ_output", "09_output")
dir.create(path.to.09.output, showWarnings = FALSE, recursive = TRUE)

sc.datasets <- c("240805_BSimons")
bulk.dataset <- c("240826_BSimons")

mouse.id <- "m1"

full.clonedf <- read.csv(file.path(path.to.04.output, "full_clonedf_with_mutation_rate.csv"))
sc.clonedf <- subset(full.clonedf, full.clonedf$dataset.name %in% sc.datasets)[, c("barcode", "id", 'V.gene', 'J.gene', 'D.gene', "aaSeqCDR3", "nSeqCDR3")]
sc.clonedf <- subset(sc.clonedf, sc.clonedf$id %in% c("M1", "P1"))
tree.seqdf.dir <- list(
  `241031_BSimons` = "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1/tree_analysis/07_output/241031_BSimons_240411_BSimons_241002_BSimons",
  `240826_BSimons` = "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1/tree_analysis/06_output/240826_BSimons_240805_BSimons"
)

all.tree.files <- Sys.glob(file.path(tree.seqdf.dir, "*"))

names(all.tree.files) <- to_vec(
  for (item in all.tree.files){
    item <- basename(item)
    item <- str_replace(item, "tree_Seqdf_", "")
    item <- str_replace(item, ".csv", "")
  }
)

input.file <- all.tree.files[["m1_IGHV1-47-01_IGHJ3-01_39_1"]]

filename <- basename(input.file)
input.V.gene <- str_split(filename, "_")[[1]][[4]]
input.J.gene <- str_split(filename, "_")[[1]][[5]]
cdr3_len <- str_split(filename, "_")[[1]][[6]] %>% as.numeric()
bulkdf <- read.csv(input.file)

tmpdf <- subset(sc.clonedf, 
       sc.clonedf$V.gene == str_split(input.V.gene, "-0")[[1]][[1]] & 
       sc.clonedf$J.gene == str_split(input.J.gene, "-0")[[1]][[1]] &
       nchar(sc.clonedf$nSeqCDR3) == cdr3_len ) %>%
  rowwise() %>%
  mutate(barcode = sprintf("%s_%s", id, barcode))

for (input.barcode in unique(tmpdf$barcode)){
  bulkdf[[input.barcode]] <- unlist(
    lapply(
      bulkdf$aaSeqCDR3, function(x){
        input.seq <- subset(tmpdf, tmpdf$barcode == input.barcode)$aaSeqCDR3
        d <- stringdist(input.seq, x, method = "lv")      
        d <- d/max(nchar(x), nchar(input.seq))
        return(d)
      }
    )
  )
}

plotdf <- bulkdf[, c("seqid", unique(tmpdf$barcode))]
plotdf.pivot <- plotdf %>% pivot_longer(!seqid, names_to = "barcode", values_to = "dist") 
plotdf.pivot %>% ggplot(aes(x = barcode, y = seqid, fill = dist)) + geom_tile(color = "white") +
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_fill_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")

