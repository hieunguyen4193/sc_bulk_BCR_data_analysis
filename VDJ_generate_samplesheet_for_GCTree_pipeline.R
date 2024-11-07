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

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####

path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
thres <- 0.85
thres.dis <- 0.15
PROJECT <- "240826_BSimons"
ref.gene <- "IMGT"
ref.gene.config <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/ref_gene_config.R"

#####----------------------------------------------------------------------#####
##### PATHS
#####----------------------------------------------------------------------#####
path.to.save.fasta <- file.path(outdir, "VDJ_output", "03_output", "FASTA_output", PROJECT,  sprintf("VDJ_output_%s", thres))
dir.create(path.to.save.fasta, showWarnings = FALSE, recursive = TRUE)

all.files <- Sys.glob(file.path(path.to.save.fasta, "*/*/*.fasta"))

planned.outdir <- "/home/hieu/outdir"
countdf <- data.frame()
for (input.file in all.files){
  filename <- str_replace(basename(input.file), ".fasta", "")
  input.case <- dirname(input.file) %>% basename()
  mouse.id <- dirname(dirname(input.file)) %>% basename()
  c <-readDNAStringSet(input.file) %>% as.data.frame() %>% nrow()
  path.to.file <- file.path(planned.outdir, 
                            "sc_bulk_BCR_data_analysis_v0.1",
                            "VDJ_output", 
                            "03_output",
                            "FASTA_output", 
                            "PROJECT",
                            sprintf("VDJ_output_%s", thres),
                            mouse.id, 
                            input.case,
                            filename
                            )
  tmpdf <- data.frame(
    filename = c(filename),
    path = c(path.to.file),
    mouse.id = c(mouse.id),
    input.case = c(input.case),
    count = c(c)
  )
  countdf <- rbind(countdf, tmpdf)
}

countdf <- countdf %>% arrange(count)
write.csv(countdf, file.path(path.to.main.src, sprintf("SampleSheet_GCTree_%s.csv", PROJECT)))
write.csv(countdf[, c("filename", "path")], file.path(path.to.main.src, sprintf("SampleSheet_GCTree_%s.nextflow.csv", PROJECT)))
