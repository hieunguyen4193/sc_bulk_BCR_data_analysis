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
PROJECT <- "220701_etc_biopsies"
ref.gene <- "IMGT"
ref.gene.config <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/ref_gene_config.R"


source(ref.gene.config) # path to the configuration file for the input BCR reference genes
if (ref.gene == "10x"){
  ref.fasta <- readDNAStringSet(ref.genes$`10x`)
} else if (ref.gene == "IMGT"){
  s.V.genes <- readDNAStringSet(ref.genes$IMGT$V.gene)
  names(s.V.genes) <- lapply(names(s.V.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    # x <- str_split(x, "[*]")[[1]][[1]]
    return(x)
  })
  s.J.genes <- readDNAStringSet(ref.genes$IMGT$J.gene)
  names(s.J.genes) <- lapply(names(s.J.genes), function(x){
    x <- str_split(x, "[|]")[[1]][[2]]
    # x <- str_split(x, "[*]")[[1]][[1]]
    return(x)
  })
}


tmpdf <- data.frame(org.gene = names(s.V.genes))
