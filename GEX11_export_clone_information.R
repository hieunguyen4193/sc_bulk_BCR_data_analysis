gc()
rm(list = ls())
#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

library(viridis)
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}
if ("ggthemes" %in% installed.packages() == FALSE){
  install.packages("ggthemes")
}

#####---------------------------------------------------------------------------#####
##### INPUT ARGS
#####---------------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNSD01/storage/all_BSimons_datasets"
path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_path_to_output.R"))
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))

outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

path.to.11.output <- file.path(outdir, "GEX_output", "11_output")
dir.create(path.to.11.output, showWarnings = FALSE, recursive = TRUE)

clonedf <- list()
for (input.dataset in names(path.to.all.s.obj)){
  if (file.exists(file.path(path.to.11.output, sprintf("%s.csv", input.dataset))) == FALSE){
    s.obj <- readRDS(path.to.all.s.obj[[input.dataset]])
    meta.data <- s.obj@meta.data %>% subset(select = c(
      CTaa,
      CTstrict,
      V.gene, 
      J.gene, 
      aaSeqCDR3, 
      nSeqCDR3, 
      VJseq.combi,
      VJ.combi,
      VJ.len.combi,        
      VJcombi_CDR3_0.85 
    )) %>% 
      subset(is.na(VJcombi_CDR3_0.85) == FALSE) %>% 
      rownames_to_column("barcode")
    write.csv(meta.data, file.path(path.to.11.output, sprintf("%s.csv", input.dataset)))
    clonedf[[input.dataset]] <- meta.data
  } else {
    print(sprintf("Clone information for dataset %s exists, reading in ...", input.dataset))
    clonedf[[input.dataset]] <- read.csv(file.path(path.to.11.output, sprintf("%s.csv", input.dataset))) %>%
      subset(select = -c(X))
  }
} 

if (file.exists(file.path(path.to.11.output, "check_SPEC_IGH_clone_definition.xlsx")) == FALSE){
  summarydf <- data.frame()
  
  for (input.dataset in names(clonedf)){
    df <- clonedf[[input.dataset]]
    countdf <- data.frame(clone = unique(df$VJcombi_CDR3_0.85))
    countdf <- countdf %>% rowwise() %>%
      mutate(num.CTstrict = subset(df, df$VJcombi_CDR3_0.85 == clone)$CTstrict %>% unique() %>% length()) %>%
      mutate(CTstrict_list = paste(subset(df, df$VJcombi_CDR3_0.85 == clone)$CTstrict %>% unique(), collapse = ";")) %>%
      arrange(desc(num.CTstrict))
    
    num1 <- subset(countdf, countdf$num.CTstrict == 1) %>% nrow()
    pct1 <- num1/nrow(countdf)
    
    tmp.summarydf <- data.frame(dataset = input.dataset,
                                numClone = nrow(countdf),
                                numCTstrict = length(unique(df$CTstrict)),
                                numSPEC = num1, 
                                pctSPEC = pct1)
    summarydf <- rbind(summarydf, tmp.summarydf)
  }
  writexl::write_xlsx(summarydf, file.path(path.to.11.output, "check_SPEC_IGH_clone_definition.xlsx"))  
} else {
  print("File exists, reading in ...")
  summarydf <- readxl::read_excel(file.path(path.to.11.output, "check_SPEC_IGH_clone_definition.xlsx"))
}


