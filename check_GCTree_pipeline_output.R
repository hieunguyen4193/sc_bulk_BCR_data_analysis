gc()
rm(list = ls())

checkdf <- list()

# input.dataset <- "240826_BSimons"
for (input.dataset in c("240826_BSimons", "241031_BSimons", "220701_etc_biopsies")){
  checkdf[[input.dataset]] <- list()
  path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
  path.to.gctree.output <- sprintf("/media/hieunguyen/HNHD01/storage/all_BSimons_datasets/%s/GCtrees/v0.2", input.dataset)
  path.to.samplesheet <- sprintf("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/SampleSheet_GCTree_nextflow_%s", input.dataset)
  
  all.samplesheets <- Sys.glob(file.path(path.to.samplesheet, "*.csv"))
  names(all.samplesheets) <- to_vec(
    for (item in all.samplesheets){
      item <- basename(item)
      item <- str_replace(item, sprintf("SampleSheet_GCTree_%s_", input.dataset), "")
      item <- str_replace(item, ".nextflow.csv", "")
      # item <- str_split(item, "_")[[1]][2:length(str_split(item, "_")[[1]])] %>% paste(collapse = "_")
    }
  )
  for (input.name in names(all.samplesheets)){
    print(sprintf("working on dataset %s and sample %s", input.dataset, input.name))
    input.file <- all.samplesheets[[input.name]]
    sampledf <- read.csv(input.file)
    sample.id <- str_split(input.name, "_")[[1]][[1]]
    sample.condition <- str_replace(input.name, sprintf("%s_", sample.id), "")
    if (input.dataset == "240826_BSimons"){
      all.trees <- Sys.glob(file.path(path.to.gctree.output, sample.condition, input.name, "*"))      
    } else {
      all.trees <- Sys.glob(file.path(path.to.gctree.output, "all", input.name, "*"))
    }

    names(all.trees) <- basename(all.trees)
    treedf <- data.frame(filename = names(all.trees), path_to_output = all.trees)
    
    sampledf <- merge(sampledf, treedf, by.x = "filename", by.y = "filename", all.x = TRUE)
    sampledf <- sampledf %>% rowwise() %>%
      mutate(check = ifelse(
        length(Sys.glob(file.path(path_to_output, "02*", "*"))) != 0, "yes", "no"
      ))
    if (nrow(sampledf) != length(all.trees)){
      print(input.name)
    }
    checkdf[[input.dataset]][[input.name]] <- subset(sampledf, sampledf$check == "no")
  }
}

