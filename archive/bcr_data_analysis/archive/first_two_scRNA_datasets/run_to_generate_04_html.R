gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### MAIN
#####----------------------------------------------------------------------#####
path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/BSimons/CRC1382_BSimons_project/BSimons_scRNAseq_analysis"

source(file.path(path.to.project.src, "00_helper_functions.R"))

outdir <- "/media/hieunguyen/HD0/outdir/CRC1382"

path.to.rmd <- file.path(path.to.project.src, "04_trajectory_analysis_with_monocle2.Rmd")

path.to.save.html <- file.path(outdir, "/BSimons/OUTPUT/THESIS_OUTPUT_20231026/html_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)

candidate.root.node <- NULL

for (dataset.name in c("1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1",
                       "2nd_dataset_removed_5_6.without_reInt.res1")){
  save.html.name <- sprintf("%s.monocle2.html", dataset.name)
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    rmarkdown::render(input = path.to.rmd, 
                      output_file = save.html.name, 
                      output_dir = path.to.save.html,
                      params = list(
                        dataset.name = dataset.name
                      ))  
  }  
}


