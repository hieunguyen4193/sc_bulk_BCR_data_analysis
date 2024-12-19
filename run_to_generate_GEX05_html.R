gc()
rm(list = ls())

path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))

scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))

all.integration.cases <- list(
  Dataset1_2 = c("1st_dataset", "2nd_dataset"),
  Dataset1_2_Bonn = c("1st_dataset", "2nd_dataset", "BonnData"),
  Dataset1_Bonn = c("1st_dataset", "BonnData"),
  Dataset2_Bonn = c("2nd_dataset", "BonnData")
)
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.rmd <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/GEX05_data_analysis_for_integrated_datasets.Rmd"
path.to.save.html <- file.path(outdir, "GEX_output", "05_output", "HTML_outputs")
dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
for (integration.case in names(all.integration.cases)){
  save.html.name <- str_replace(basename(path.to.rmd), ".Rmd", sprintf(".%s.html", integration.case))
  if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
    print(sprintf("Generate html file for the case %s", integration.case))
    rmarkdown::render(input = path.to.rmd, 
                      output_file = save.html.name, 
                      output_dir = path.to.save.html, 
                      params = list(integration.case = integration.case))
  } else {
    print(sprintf("HTML output for %s exists", integration.case))
  }
}