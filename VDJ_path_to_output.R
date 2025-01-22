outdir <- "/media/hieunguyen/GSHD_HN01/outdir"
path.to.VDJ.output <- file.path(outdir, "sc_bulk_BCR_data_analysis_v0.1", "VDJ_output")

all.integration.cases <- list(
  Dataset1_2 = c("1st_dataset", "2nd_dataset"),
  Dataset1_2_Bonn = c("1st_dataset", "2nd_dataset", "BonnData"),
  Dataset1_Bonn = c("1st_dataset", "BonnData"),
  Dataset2_Bonn = c("2nd_dataset", "BonnData")
)

path.to.all.VDJ.output <- list(
  `2nd_dataset` = file.path(path.to.VDJ.output, "1st_2nd_BSimons_Datasets", "VDJ_output_0.85", "1st_2nd_BSimons_Datasets_combined_contigs.rds"),
  `1st_dataset` = file.path(path.to.VDJ.output, "1st_2nd_BSimons_Datasets", "VDJ_output_0.85", "1st_2nd_BSimons_Datasets_combined_contigs.rds"),
  `240805_BSimons` = file.path(path.to.VDJ.output, "240805_BSimons", "VDJ_output_0.85", "240805_BSimons_combined_contigs.rds"),
  `241104_BSimons` = file.path(path.to.VDJ.output, "241104_BSimons", "VDJ_output_0.85", "241104_BSimons_combined_contigs.rds"),
  `241002_BSimons` = file.path(path.to.VDJ.output, "241002_BSimons", "VDJ_output_0.85", "241002_BSimons_combined_contigs.rds"),
  BonnData = file.path(path.to.VDJ.output, "BonnData", "VDJ_output_0.85", "BonnData_combined_contigs.rds"),
  `240805_BSimons_filterHT` = file.path(path.to.VDJ.output, "240805_BSimons", "VDJ_output_0.85", "240805_BSimons_combined_contigs.rds"),
  `240805_BSimons_filterHT_cluster` = file.path(path.to.VDJ.output, "240805_BSimons", "VDJ_output_0.85", "240805_BSimons_combined_contigs.rds"),
  `240805_BSimons_filterHT_cluster_renamed` = file.path(path.to.VDJ.output, "240805_BSimons", "VDJ_output_0.85", "240805_BSimons_combined_contigs.rds"),
  `241002_241104_BSimons` = file.path(path.to.VDJ.output, "241002_241104_BSimons", "VDJ_output_0.85", "241002_241104_BSimons_combined_contigs.rds")
)

dataset1_sample_list <- c("17_MM9_Ecoli",
                          "20_MM9_Ecoli",
                          "21_MM9_Ecoli_SPF",
                          "MM9_S2",
                          "MM9_S4",          
                          "MM9_SPF_S3",
                          "MM9_SPF_S9")

dataset2_sample_list <- c( "Sample_132",
                           "Sample_133")

##### process VDJ merge output for integrated dataset.
##### 241002_241104_BSimons
if (file.exists(path.to.all.VDJ.output[["241002_241104_BSimons"]]) == FALSE){
  dir.create(dirname(path.to.all.VDJ.output[["241002_241104_BSimons"]]), 
             showWarnings = FALSE, recursive = TRUE)
  vdj.output1 <- readRDS(path.to.all.VDJ.output[["241002_BSimons"]])
  vdj.output2 <- readRDS(path.to.all.VDJ.output[["241104_BSimons"]])
  
  vdj.output <- c(vdj.output1, vdj.output2)
  saveRDS(vdj.output, path.to.all.VDJ.output[["241002_241104_BSimons"]])  
}


