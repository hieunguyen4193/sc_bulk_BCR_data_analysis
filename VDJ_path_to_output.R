outdir <- "/media/hieunguyen/GSHD_HN01/outdir"
path.to.VDJ.output <- file.path(outdir, "sc_bulk_BCR_data_analysis_v0.1", "VDJ_output")
path.to.all.VDJ.output <- list(
  `2nd_dataset` = file.path(path.to.VDJ.output, "1st_2nd_BSimons_Datasets", "VDJ_output_0.85", "1st_2nd_BSimons_Datasets_combined_contigs.rds"),
  `1st_dataset` = file.path(path.to.VDJ.output, "1st_2nd_BSimons_Datasets", "VDJ_output_0.85", "1st_2nd_BSimons_Datasets_combined_contigs.rds"),
  `240805_BSimons` = file.path(path.to.VDJ.output, "240805_BSimons", "VDJ_output_0.85", "240805_BSimons_combined_contigs.rds"),
  `241104_BSimons` = file.path(path.to.VDJ.output, "241104_BSimons", "VDJ_output_0.85", "241104_BSimons_combined_contigs.rds"),
  `241002_BSimons` = file.path(path.to.VDJ.output, "241002_BSimons", "VDJ_output_0.85", "241002_BSimons_combined_contigs.rds")
)

