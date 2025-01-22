path.to.00.output <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1/GEX_output/00_output"
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"
path.to.all.s.obj <- list(
  `1st_dataset` = file.path(path.to.00.output, "1st_dataset", "1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1.addedCloneInfo.rds"),
  `2nd_dataset` = file.path(path.to.00.output, "2nd_dataset", "2nd_dataset_removed_5_6.without_reInt.res1.addedCloneInfo.rds"),
  `240805_BSimons` = file.path(path.to.00.output, "240805_BSimons", "240805_BSimons.output.s8.addedCloneInfo.rds"),
  `241002_BSimons` = file.path(path.to.00.output, "241002_BSimons", "GEX_sample_PP3_seurat_object.addedCloneInfo.rds"),
  `241104_BSimons` = file.path(path.to.00.output, "241104_BSimons", "GEX_sample_PP7_seurat_object.addedCloneInfo.rds"),
  BonnData = file.path( path.to.00.output, "BonnData", "GEX_sample_BonnData_seurat_object.addedSampleID.addedCloneInfo.rds"),
  Dataset1_2 = file.path(outdir, "GEX_output/04_output", "Dataset1_2", sprintf("s8_output/%s.output.s8.rds", "Dataset1_2")),
  Dataset1_2_Bonn = file.path(outdir, "GEX_output/04_output", "Dataset1_2_Bonn", sprintf("s8_output/%s.output.s8.rds", "Dataset1_2_Bonn")),
  Dataset1_Bonn = file.path(outdir, "GEX_output/04_output", "Dataset1_Bonn", sprintf("s8_output/%s.output.s8.rds", "Dataset1_Bonn")),
  Dataset2_Bonn = file.path(outdir, "GEX_output/04_output", "Dataset2_Bonn", sprintf("s8_output/%s.output.s8.rds", "Dataset2_Bonn")),
  `240805_BSimons_filterHT` = file.path(path.to.00.output, "240805_BSimons_filterHT/240805_BSimons.output.s8.addedCloneInfo.rds"),
  `240805_BSimons_filterHT_cluster` = file.path(path.to.00.output, "240805_BSimons_filterHT_cluster/240805_BSimons.output.s8.addedCloneInfo.rds"),
  `240805_BSimons_filterHT_cluster_renamed` = file.path(path.to.00.output, "240805_BSimons_filterHT_cluster_renamed/240805_BSimons.renamedClusters.output.s8.addedCloneInfo.rds"),
  `241002_241104_BSimons` = file.path(path.to.00.output, "241002_241104_BSimons", "241002_241104_BSimons.output.s8.addedCloneInfo.rds")
)

