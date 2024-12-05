gc()
rm(list = ls())

my_random_seed <- 42
set.seed(my_random_seed)
config.version <- "default"
# __________VDJ DATA ANYLYSIS PIPELINE__________

if (packageVersion("Matrix") != "1.5.4.1"){
  install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", type =  "source", repos = NULL)  
}
if (packageVersion("scRepertoire") != "1.11.0"){
  devtools::install_github("ncborcherding/scRepertoire@v1", upgrade = FALSE)
}

source("/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")

path.to.storage <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"

outdir <- "/media/hieunguyen/HNSD_mini/outdir/sc_bulk_BCR_data_analysis_v0.1"
PROJECT <- "BonnData"
thres <- 0.85

path.to.main.input <- file.path(path.to.storage, PROJECT)

path.to.main.output <- file.path(outdir, PROJECT)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

path.to.VDJ.input <- file.path(path.to.main.input, "VDJ")
path.to.VDJ.output <- file.path(outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres))
dir.create(path.to.VDJ.output, showWarnings = FALSE, recursive = TRUE)

path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source(file.path(path2src, "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

summarize_vdj_data(path.to.input = path.to.VDJ.input, 
                   path.to.output = path.to.VDJ.output, 
                   PROJECT = PROJECT, 
                   removeNA = FALSE, 
                   removeMulti = FALSE, 
                   T_or_B = "B",
                   threshold = thres)

path2input <- file.path(path.to.main.input, "GEX")

path.to.project.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/BonnData"
source(file.path(path.to.project.src, "config.R"))
analysis.round <- "1st"
# _____stage lst for single sample_____
stage_lst <- c(BonnData = "BonnData")

MINCELLS  = 5
MINGENES  = 50

save.RDS <- list(s1 = TRUE,
                 s2 = TRUE,
                 s3 = TRUE,
                 s4 = TRUE,
                 s5 = TRUE,
                 s6 = TRUE,
                 s7 = TRUE,
                 s8 = TRUE,
                 s8a = TRUE,
                 s9 = FALSE)

sw <- list(s1 = "on",
           s2 = "on",
           s3 = "on",
           s4 = "on",
           s5 = "on",
           s6 = "on",
           s7 = "on",
           s8 = "off",
           s8a = "on",
           s9 = "off")

rerun <- list(s1 = FALSE, 
              s2 = FALSE,
              s3 = FALSE,
              s4 = FALSE,
              s5 = FALSE,
              s6 = FALSE,
              s7 = FALSE,
              s8 = FALSE,
              s8a = FALSE,
              s9 = FALSE)

filter.thresholds <- filter.config.params[[config.version]]

remove_doublet <- FALSE
path.to.10X.doublet.estimation <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/DoubletEstimation10X.csv"

#####--------------------------------------------------------------------#####
##### IMPORTANT INPUT PARAMS
#####--------------------------------------------------------------------#####
for (sample in names(stage_lst)){
  print(sprintf("Working on sample %s", sample))
  path.to.output <- file.path(path.to.main.output, sprintf("%s_round", analysis.round))
  dir.create(path.to.output, showWarnings = FALSE)
  path.to.anno.contigs <- file.path(path.to.VDJ.output, sprintf("annotated_contigs_clonaltype_%s.csv", sample))
  path.to.count.clonaltype <- file.path(path.to.VDJ.output, sprintf("count_clonaltype_%s.csv", sample))
  filtered.barcodes <- NULL
  path.to.s3a <- NULL
  input.method <- "filterIG_with_hashtagAntibody"
  s.obj <- run_pipeline_GEX(path2src=path2src,
                            path2input=file.path(path2input, sample),
                            path.to.logfile.dir=file.path(path.to.output, sprintf("%s_%s_round", sample, analysis.round), "logs"),
                            stage_lst=stage_lst[[sample]],
                            path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                            MINCELLS=MINCELLS,
                            MINGENES=MINGENES,
                            PROJECT=PROJECT,
                            remove_doublet=remove_doublet,
                            save.RDS=save.RDS,
                            path.to.output=file.path(path.to.output, sprintf("%s_%s_round", sample, analysis.round)),
                            rerun=rerun,
                            DE.test="wilcox",
                            num.PCA=num.PCA,
                            num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                            use.sctransform=use.sctransform,
                            filtered.barcodes=filtered.barcodes,
                            filter.thresholds=filter.thresholds,
                            path.to.anno.contigs=path.to.anno.contigs,
                            path.to.count.clonaltype=path.to.count.clonaltype,
                            input.method = input.method,
                            my_random_seed = my_random_seed,
                            with.VDJ = TRUE, 
                            path.to.s3a.source = path.to.s3a, 
                            regressOut_mode = regressOut_mode,
                            features_to_regressOut = features_to_regressOut,
                            sw = sw,
                            cluster.resolution = cluster.resolution)
}
#### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))