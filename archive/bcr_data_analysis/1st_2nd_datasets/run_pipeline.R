gc()
rm(list = ls())

# we must re-install the dev version from github/scRepertoire to use 
# the function combineBCR with threshold option. 
# remove.packages("scRepertoire")
# install.packages("devtools")
# devtools::install_github("ncborcherding/scRepertoire")

# Run all samples from the first and second datasets, each samples is processed
# separately.

# install.packages("https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.5-4.1.tar.gz", type =  "source", repos = NULL)
# devtools::install_github("ncborcherding/scRepertoire@v1")

# __________VDJ DATA ANYLYSIS PIPELINE__________
path.to.storage <- "/media/hieunguyen/GSHD_HN01/storage"
outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "1st_2nd_BSimons_Datasets"

path.to.main.input <- file.path(path.to.storage, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)

source("/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")
path.to.VDJ.input <- file.path(path.to.main.input, "VDJ")
path.to.VDJ.output <- file.path(path.to.main.output, "VDJ_output")

dir.create(path.to.VDJ.output, showWarnings = FALSE)

clone_similarity_threshold <- 0.95

summarize_vdj_data(path.to.VDJ.input, 
                   path.to.VDJ.output, 
                   PROJECT, 
                   removeNA=FALSE, 
                   removeMulti=FALSE, 
                   T_or_B = "B",
                   threshold = clone_similarity_threshold)

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.pipeline.src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

set.seed(411)

source(file.path(path2src, "import_libraries.R"))

source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))

for (analysis.round in c("1st")){
  
  path2input <- file.path(path.to.main.input, "GEX")
  
  stage_lst <- hash()
  stage_lst[["Sample_132"]] <- c(Sample_132 = "Sample_132")
  stage_lst[["Sample_133"]] <- c(Sample_133 = "Sample_133")
  stage_lst[["17_MM9_Ecoli"]] <- c(`17_MM9_Ecoli` = "MM9_Ecoli")
  stage_lst[["20_MM9_Ecoli"]] <- c(`20_MM9_Ecoli` = "MM9_Ecoli")
  stage_lst[["21_MM9_Ecoli_SPF"]] <- c(`21_MM9_Ecoli_SPF` = "MM9_Ecoli_SPF")
  stage_lst[["MM9_S2"]] <- c(MM9_S2 = "MM9")
  stage_lst[["MM9_S4"]] <- c(MM9_S4 = "MM9")
  stage_lst[["MM9_SPF_S3"]] <- c(MM9_SPF_S3 = "MM9_SPF")
  stage_lst[["MM9_SPF_S9"]] <- c(MM9_SPF_S9 = "MM9_SPF")
  
  
  MINCELLS  = 5
  MINGENES  = 50
  
  save.RDS <- list(s1 = TRUE,
                   s2 = TRUE,
                   s3 = TRUE,
                   s4 = TRUE,
                   s5 = TRUE,
                   s6 = TRUE,
                   s7 = FALSE,
                   s8 = TRUE,
                   s8a = TRUE,
                   s9 = TRUE)
  
  sw <- list(s1 = "on",
             s2 = "on",
             s3 = "on",
             s4 = "on",
             s5 = "on",
             s6 = "on",
             s7 = "off",
             s8 = "off",
             s8a = "on",
             s9 = "on")
  
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
  
  filter.thresholds <- list(nFeatureRNAfloor = NULL,
                            nFeatureRNAceiling = NULL,
                            nCountRNAfloor = NULL, 
                            nCountRNAceiling = NULL,
                            pct_mitofloor = NULL, 
                            pct_mitoceiling = 10,
                            pct_ribofloor = NULL, 
                            pct_riboceiling = NULL,
                            ambientRNA_thres = 0.5)
  
  remove_doublet <- TRUE
  path.to.10X.doublet.estimation <- file.path(path.to.storage, "DoubletEstimation10X.csv")
  
  num.PCA <- 25
  num.PC.used.in.UMAP <- 25
  num.PC.used.in.Clustering <- 25
  
  for (sample.id in names(stage_lst)){
    path.to.output <- file.path(path.to.main.output, sprintf("%s_round", analysis.round))
    dir.create(path.to.output, showWarnings = FALSE)
    
    path.to.anno.contigs <- file.path(path.to.VDJ.output, sprintf("annotated_contigs_clonaltype_%s.csv", sample.id))
    path.to.count.clonaltype <- file.path(path.to.VDJ.output, sprintf("count_clonaltype_%s.csv", sample.id))
    
    filtered.barcodes <- NULL
    
    s.obj <- run_pipeline_GEX(path2src=path2src,
                              path2input=file.path(path2input, sample.id),
                              path.to.logfile.dir=file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round), "logs"),
                              stage_lst=stage_lst[[sample.id]],
                              path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                              MINCELLS=MINCELLS,
                              MINGENES=MINGENES,
                              PROJECT=PROJECT,
                              remove_doublet=remove_doublet,
                              save.RDS=save.RDS,
                              path.to.output=file.path(path.to.output, sprintf("%s_%s_round", sample.id, analysis.round)),
                              rerun=rerun,
                              DE.test="wilcox",
                              num.PCA=num.PCA,
                              num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                              num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                              use.sctransform=FALSE,
                              filtered.barcodes=filtered.barcodes,
                              filter.thresholds=filter.thresholds,
                              path.to.anno.contigs=path.to.anno.contigs,
                              path.to.count.clonaltype=path.to.count.clonaltype,
                              input.method = "filterIg",
                              sw = sw)
  }
}



