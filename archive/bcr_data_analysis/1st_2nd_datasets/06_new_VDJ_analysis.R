gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/00_helper_functions.R")
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
} 
library("Biostrings")
new.pkgs <- c("msa")
for (pkg in new.pkgs){
  if (pkg %in% installed.packages() ==  FALSE){
    BiocManager::install(pkg, update = FALSE)
  }
}
library("msa")

#####----------------------------------------------------------------------#####
##### GENERATE VDJ DATA WITH THRESHOLD 0.85
#####----------------------------------------------------------------------#####
outdir <- "/media/hieunguyen/HNSD_mini/outdir"
PROJECT <- "1st_2nd_BSimons_Datasets"

path.to.storage <- "/media/hieunguyen/GSHD_HN01/storage"
thres <- 0.85

path.to.main.input <- file.path(path.to.storage, PROJECT)
path.to.main.output <- file.path(outdir, PROJECT, "data_analysis")
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
path.to.VDJ.output <- file.path(path.to.main.output, sprintf("VDJ_output_%s", thres))
dir.create(path.to.VDJ.output, showWarnings = FALSE, recursive = TRUE)

if (file.exists(file.path(path.to.VDJ.output, sprintf("finished_saving_VDJ_thres_%s", thres))) == FALSE){
  print(sprintf("Generating VDJ data analysis with threshold %s", thres))
  source("/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_VDJ_pipeline/main_VDJ_pipeline.R")
  path.to.VDJ.input <- file.path(path.to.main.input, "VDJ")
  
  summarize_vdj_data(path.to.VDJ.input, 
                     path.to.VDJ.output, 
                     PROJECT, 
                     removeNA=FALSE, 
                     removeMulti=FALSE, 
                     T_or_B = "B",
                     threshold = thres)
  write.csv(data.frame(status = c("finished")), 
            file.path(path.to.VDJ.output, sprintf("finished_saving_VDJ_thres_%s", thres)))
} 

source("/media/hieunguyen/HNSD01/src/bcr_data_analysis/1st_2nd_datasets/get_germline_sequences.R")
names(s.V.genes) <- to_vec(
  for (item in names(s.V.genes)) str_split(item, "[*]")[[1]][[1]]
)
names(s.J.genes) <- to_vec(
  for (item in names(s.J.genes)) str_split(item, "[*]")[[1]][[1]]
)

names(cellranger.ref) <- to_vec(
  for (item in names(cellranger.ref)) str_split(str_split(item, " ")[[1]][[1]], "[|]")[[1]][[2]]
)
#####----------------------------------------------------------------------#####
##### read file from VDJ output, after pre-processing. 
#####----------------------------------------------------------------------#####
all.VDJ.files <- Sys.glob(file.path(path.to.VDJ.output, "annotated_contigs*.csv"))
names(all.VDJ.files) <- to_vec(
  for (item in all.VDJ.files) str_replace(str_replace(basename(item), "annotated_contigs_clonaltype_", ""), ".csv", "")
)
path.to.output <- file.path(outdir, PROJECT, "data_analysis")
path.to.save.output <- file.path(path.to.output, "VDJ_data_output")
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
##### read file "filtered_contig_annotations.csv" from the raw CellRanger pipeline
#####----------------------------------------------------------------------#####
path.to.main.input <- "/media/hieunguyen/GSHD_HN01/storage/1st_2nd_BSimons_Datasets/VDJ"
all.raw.VDJ.files <- Sys.glob(file.path(path.to.main.input, "*", "filtered_contig_annotations.csv"))
names(all.raw.VDJ.files) <- to_vec(
  for (item in all.raw.VDJ.files) basename(dirname(item))
)
all.raw.VDJ.files <- all.raw.VDJ.files[names(all.VDJ.files)]
all.raw.VDJ.fasta <- Sys.glob(file.path(path.to.main.input, "*", "filtered_contig.fasta"))
names(all.raw.VDJ.fasta) <- to_vec(
  for (item in all.raw.VDJ.fasta) basename(dirname(item))
)
all.raw.VDJ.fasta <- all.raw.VDJ.fasta[names(all.VDJ.files)]

#####----------------------------------------------------------------------#####
##### get metadata sheets and single cell data from the 1st and 2nd datasets
#####----------------------------------------------------------------------#####
path.to.06.output <- file.path(path.to.main.output, "data_analysis", "06_output", sprintf("VDJ_thres_%s", thres))
dir.create(path.to.06.output, showWarnings = FALSE, recursive = TRUE)
s.obj.1st <- readRDS(file.path(path.to.main.output, 
                               "01_output",
                               "1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1", 
                               "1st_dataset_removed_7_9_removed_16_removed_9.with_reInt.res1.rds"))
metadata.1st <- s.obj.1st@meta.data %>% rownames_to_column("barcode")

s.obj.2nd <- readRDS(file.path(path.to.main.output, 
                               "01_output", 
                               "2nd_dataset_removed_5_6.without_reInt.res1", 
                               "2nd_dataset_removed_5_6.without_reInt.res1.rds"))
metadata.2nd <- s.obj.2nd@meta.data %>% rownames_to_column("barcode")

all.samples <- c(unique(s.obj.1st$name), unique(s.obj.2nd$name))

#####----------------------------------------------------------------------#####
##### Generate clone dataframe and calculate mutation rate for each clone/cell
#####----------------------------------------------------------------------#####
path.to.clone.dfs <- file.path(path.to.main.output, "VDJ_data_output")

if (file.exists(file.path(path.to.06.output, "clonedf_with_mutation_rate.xlsx")) == FALSE){
  clonedf <- data.frame()
  for (sample.id in all.samples){
    tmpdf <- readxl::read_excel(file.path(path.to.clone.dfs, sprintf("%s.xlsx", sample.id)))
    tmpdf$SampleID <- sample.id
    clonedf <- rbind(clonedf, tmpdf)
  }
  
  nt_cols <- to_vec(
    for (item in colnames(clonedf)) if (grepl("_nt", item) == TRUE) item
  )
  clonedf.raw <- clonedf 
  
  ##### keep IGH clonedf only
  clonedf <- subset(clonedf, clonedf$chain == "IGH")
  
  clonedf$num.mutation <- unlist(lapply(
    seq(1, nrow(clonedf)), function(i){
      V.gene <- clonedf[i, ]$v_gene
      J.gene <- clonedf[i, ]$j_gene
      CDR3.seq <- clonedf[i, ]$cdr3_nt
      CDR3.length <- nchar(CDR3.seq)
      clone.seq <- clonedf[i, ]$prep.full.seq
      
      GL.V.gene.IMGT <- s.V.genes[[V.gene]] %>% as.character()
      GL.J.gene.IMGT <- s.J.genes[[J.gene]] %>% as.character()
      
      GL.V.gene.10x <- cellranger.ref[[V.gene]] %>% as.character()
      GL.J.gene.10x <- cellranger.ref[[J.gene]] %>% as.character()
      
      repN.seq <- paste(replicate(n = CDR3.length, expr = "N"), collapse = "")
      
      GL.seq.IMGT <- sprintf("%s%s%s", GL.V.gene.IMGT, repN.seq, GL.J.gene.IMGT)
      input.seqs.IMGT <- c(clone.seq, GL.seq.IMGT) %>% DNAStringSet()
      msa.IMGT <- msa(input.seqs.IMGT, method = "Muscle")
      aligned.clone.seq.IMGT <- toString(unmasked(msa.IMGT)[[1]])
      aligned.GL.seq.IMGT <- toString(unmasked(msa.IMGT)[[2]])
      
      GL.seq.10x <- sprintf("%s%s%s", GL.V.gene.10x, repN.seq, GL.J.gene.10x)
      input.seqs.10x <- c(clone.seq, GL.seq.10x) %>% DNAStringSet()
      msa.10x <- msa(input.seqs.10x, method = "Muscle")
      aligned.clone.seq.10x <- toString(unmasked(msa.10x)[[1]])
      aligned.GL.seq.10x <- toString(unmasked(msa.10x)[[2]])
      
      s1 <- aligned.clone.seq.10x
      s2 <- aligned.GL.seq.10x
      sdf <- data.frame(s1 = str_split(s1, "")[[1]],
                        s2 = str_split(s2, "")[[1]])
      
      check_mutation <- function(x1, x2){
        if (x1 == "-" | x2 == "-"){
          output <- "undefined"
        } else if (x1 == "N" | x2 == "N"){
          output <- "undefined"
        } else {
          if (x1 == x2){
            output <- "matched"
          } else {
            output <- "mutation"
          }
        }
        return(output)
      }
      
      sdf <- sdf %>% rowwise() %>% mutate(status = check_mutation(s1, s2))
      count.mutation <- table(sdf$status)
      if ("mutation" %in% names(count.mutation) == FALSE){
        return(0)
      } else {
        return(count.mutation[["mutation"]])      
      }
    }
  ))
  
  writexl::write_xlsx(clonedf, file.path(path.to.06.output, "clonedf_with_mutation_rate.xlsx"))
} else {
  clonedf <- readxl::read_excel(file.path(path.to.06.output, "clonedf_with_mutation_rate.xlsx"))
}

#####----------------------------------------------------------------------#####
##### 06.1: SCATTER PLOT: mutation rate and clone size
#####----------------------------------------------------------------------#####
for (sample.id in unique(clonedf$SampleID)){
  tmpdf <- clonedf[, c("barcode", "CTstrict", "num.mutation", "SampleID")] %>%
    subset(SampleID == sample.id)
  
  mean.mutationdf <- tmpdf %>% group_by(CTstrict) %>% 
    summarise(mean.num.mutation = mean(num.mutation, na.rm = TRUE)) 
  colnames(mean.mutationdf) <- c("clone", "avg_num_mutation")
  mean.mutationdf <- mean.mutationdf %>% rowwise() %>% 
    mutate(clone.size = nrow(subset(tmpdf, tmpdf$CTstrict == clone)))
  
  p <- mean.mutationdf %>% ggplot(aes(x = clone.size, y = avg_num_mutation)) + geom_point()
  dir.create(file.path(path.to.06.output, "06.1_scatter_plot"), showWarnings = FALSE, recursive = TRUE)
  ggsave(plot = p, 
         filename = sprintf("scatterplot_mutatino_rate_vs_clone_size.%s.svg", sample.id), 
         path = file.path(path.to.06.output, "06.1_scatter_plot"),
         device = "svg", 
         dpi = 300, 
         width = 14, 
         height = 10)
}

#####----------------------------------------------------------------------#####
##### 06.1: SCATTER PLOT: mutation rate and clone size
#####----------------------------------------------------------------------#####