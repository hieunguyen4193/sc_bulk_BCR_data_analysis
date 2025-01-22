gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### packages
#####----------------------------------------------------------------------#####
path.to.main.src <- "/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis"
source(file.path(path.to.main.src, "VDJ_helper_functions.R"))
scrna_pipeline_src <- "/media/hieunguyen/HNSD01/src/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(scrna_pipeline_src, "import_libraries.R"))
source(file.path(scrna_pipeline_src, "helper_functions.R"))
source(file.path(scrna_pipeline_src, "s8_integration_and_clustering.R"))
library(circlize)

#####----------------------------------------------------------------------#####
##### INPUT ARGS
#####----------------------------------------------------------------------#####
path.to.storage <- "/media/hieunguyen/HNHD01/storage/all_BSimons_datasets"
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.addedClone.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE

# circos.group.type <- "VJnt"
circos.group.type <- "VJaa"

#####----------------------------------------------------------------------#####
##### READ METADATA
#####----------------------------------------------------------------------#####
bulk.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241031_BSimons/241031_sample_sheet.xlsx") %>%
  rowwise() %>% 
  mutate(MID = sprintf("MID%s", MID))
sc.projects <- c("241002_241104_BSimons")
bulk.projects <- c("241031_BSimons")
list.of.PROJECT <- c(sc.projects, bulk.projects)

path.to.05.output <- file.path(outdir, "VDJ_output", "05_output", paste(list.of.PROJECT, collapse = "_"))
dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(path.to.05.output, circos.group.type), showWarnings = FALSE, recursive = TRUE)

#####----------------------------------------------------------------------#####
#### read the Seurat object of single cell dataset
######----------------------------------------------------------------------#####
s.obj <- list()
sc.meta.data <- list()
list.of.ht <- list()
for (PROJECT in sc.projects){
  s.obj[[PROJECT]] <- readRDS(path.to.all.s.obj[[PROJECT]])
  sc.meta.data[[PROJECT]] <- s.obj[[PROJECT]]@meta.data %>% rownames_to_column("barcode")
  list.of.ht[[PROJECT]] <- list()
  for (sample.id in unique(sc.meta.data[[PROJECT]]$name)){
    list.of.ht[[PROJECT]][[sample.id]] <- unique(subset(sc.meta.data[[PROJECT]], sc.meta.data[[PROJECT]]$name == sample.id)$HTO_classification)
  }
}

#####----------------------------------------------------------------------#####
##### READ CLONE DATA -----> NEW DATA
#####----------------------------------------------------------------------#####
if (file.exists(file.path(path.to.05.output, circos.group.type, "all_data.rds")) == FALSE){
  print("GENERATING NEW DATA!!!!!")
  all.data <- list()
  for (PROJECT in list.of.PROJECT){
    path.to.VDJ.output <- file.path( outdir, "VDJ_output", PROJECT, sprintf("VDJ_output_%s", thres), "preprocessed_files")
    
    path.to.mid.output <- file.path(path.to.storage, PROJECT, "mixcr_pipeline_output/v0.2/mid_based_output")
    
    if (PROJECT %in% bulk.projects){
      clone.obj <- run_preprocessing_all_bulk_VDJ_data(path.to.mid.output = path.to.mid.output,
                                                       path.to.save.output = path.to.05.output,
                                                       PROJECT = PROJECT,
                                                       thres = thres, 
                                                       thres.dis = thres.dis,
                                                       savefile = savefile,
                                                       verbose = verbose,
                                                       rerun = rerun,
                                                       define.clone.clusters = define.clone.clusters) 
    } else if (PROJECT %in% sc.projects){
      clone.obj <- run_preprocessing_all_sc_data( path.to.VDJ.output = path.to.VDJ.output, 
                                                  path.to.save.output = path.to.05.output, 
                                                  PROJECT = PROJECT,
                                                  thres = thres, 
                                                  thres.dis = thres.dis,
                                                  savefile = savefile,
                                                  rerun = rerun,
                                                  define.clone.clusters =  define.clone.clusters)
    }
    
    full.clonedf <- clone.obj$clonesets
    full.clonedf <- full.clonedf %>%
      rowwise() %>%
      mutate(V.gene = ifelse(grepl("[*]", V.gene) == TRUE, str_split(V.gene, "[*]")[[1]][[1]], V.gene )) %>%
      mutate(J.gene = ifelse(grepl("[*]", J.gene) == TRUE, str_split(J.gene, "[*]")[[1]][[1]], J.gene ))
    full.clonedf <- subset(full.clonedf, full.clonedf$aaSeqCDR3 != "region_not_covered")
    if (PROJECT %in% sc.projects){
      full.clonedf <- full.clonedf %>% rowwise() %>%
        mutate(barcode.full = sprintf("%s_%s", id, barcode)) %>%
        # keep only cells that have hashtag information and remain in the 
        # filtered seurat object data.
        mutate(HT = ifelse(barcode.full %in% sc.meta.data[[PROJECT]]$barcode, 
                           subset(sc.meta.data[[PROJECT]], sc.meta.data[[PROJECT]]$barcode == barcode.full)$HTO_classification,
                           NA)) %>%
        subset(is.na(HT) == FALSE) %>%
        mutate(id.origin = id) %>%
        mutate(id = sprintf("%s_%s", id, HT))
    }
    
    dir.create(file.path(path.to.05.output, circos.group.type), showWarnings = FALSE, recursive = TRUE)
    
    ##### split the full clone dataframe to smaller dataframe for each MID/each single cell sample hashtag
    for (mid in unique(full.clonedf$id)){
      if (file.exists(file.path(path.to.05.output, 
                                circos.group.type, 
                                sprintf("%s.simplified.csv", mid))) == FALSE){
        
        print(sprintf("Working on sample MID %s", mid))
        clonedf <- subset(full.clonedf, full.clonedf$id == mid)
        input.circos <- subset(clonedf, select = c(V.gene, 
                                                   J.gene, 
                                                   aaSeqCDR3, 
                                                   nSeqCDR3, 
                                                   uniqueMoleculeCount)) %>%
          rowwise() %>%
          mutate(VJnt = sprintf("%s_%s_%s", V.gene, J.gene, nSeqCDR3)) %>%
          mutate(VJaa = sprintf("%s_%s_%s", V.gene, J.gene, aaSeqCDR3))
        if (PROJECT %in% sc.projects){
          input.circos <- data.frame(table(input.circos[[circos.group.type]]))
          colnames(input.circos) <- c("id", "cloneCount")          
        } else if (PROJECT %in% bulk.projects){
          input.circos <- input.circos[, c(circos.group.type, "uniqueMoleculeCount")] 
          colnames(input.circos) <- c("id", "uniqueMoleculeCount")      
          input.circos <- input.circos %>% 
            group_by(id) %>%
            summarise(cloneCount = sum(uniqueMoleculeCount))
          input.circos <- subset(input.circos, select = c(id, cloneCount)) 
        }
        
        input.circos <- input.circos %>% rowwise() %>%
          mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
          mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
          mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]]) %>%
          arrange(desc(cloneCount))
        write.table(input.circos, 
                    file.path(path.to.05.output, 
                              circos.group.type, 
                              sprintf("%s.simplified.csv", mid)), 
                    quote = FALSE, 
                    sep = "\t", 
                    row.names = FALSE) 
        all.data[[sprintf("%s_%s", PROJECT, mid)]] <- input.circos
      } else {
        print(sprintf("File exists at %s, reading in ...", 
                      file.path(path.to.05.output, circos.group.type, sprintf("%s.simplified.csv", mid))))
        all.data[[sprintf("%s_%s", PROJECT, mid)]] <- read.csv(file.path(path.to.05.output, circos.group.type, sprintf("%s.simplified.csv", mid)), sep = "\t")
      }
    }
    
    ##### if the dataset is a single cell dataset, get clonedf for all hashtag in a mouse sample. 
    if (PROJECT %in% sc.projects){
      for (mid in unique(full.clonedf$id.origin)){
        if (file.exists(file.path(path.to.05.output, 
                                  circos.group.type, 
                                  sprintf("%s.simplified.csv", mid))) == FALSE){
          
          print(sprintf("Working on sample MID %s", mid))
          clonedf <- subset(full.clonedf, full.clonedf$id.origin == mid)
          input.circos <- subset(clonedf, select = c(V.gene, 
                                                     J.gene, 
                                                     aaSeqCDR3, 
                                                     nSeqCDR3, 
                                                     uniqueMoleculeCount)) %>%
            rowwise() %>%
            mutate(VJnt = sprintf("%s_%s_%s", V.gene, J.gene, nSeqCDR3)) %>%
            mutate(VJaa = sprintf("%s_%s_%s", V.gene, J.gene, aaSeqCDR3)) 
          
          input.circos <- data.frame(table(input.circos[[circos.group.type]]))
          colnames(input.circos) <- c("id", "cloneCount")
          input.circos <- input.circos %>% rowwise() %>%
            mutate(bestVHit = str_split(id, "_")[[1]][[1]]) %>%
            mutate(bestJHit = str_split(id, "_")[[1]][[2]]) %>%
            mutate(nSeqCDR3 = str_split(id, "_")[[1]][[3]]) %>%
            arrange(desc(cloneCount))
          write.table(input.circos, 
                      file.path(path.to.05.output, 
                                circos.group.type, 
                                sprintf("%s.simplified.csv", mid)), 
                      quote = FALSE, 
                      sep = "\t", 
                      row.names = FALSE) 
          all.data[[sprintf("%s_%s", PROJECT, mid)]] <- input.circos
        } else {
          print(sprintf("File exists at %s, reading in ...", 
                        file.path(path.to.05.output, circos.group.type, sprintf("%s.simplified.csv", mid))))
          all.data[[sprintf("%s_%s", PROJECT, mid)]] <- read.csv(file.path(path.to.05.output, circos.group.type, sprintf("%s.simplified.csv", mid)), sep = "\t")
        }
      }
    }
  }
  saveRDS(all.data, file.path(path.to.05.output, circos.group.type, "all_data.rds"))
  meta.data <- data.frame(MID = names(all.data))
  meta.data <- meta.data  %>% rowwise() %>% 
    mutate(PROJECT = ifelse(grepl("_HT", MID), 
                            paste(str_split(MID, "_")[[1]][1: (length(str_split(MID, "_")[[1]]) - 2)], collapse = "_"),
                            paste(str_split(MID, "_")[[1]][1: (length(str_split(MID, "_")[[1]]) - 1)], collapse = "_")) ) %>%
    mutate(SampleID = str_replace(MID, sprintf("%s_", PROJECT), "")) %>%
    mutate(mouse = ifelse(
      PROJECT %in% sc.projects,
      sprintf("m%s", str_split(SampleID, "")[[1]][[3]]),
      subset(bulk.metadata, bulk.metadata$MID == SampleID)$mouse
    )) %>%
    mutate(organ = ifelse(
      PROJECT %in% sc.projects,
      paste0(str_split(SampleID, "")[[1]][1:2], collapse = ""),
      subset(bulk.metadata, bulk.metadata$MID == SampleID)$organ
    ))
  writexl::write_xlsx(meta.data, file.path(path.to.05.output, "all_data_metadata.xlsx"))
} else {
  print(sprintf("All samples data exists, reading in ..."))
  all.data <- readRDS(file.path(path.to.05.output, circos.group.type, "all_data.rds"))
  meta.data <- readxl::read_excel(file.path(path.to.05.output, "all_data_metadata.xlsx"))
}

#####----------------------------------------------------------------------#####
##### PREPROCESS FILE: ASSIGN CLONES TO CLUSTERS 0.85 SIMILARITY
#####----------------------------------------------------------------------#####
all.input.files <- Sys.glob(file.path(path.to.05.output,
                                      circos.group.type,
                                      "*"))

input.metadata <- data.frame(
  path = all.input.files,
  SampleID = to_vec(for (item in all.input.files){
    str_replace(basename(item), ".simplified.csv", "") 
  }))

all.input.files <- input.metadata$path
names(all.input.files) <- input.metadata$SampleID

##### generate circos plot for all hashtags
exclude.samples <- c("PP3", "PP7")

meta.data.splitted.or.not <- list(
  with_hashtags = subset(meta.data, meta.data$SampleID %in% exclude.samples == FALSE),
  without_hashtags = subset(meta.data, grepl("_", meta.data$SampleID) == FALSE)
)


for (meta.data.name in names(meta.data.splitted.or.not)){
  tmp.metadata <- meta.data.splitted.or.not[[meta.data.name]]
  dir.create(file.path(path.to.05.output, sprintf("VJcombi_CDR3_%s", thres), meta.data.name), showWarnings = FALSE, recursive = TRUE)
  
  if (file.exists(file.path(path.to.05.output, sprintf("VJcombi_CDR3_%s", thres), meta.data.name, "all.data.VJ.combi.rds")) == FALSE){
    all.data.VJ.combi <- list()
    for (mouse.id in c("m3", "m7")){
      selected.mids <- subset(tmp.metadata, tmp.metadata$mouse == mouse.id)$SampleID
      input.files <- all.input.files[selected.mids]
      
      fileAliases <- to_vec(
        for (item in names(input.files)){
          sprintf("%s (%s)", item, subset(tmp.metadata, tmp.metadata$SampleID == item)$organ)
        }
      )
      saveFileName <- sprintf("%s_hashtags_circos.svg", mouse.id)
      outputdir <- file.path(path.to.05.output,
                             circos.group.type,
                             "circos_plot")
      filter.clone <- FALSE
      filter.clone.cutoff <- NA
      source(file.path(path.to.main.src, "circos_helper.R"))
      thres.dis <- 0.15
      thres <- 0.85
      
      clonesets <- read_tsv(input.files, id = "fileName") %>%
        rowwise() %>%
        mutate(fileName = basename(fileName) %>% str_replace(".simplified.csv", "")) %>%
        # mutate(bestVHit = str_replace_all(bestVHit, "[*]", "-")) %>%
        # mutate(bestJHit = str_replace_all(bestJHit, "[*]", "-")) %>%
        mutate(bestVHit = str_split(bestVHit, "[*]")[[1]][[1]]) %>%
        mutate(bestJHit = str_split(bestJHit, "[*]")[[1]][[1]]) %>%
        mutate(len = nchar(nSeqCDR3)) %>%
        mutate(VJ.len.combi = sprintf("%s_%s_%s", bestVHit, bestJHit, len ))
      
      colnames(clonesets) <- c("fileName", "id", "cloneCount", "bestVHit", "bestJHit", "seq", "len", "VJ.len.combi")
      ##### Group sequences + Gene usages to clones
      new.clonesets <- data.frame()
      for (input.VJ.combi in unique(clonesets$VJ.len.combi)){
        tmpdf <- subset(clonesets, clonesets$VJ.len.combi == input.VJ.combi)
        seqs <- unique(tmpdf$seq)
        print(sprintf("VJ.len.combi: %s, num seqs: %s", input.VJ.combi, length(seqs)))
        if (length(seqs) >= 2){
          cluster.output <- assign_clusters_to_sequences(seqs = seqs, threshold = thres.dis)$res
          tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- unlist(lapply(
            tmpdf$seq, function(x){
              return(sprintf("%s_%s", input.VJ.combi, subset(cluster.output, cluster.output$seq == x)$cluster))
            }
          ))    
        } else {
          tmpdf[[sprintf("VJcombi_CDR3_%s", thres)]] <- sprintf("%s_1", input.VJ.combi)
        }
        new.clonesets <- rbind(new.clonesets, tmpdf)
      }
      
      for (input.mid in unique(new.clonesets$fileName)){
        tmpdf <- subset(new.clonesets, new.clonesets$fileName == input.mid)[, c(sprintf("VJcombi_CDR3_%s", thres),
                                                                                "cloneCount", 
                                                                                "bestVHit",
                                                                                "bestJHit",
                                                                                "seq")]
        colnames(tmpdf) <- c("id", "cloneCount", "bestVHit", "bestJHit", "nSeqCDR3")
        write.table(tmpdf, 
                    file.path(path.to.05.output, 
                              sprintf("VJcombi_CDR3_%s", thres), 
                              meta.data.name,
                              sprintf("%s.simplified.csv", input.mid)), 
                    quote = FALSE, 
                    sep = "\t", 
                    row.names = FALSE) 
        all.data.VJ.combi[[input.mid]] <- tmpdf
      }
    }
    saveRDS(all.data.VJ.combi, file.path(path.to.05.output, sprintf("VJcombi_CDR3_%s", thres), meta.data.name, "all.data.VJ.combi.rds"))
  }
}

for (meta.data.name in names(meta.data.splitted.or.not)){
  tmp.metadata <- meta.data.splitted.or.not[[meta.data.name]]
  all.input.files <- Sys.glob(file.path(path.to.05.output, 
                                        sprintf("VJcombi_CDR3_%s", thres), 
                                        meta.data.name,
                                        "*.simplified.csv"))
  
  input.metadata <- data.frame(
    path = all.input.files,
    SampleID = to_vec(for (item in all.input.files){
      str_replace(basename(item), ".simplified.csv", "") 
    }),
    PROJECT = to_vec(for (item in all.input.files){
      str_split(item, "/")[[1]][[8]]
    })
  ) 
  
  all.input.files <- input.metadata$path
  names(all.input.files) <- input.metadata$SampleID
  
  for (mouse.id in c("m3", "m7")){
    selected.mids <- subset(tmp.metadata, tmp.metadata$mouse == mouse.id)$SampleID
    
    if (meta.data.name == "with_hashtags"){
      if (mouse.id == "m3"){
        p.sample <- c("PP3_HT1", "PP3_HT2", "PP3_HT3")
      } else if (mouse.id == "m7"){
        p.sample <- c("PP7_HT3", "PP7_HT1", "PP7_HT2")
      }
    } else if (meta.data.name == "without_hashtags"){
      p.sample <- selected.mids[grepl("P", selected.mids)]
      p.sample <- sort(p.sample, decreasing = TRUE)
    }
    
    MID.samples <- selected.mids[grepl("MID", selected.mids) == TRUE]
    
    ordered.selected.mids <- c(p.sample, MID.samples)
    input.files <- all.input.files[ordered.selected.mids]
    
    group.to.highlight1 <- subset(tmp.metadata, tmp.metadata$mouse == mouse.id & organ %in% c("PP"))$SampleID
    group.to.highlight2 <- subset(tmp.metadata, tmp.metadata$mouse == mouse.id & organ %in% c("PP") == FALSE)$SampleID
    
    fileAliases <- to_vec(
      for (item in names(input.files)){
        sprintf("%s (%s)", item, subset(tmp.metadata, tmp.metadata$SampleID == item)$organ)
        # convert.sampleID[[item]]
      }
    )
    
    names(fileAliases) <- names(input.files)
    if (meta.data.name == "with_hashtags"){
      saveFileName <- sprintf("%s_hashtags_circos.svg", mouse.id)
    } else {
      saveFileName <- sprintf("%s_circos.svg", mouse.id)
    }
    outputdir <- file.path(path.to.05.output,
                           sprintf("VJcombi_CDR3_%s", thres),
                           "circos_plot")
    filter.clone <- FALSE
    filter.clone.cutoff <- NA
    source(file.path(path.to.main.src, "circos_helper.R"))
    
    if (file.exists(file.path(outputdir, saveFileName)) == FALSE){
      generate_circos(
        input = input.files,
        fileAliases = fileAliases,
        saveFileName = saveFileName,
        outputdir = outputdir,
        filter.clone = filter.clone,
        filter.clone.cutoff = filter.clone.cutoff,
        group.to.highlight1 = group.to.highlight1,
        group.to.highlight2 = group.to.highlight2,
        linkColor1 = "#FF000080",
        linkColor2 = "lightgray",
        ordered.samples = ordered.selected.mids
      )
    }
  }
}
