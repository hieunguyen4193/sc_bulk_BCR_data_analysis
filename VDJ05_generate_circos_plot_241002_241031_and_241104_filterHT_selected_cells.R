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
source(file.path(path.to.main.src, "GEX_path_to_seurat_obj.R"))
outdir <- "/media/hieunguyen/GSHD_HN01/outdir/sc_bulk_BCR_data_analysis_v0.1"

thres <- 0.85
thres.dis <- 0.15
savefile <- TRUE
verbose <- TRUE
rerun <- FALSE
define.clone.clusters <- FALSE
cut.off.ht <- 40
for (cell.group in c("group1", "group2", "group3")){
  for (circos.group.type in c("VJnt", "VJaa")){
    for (filter.ht.type in c("remove_last_ht", sprintf("cutoff_%s", cut.off.ht))){
      #####----------------------------------------------------------------------#####
      ##### READ METADATA
      #####----------------------------------------------------------------------#####
      bulk.metadata <- readxl::read_excel("/media/hieunguyen/HNSD01/src/sc_bulk_BCR_data_analysis/preprocessing/241031_BSimons/241031_sample_sheet.xlsx") %>%
        rowwise() %>%
        mutate(MID = sprintf("MID%s", MID))
      sc.projects <- c("241002_BSimons", "241104_BSimons")
      bulk.projects <- c("241031_BSimons")
      list.of.PROJECT <- c(sc.projects, bulk.projects)
      
      path.to.05.output <- file.path(outdir, "VDJ_output", 
                                     sprintf("05_output_filterHT_%s_selected_cells_%s", filter.ht.type, cell.group), 
                                     paste(list.of.PROJECT, collapse = "_"))
      dir.create(path.to.05.output, showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(path.to.05.output, circos.group.type), showWarnings = FALSE, recursive = TRUE)
      
      #####----------------------------------------------------------------------#####
      #### read the Seurat object of single cell dataset
      ######----------------------------------------------------------------------#####
      s.obj <- list()
      sc.meta.data <- list()
      keep.ht <- list()
      keep.barcodes <- list()
      for (PROJECT in sc.projects){
        keep.ht[[PROJECT]] <- list()
        keep.ht[[PROJECT]][["remove_last_ht"]] <- c()
        s.obj[[PROJECT]] <- readRDS(path.to.all.s.obj[[PROJECT]])
        tmp.metadata <- s.obj[[PROJECT]]@meta.data %>% rownames_to_column("barcode") %>%
          rowwise() %>%
          mutate(sample_HT = sprintf("%s_%s", name, HTO_classification))
        count.ht <- table(tmp.metadata$sample_HT)
        list.of.ht <- list()
        for (sampleid in unique(tmp.metadata$name)){
          list.of.ht[[sampleid]] <- subset(tmp.metadata, tmp.metadata$name == sampleid)$HTO_classification %>% unique() %>% sort()
          keep.ht[[PROJECT]][["remove_last_ht"]] <- c(keep.ht[[PROJECT]][["remove_last_ht"]],
                                                      to_vec( 
                                                        for (item in list.of.ht[[sampleid]][1:length(list.of.ht[[sampleid]]) - 1]){
                                                          sprintf("%s_%s", sampleid, item)
                                                        }
                                                      ))
          keep.barcodes[[PROJECT]][[cell.group]] <- read.csv(file.path(outdir, 
                                                                       "GEX_output", 
                                                                       "06_output", 
                                                                       PROJECT, 
                                                                       sprintf("cell_group_%s.csv", cell.group)))$barcode
        }
        sc.meta.data[[PROJECT]] <- tmp.metadata
        keep.ht[[PROJECT]][[sprintf("cutoff_%s", cut.off.ht)]] <- count.ht[count.ht >= cut.off.ht] %>% names()
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
                                                             path.to.save.output = path.to.VDJ.output,
                                                             PROJECT = PROJECT,
                                                             thres = thres, 
                                                             thres.dis = thres.dis,
                                                             savefile = savefile,
                                                             verbose = verbose,
                                                             rerun = rerun,
                                                             define.clone.clusters = define.clone.clusters) 
          } else if (PROJECT %in% sc.projects){
            clone.obj <- run_preprocessing_all_sc_data( path.to.VDJ.output = path.to.VDJ.output, 
                                                        path.to.save.output = path.to.VDJ.output, 
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
              mutate(id = sprintf("%s_%s", id, HT)) %>%
              subset(id %in% keep.ht[[PROJECT]][[filter.ht.type]]) %>%
              subset(barcode.full %in% keep.barcodes[[PROJECT]][[cell.group]])
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
                                                         uniqueMoleculeCount,
                                                         VJ.len.combi)) %>%
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
                                                           uniqueMoleculeCount,
                                                           VJ.len.combi)) %>%
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
          mutate(PROJECT = paste(str_split(MID, "_")[[1]][1:2], collapse = "_") ) %>%
          mutate(SampleID = str_replace(MID, sprintf("%s_", PROJECT), "")) %>%
          mutate(mouse = ifelse(
            PROJECT %in% sc.projects,
            sprintf("m%s", str_split(SampleID, "")[[1]][[3]]),
            subset(bulk.metadata, bulk.metadata$MID == SampleID)$mouse
          )) %>%
          mutate(organ = ifelse(
            PROJECT %in% sc.projects,
            str_split(SampleID, "")[[1]][[1]],
            subset(bulk.metadata, bulk.metadata$MID == SampleID)$organ
          ))
        writexl::write_xlsx(meta.data, file.path(path.to.05.output, "all_data_metadata.xlsx"))
      } else {
        print(sprintf("All samples data exists, reading in ..."))
        all.data <- readRDS(file.path(path.to.05.output, circos.group.type, "all_data.rds"))
        meta.data <- readxl::read_excel(file.path(path.to.05.output, "all_data_metadata.xlsx"))
      }
      
      #####----------------------------------------------------------------------#####
      ##### MAIN FUNCTIONS GENERATE CIRCOS PLOT
      #####----------------------------------------------------------------------#####
      all.input.files <- Sys.glob(file.path(path.to.05.output,
                                            circos.group.type,
                                            "*"))
      
      input.metadata <- data.frame(
        path = all.input.files,
        SampleID = to_vec(for (item in all.input.files){
          str_replace(basename(item), ".simplified.csv", "") 
        }),
        PROJECT = to_vec(for (item in all.input.files){
          str_split(item, "/")[[1]][[8]]
        }))
      
      all.input.files <- input.metadata$path
      names(all.input.files) <- input.metadata$SampleID
      
      ##### generate circos plot for all hashtags
      exclude.samples <- c("PP3", "PP7")
      meta.data.splitted <- subset(meta.data, meta.data$SampleID %in% exclude.samples == FALSE)
      meta.data.non.splitted <- subset(meta.data, grepl("_", meta.data$SampleID) == FALSE)
      
      for (mouse.id in c("m3", "m7")){
        selected.mids <- subset(meta.data.splitted, meta.data.splitted$mouse == mouse.id)$SampleID
        input.files <- all.input.files[selected.mids]
        
        fileAliases <- to_vec(
          for (item in names(input.files)){
            sprintf("%s (%s)", item, subset(meta.data.splitted, meta.data.splitted$SampleID == item)$organ)
          }
        )
        names(fileAliases) <- names(input.files)
        saveFileName <- sprintf("%s_hashtags_circos.svg", mouse.id)
        outputdir <- file.path(path.to.05.output,
                               circos.group.type,
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
            filter.clone.cutoff = filter.clone.cutoff
          )
        }
      }
      
      ##### generate circos plot for mice only, no hashtag information
      for (mouse.id in c("m3", "m7")){
        selected.mids <- subset(meta.data.non.splitted, meta.data.non.splitted$mouse == mouse.id)$SampleID
        input.files <- all.input.files[selected.mids]
        
        fileAliases <- to_vec(
          for (item in names(input.files)){
            sprintf("%s (%s)", item, subset(meta.data.non.splitted, meta.data.non.splitted$SampleID == item)$organ)
          }
        )
        names(fileAliases) <- names(input.files)
        saveFileName <- sprintf("%s_circos.svg", mouse.id)
        outputdir <- file.path(path.to.05.output,
                               circos.group.type,
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
            filter.clone.cutoff = filter.clone.cutoff
          )
        }
      }
    }
  }
}

