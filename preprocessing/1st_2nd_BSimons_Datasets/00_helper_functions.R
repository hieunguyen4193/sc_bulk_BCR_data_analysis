`%ni%` = Negate(`%in%`)

#####----------------------------------------------------------------------#####
# function to create an interactive data table
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

#####----------------------------------------------------------------------#####
# convert a gene from all upper-case to the format in our data
#####----------------------------------------------------------------------#####

to_lower_gene_symbol <- function(s){
  list_string <- unlist(strsplit(s, ""))
  
  list_string[1] <- toupper(list_string[1])
  list_string[2:length(list_string)] <- tolower(list_string[2:length(list_string)])
  new.string <- paste(list_string, collapse = "")
  return(new.string)
}


#####----------------------------------------------------------------------#####
# Function to remove any clusters from any input Seurat R object and run the 
# downstream analysis. 
#####----------------------------------------------------------------------#####

remove_clusters_from_an_object <- function(s.obj, 
                                           clusters.to.be.removed, 
                                           with.re.integration, 
                                           save.dataset.name,
                                           path.to.output, 
                                           cluster.resolution = 1,
                                           generate.html = FALSE){
  if (with.re.integration == TRUE){
    status <- "with_reInt"
  } else if (with.re.integration == FALSE){
    status <- "without_reInt"
  } else {
    stop("The parameter with.re.integration must be either TRUE or FALSE")
  }
  
  save.obj.name <- sprintf("%s_removed_%s.%s.res%s.rds", save.dataset.name, paste(clusters.to.be.removed, collapse = "_"), status, cluster.resolution)
  save.html.name <- str_replace(save.obj.name, ".rds", ".html")
  
  path.to.output <- file.path(path.to.output, str_replace(save.obj.name, ".rds", ""))
  dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
  
  # Fixed parameters
  chosen.seed <- 42
  num.dim.integration <- 25 
  num.PCA <- 25
  num.dim.cluster <- 25
  num.PC.used.in.Clustering <- 25
  
  #####----------------------------------------------------------------------#####
  # MAIN SCRIPTS
  #####----------------------------------------------------------------------#####
  
  ##### 1. Processing the Seurat object
  s.obj.removed <- subset(s.obj, seurat_clusters %in% clusters.to.be.removed == FALSE)
  DefaultAssay(s.obj.removed) <- "RNA"
  
  if (with.re.integration == TRUE){
    
    data.list <- SplitObject(s.obj.removed, split.by = "name")
    data.list <- lapply(X = data.list, FUN = function(x) {
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)})
    
    k.filter <- 200
    
    anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:num.dim.integration, scale=F,
                                      k.filter = k.filter) ## THIS IS CCA DIMENSIONS
    
    s.obj_inte <- IntegrateData(anchorset = anchors, dims = 1:num.dim.integration, k.weight = k.filter) ## THIS IS PCA DIMENSION
    
    s.obj_inte <- s.obj_inte[, colnames(s.obj.removed)]
    
    s.obj.removed[['integrated']] <- s.obj_inte[['integrated']]
    
    s.obj.removed@commands <- c(s.obj.removed@commands, s.obj_inte@commands)
    
    s.obj.removed@tools <- c(s.obj.removed@tools, s.obj_inte@tools)
    
    DefaultAssay(s.obj.removed) <- "integrated"
    
    remove.YFP.VariableFeatures <- to_vec(for (item in VariableFeatures(s.obj.removed)) if (item != "YFP") item)
    
    s.obj.removed <- ScaleData(s.obj.removed, verbose = FALSE, features = remove.YFP.VariableFeatures)
    
    s.obj.removed <- RunPCA(s.obj.removed, npcs = num.PCA, verbose = FALSE, reduction.name="INTE_PCA")
    
    s.obj.removed <- RunUMAP(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.PCA, reduction.name="INTE_UMAP")
    
    s.obj.removed <- FindNeighbors(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.dim.cluster)
    
    s.obj.removed <- FindClusters(s.obj.removed, resolution = cluster.resolution)
    
  } else {
    
    s.obj.removed <- NormalizeData(s.obj.removed) # ---> use Log Normalized
    s.obj.removed <- FindVariableFeatures(s.obj.removed, selection.method = "vst")
    
    remove.YFP.VariableFeatures <- to_vec(for (item in VariableFeatures(s.obj.removed)) if (item != "YFP") item)
    
    s.obj.removed <- ScaleData(s.obj.removed, features = to_vec(for (item in rownames(s.obj.removed)) if(item != "YFP") item))
    
    s.obj.removed <- RunPCA(s.obj.removed, npcs = num.PCA, verbose = FALSE, reduction.name="INTE_PCA", features = remove.YFP.VariableFeatures)
    
    s.obj.removed <- RunUMAP(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.PCA, reduction.name="INTE_UMAP", seed.use = chosen.seed)
    s.obj.removed <- FindNeighbors(s.obj.removed, reduction = "INTE_PCA", dims = 1:num.PC.used.in.Clustering)
    s.obj.removed <- FindClusters(s.obj.removed, resolution = cluster.resolution, random.seed = chosen.seed)
  }
  
  saveRDS(s.obj.removed, file.path(path.to.output, save.obj.name))
  
  if (generate.html == TRUE){
    ##### 2. Generate report from defined templates
    rmarkdown::render(input = file.path(path.to.template, "01_generate_cluster_DE_genes.Rmd"), 
                      params = list(
                        path.to.input.sobj = file.path(path.to.output, save.obj.name),
                        path.to.wd = wd,
                        dataset_name = save.obj.name,
                        path.to.output = path.to.output
                      ),
                      output_file = save.html.name,
                      output_dir = path.to.html)  
  }
  
}


calculate_shannon_entropy <- function(barcodes, input.s.obj){
  cell.cluster.data <- input.s.obj@meta.data %>%
    rownames_to_column("barcode") %>%
    subset(select = c(barcode, seurat_clusters)) %>%
    subset(barcode %in% barcodes)
  
  row.names(cell.cluster.data) <- NULL
  
  N <- length(unique(input.s.obj$seurat_clusters))
  
  count.cell.in.cluster <- cell.cluster.data %>% column_to_rownames("barcode") %>% table() %>% as.data.frame() 
  count.cell.in.cluster <- count.cell.in.cluster %>%
    rowwise() %>%
    mutate(Prob = Freq / sum(count.cell.in.cluster$Freq)) 
  
  count.cell.in.cluster <- subset(count.cell.in.cluster, count.cell.in.cluster$Prob != 0)
  
  shannon_entropy <- -sum(count.cell.in.cluster$Prob * log2(count.cell.in.cluster$Prob))/log2(N)
  
  return(shannon_entropy)
}

#####----------------------------------------------------------------------#####
##### RUN MONOCLE2 FROM S.OBJ
#####----------------------------------------------------------------------#####
run_monocle2 <- function(s.obj, path.to.save.monocle.obj){
  my_random_seed <- 411
  library(monocle)
  data <- GetAssayData(s.obj, slot = "count", assay = "RNA")
  
  pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
  
  fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fd)
  
  monocle.obj <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  monocle.obj <- estimateSizeFactors(monocle.obj)
  monocle.obj <- estimateDispersions(monocle.obj)
  
  monocle.obj <- detectGenes(monocle.obj, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(monocle.obj),
                                      num_cells_expressed >= 100))
  
  fData(monocle.obj)$use_for_ordering <-
    fData(monocle.obj)$num_cells_expressed > 0.05 * ncol(monocle.obj)
  
  ordering.genes <- subset(fData(monocle.obj), fData(monocle.obj)$use_for_ordering == TRUE)$gene_short_name
  ordering.genes <- intersect(ordering.genes, expressed_genes)
  monocle.obj <- monocle.obj[ordering.genes,]
  
  set.seed(my_random_seed)
  monocle.obj <- reduceDimension(monocle.obj,
                                 max_components = 2,
                                 norm_method = 'log',
                                 num_dim = 3,
                                 reduction_method = 'tSNE',
                                 verbose = T, 
                                 random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- clusterCells(monocle.obj, verbose = F, random_seed = my_random_seed)
  
  clustering_DEG_genes <-
    differentialGeneTest(monocle.obj[ordering.genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 20)
  
  HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  
  monocle.obj <-
    setOrderingFilter(monocle.obj,
                      ordering_genes = HSMM_ordering_genes)
  set.seed(my_random_seed)
  monocle.obj <-
    reduceDimension(monocle.obj, method = 'DDRTree', random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- orderCells(monocle.obj)
  
  p <- plot_cell_trajectory(monocle.obj)
  
  saveRDS(monocle.obj, file.path(path.to.save.monocle.obj, sprintf("monocle_obj.rds")))
  return(monocle.obj)
}

`%ni%` = Negate(`%in%`)

#####----------------------------------------------------------------------#####
# function to create an interactive data table
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

#####
##### RUN MONOCLE2 FROM S.OBJ
#####
run_monocle2 <- function(s.obj, path.to.save.monocle.obj){
  library(monocle)
  data <- GetAssayData(s.obj, slot = "count", assay = "RNA")
  
  pd <- new('AnnotatedDataFrame', data = s.obj@meta.data)
  
  fd <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fd)
  
  monocle.obj <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())
  monocle.obj <- estimateSizeFactors(monocle.obj)
  monocle.obj <- estimateDispersions(monocle.obj)
  
  monocle.obj <- detectGenes(monocle.obj, min_expr = 0.1)
  expressed_genes <- row.names(subset(fData(monocle.obj),
                                      num_cells_expressed >= 100))
  
  fData(monocle.obj)$use_for_ordering <-
    fData(monocle.obj)$num_cells_expressed > 0.05 * ncol(monocle.obj)
  
  ordering.genes <- subset(fData(monocle.obj), fData(monocle.obj)$use_for_ordering == TRUE)$gene_short_name
  ordering.genes <- intersect(ordering.genes, expressed_genes)
  monocle.obj <- monocle.obj[ordering.genes,]
  
  set.seed(my_random_seed)
  monocle.obj <- reduceDimension(monocle.obj,
                                 max_components = 2,
                                 norm_method = 'log',
                                 num_dim = 3,
                                 reduction_method = 'tSNE',
                                 verbose = T, 
                                 random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- clusterCells(monocle.obj, verbose = F, random_seed = my_random_seed)
  
  clustering_DEG_genes <-
    differentialGeneTest(monocle.obj[ordering.genes,],
                         fullModelFormulaStr = '~Cluster',
                         cores = 20)
  
  HSMM_ordering_genes <-
    row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]
  
  monocle.obj <-
    setOrderingFilter(monocle.obj,
                      ordering_genes = HSMM_ordering_genes)
  set.seed(my_random_seed)
  monocle.obj <-
    reduceDimension(monocle.obj, method = 'DDRTree', random_seed = my_random_seed)
  set.seed(my_random_seed)
  monocle.obj <- orderCells(monocle.obj)
  
  p <- plot_cell_trajectory(monocle.obj)
  
  saveRDS(monocle.obj, file.path(path.to.save.monocle.obj, sprintf("monocle_obj.rds")))
  return(monocle.obj)
}

#####----------------------------------------------------------------------#####
##### extract monocle information from monocle object, helper function
#####----------------------------------------------------------------------#####
extract_monocle_info <- function(cds) {
  if (cds@dim_reduce_type != "DDRTree") {
    stop(paste0("For now tradeSeq only support Monocle with DDRTree",
                "reduction. If you want to use another type",
                "please use another format for tradeSeq inputs."))
  }
  # Get the reduced dimension of DDRT
  rd <- t(monocle::reducedDimS(cds)) %>% as.data.frame()
  
  # Get the various lineages info for weights and pseudotime
  y_to_cells <- cds@auxOrderingData[["DDRTree"]]
  y_to_cells <- y_to_cells$pr_graph_cell_proj_closest_vertex %>%
    as.data.frame()
  y_to_cells$cells <- rownames(y_to_cells)
  y_to_cells$Y <- y_to_cells$V1
  root <- cds@auxOrderingData[[cds@dim_reduce_type]]$root_cell
  root <- y_to_cells$Y[y_to_cells$cells == root]
  mst <- monocle::minSpanningTree(cds)
  endpoints <- names(which(igraph::degree(mst) == 1))
  endpoints <- endpoints[endpoints != paste0("Y_", root)]
  cellWeights <- lapply(endpoints, function(endpoint) {
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    df <- y_to_cells[y_to_cells$Y %in% path, ]
    df <- data.frame(weights = as.numeric(colnames(cds) %in% df$cells))
    colnames(df) <- endpoint
    return(df)
  }) %>% do.call(what = 'cbind', args = .)
  pseudotime <- sapply(cellWeights, function(w) cds$Pseudotime)
  rownames(cellWeights) <- rownames(pseudotime) <- colnames(cds)
  # Get the lineages representation
  edges_rd <- t(monocle::reducedDimK(cds)) %>% as.data.frame()
  rd_lineages <- lapply(endpoints, function(endpoint){
    path <- igraph::shortest_paths(mst, root, endpoint)$vpath[[1]]
    path <- as.character(path)
    path <- paste("Y", path, sep = "_")
    return(edges_rd[path, ])
  })
  return(list("pseudotime" = pseudotime,
              "cellWeights" = as.matrix(cellWeights)))
}


#####----------------------------------------------------------------------#####
##### Get earliest principal node 
#####----------------------------------------------------------------------#####
get_earliest_principal_node <- function(cds, cluster.id){
  cell_ids <- which(colData(cds)[, "seurat_clusters"] == cluster.id)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}