`%ni%` = Negate(`%in%`)

#####----------------------------------------------------------------------#####
##### function to create an interactive data table
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
##### Compute distance matrix for each pair of AA sequences
#####----------------------------------------------------------------------#####
compute_distance_matrix <- function(seqs) {
  n <- length(seqs)
  dist_matrix <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      dist <- stringdist(seqs[i], seqs[j], method = "lv")
      dist_matrix[i, j] <- dist
      dist_matrix[j, i] <- dist
    }
  }
  return(dist_matrix)
}

assign_clusters_to_sequences <- function(seqs, threshold = 0.15){
  distance_matrix <- compute_distance_matrix(seqs)
  max_length <- max(nchar(seqs))
  normalized_matrix <- distance_matrix / max_length
  graph <- graph_from_adjacency_matrix(normalized_matrix < threshold, mode = "undirected", diag = FALSE)
  cl <- cluster_optimal(graph)
  clus.res <- data.frame(seq = seqs, cluster = cl$membership)
  # function to plot the graph and its node community
  # plot(cl, graph, layout=layout_with_fr(graph))
  return(list(
    graph = graph,
    cl = cl, 
    res = clus.res
  ))
}