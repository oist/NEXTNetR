#' @export
graph_outdegree <- function(nw, nodes) {
  episimR_graph_outdegree(nw, as.integer(nodes))
}

#' @export
graph_neighbour <- function(nw, nodes, indices) {
  episimR_graph_neighbour(nw, as.integer(nodes), as.integer(indices))
}

#' @export
graph_adjacencylist <- function(nw) {
  episimR_graph_adjacencylist(nw)
}

#' @export
erdos_reyni_graph <- function(size, avg_degree) {
  episimR_erdos_reyni_graph(as.integer(size), as.double(avg_degree))
}

#' @export
fully_connected_graph <- function(size) {
  episimR_fully_connected_graph(as.integer(size))
}

#' @export
acyclic_graph <- function(size, avg_degree, reduced_root_degree) {
  episimR_acyclic_graph(as.integer(size), as.double(avg_degree), as.logical(reduced_root_degree))
}

#' @export
configmodel_graph <- function(degrees) {
  episimR_configmodel_graph(as.integer(degrees))
}

#' @export
scalefree_graph <- function(size) {
  episimR_scalefree_graph(as.integer(size))
}

#' @export
stored_graph <- function(filename) {
  episimR_stored_graph(as.character(filename))
}
