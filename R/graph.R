#' @export
graph_size <- function(nw) {
  episimR_graph_size(nw)
}

#' @export
graph_is_undirected <- function(nw) {
  episimR_graph_is_undirected(nw)
}

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
graph_coordinates <- function(nw, nodes) {
  episimR_graph_coordinates(nw, as.integer(nodes))
}

#' @export
graph_bounds <- function(nw) {
  episimR_graph_bounds(nw)
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
configmodel_clustered_graph <- function(degrees, alpha_or_ck_or_triangles, beta) {
  if (is.numeric(alpha_or_ck_or_triangles) && (length(alpha_or_ck_or_triangles) == 1))
    episimR_configmodel_clustered_alpha_graph(degrees, episimR_configmodel_clustered_alpha_graph, beta)
  else if (is.function(alpha_or_ck_or_triangles))
    episimR_configmodel_clustered_ck_graph(degrees, alpha_or_ck_or_triangles, beta)
  else
    episimR_configmodel_clustered_triangles_graph(degrees, as.integers(alpha_or_ck_or_triangles), beta)
}

#' @export
barabasialbert_graph <- function(size, m) {
  episimR_barabasialbert_graph(as.integer(size), as.integer(m))
}

#' @export
cubiclattice2d_graph <- function(length) {
  episimR_cubiclattice2d_graph(length)
}

#' @export
cubiclattice3d_graph <- function(length) {
  episimR_cubiclattice3d_graph(length)
}

#' @export
cubiclattice4d_graph <- function(length) {
  episimR_cubiclattice4d_graph(length)
}

#' @export
cubiclattice5d_graph <- function(length) {
  episimR_cubiclattice5d_graph(length)
}

#' @export
cubiclattice6d_graph <- function(length) {
  episimR_cubiclattice6d_graph(length)
}

#' @export
cubiclattice7d_graph <- function(length) {
  episimR_cubiclattice7d_graph(length)
}

#' @export
cubiclattice8d_graph <- function(length) {
  episimR_cubiclattice8d_graph(length)
}

#' @export
brownian_proximity_dyngraph <- function(size, avg_degree, radius, D, dt=NULL) {
  episimR_brownian_proximity_dyngraph(as.integer(size), as.double(avg_degree), as.double(radius),
                                      as.double(D), if (!is.null(dt)) dt else NULL)
}

#' @export
empirical_dyngraph <- function(file, finite_duration, dt) {
  episimR_empirical_dyngraph(as.character(file), as.logical(finite_duration), as.double(dt))
}

#' @export
stored_graph <- function(filename) {
  episimR_stored_graph(as.character(filename))
}

#' @export
userdefined_graph <- function(adjacencylist, is_undirected = FALSE) {
  episimR_userdefined_graph(lapply(adjacencylist, as.integer), as.logical(is_undirected))
}
