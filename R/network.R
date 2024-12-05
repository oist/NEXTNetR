#' @export
network_size <- function(nw) {
  nextnetR_network_size(nw)
}

#' @export
network_is_undirected <- function(nw) {
  nextnetR_network_is_undirected(nw)
}

#' @export
network_outdegree <- function(nw, nodes) {
  nextnetR_network_outdegree(nw, as.integer(nodes))
}

#' @export
network_neighbour <- function(nw, nodes, indices) {
  nextnetR_network_neighbour(nw, as.integer(nodes), as.integer(indices))
}

#' @export
network_adjacencylist <- function(nw) {
  nextnetR_network_adjacencylist(nw)
}

#' @export
network_coordinates <- function(nw, nodes) {
  nextnetR_network_coordinates(nw, as.integer(nodes))
}

#' @export
network_bounds <- function(nw) {
  nextnetR_network_bounds(nw)
}

#' @export
erdos_reyni_network <- function(size, avg_degree) {
  nextnetR_erdos_reyni_network(as.integer(size), as.double(avg_degree))
}

#' @export
fully_connected_network <- function(size) {
  nextnetR_fully_connected_network(as.integer(size))
}

#' @export
acyclic_network <- function(size, avg_degree, reduced_root_degree) {
  nextnetR_acyclic_network(as.integer(size), as.double(avg_degree), as.logical(reduced_root_degree))
}

#' @export
configmodel_network <- function(degrees) {
  nextnetR_configmodel_network(as.integer(degrees))
}

#' @export
configmodel_clustered_network <- function(degrees, alpha_or_ck_or_triangles, beta) {
  if (is.numeric(alpha_or_ck_or_triangles) && (length(alpha_or_ck_or_triangles) == 1))
    nextnetR_configmodel_clustered_alpha_network(degrees, nextnetR_configmodel_clustered_alpha_network, beta)
  else if (is.function(alpha_or_ck_or_triangles))
    nextnetR_configmodel_clustered_ck_network(degrees, alpha_or_ck_or_triangles, beta)
  else
    nextnetR_configmodel_clustered_triangles_network(degrees, as.integers(alpha_or_ck_or_triangles), beta)
}

#' @export
barabasialbert_network <- function(size, m) {
  nextnetR_barabasialbert_network(as.integer(size), as.integer(m))
}

#' @export
cubiclattice2d_network <- function(length) {
  nextnetR_cubiclattice2d_network(length)
}

#' @export
cubiclattice3d_network <- function(length) {
  nextnetR_cubiclattice3d_network(length)
}

#' @export
cubiclattice4d_network <- function(length) {
  nextnetR_cubiclattice4d_network(length)
}

#' @export
cubiclattice5d_network <- function(length) {
  nextnetR_cubiclattice5d_network(length)
}

#' @export
cubiclattice6d_network <- function(length) {
  nextnetR_cubiclattice6d_network(length)
}

#' @export
cubiclattice7d_network <- function(length) {
  nextnetR_cubiclattice7d_network(length)
}

#' @export
cubiclattice8d_network <- function(length) {
  nextnetR_cubiclattice8d_network(length)
}

#' @export
brownian_proximity_temporalnetwork <- function(size, avg_degree, radius, D0, D1=NULL,
                                        gamma=0.0, dt=NULL)
{
  if (is.null(D1))
    D1 <- D0  
  nextnetR_brownian_proximity_tempornextnetR_(as.integer(size), as.double(avg_degree), as.double(radius),
                                      as.double(D0), as.double(D1), as.double(gamma),
                                      if (!is.null(dt)) dt else NULL)
}

#' @export
empirical_temporalnetwork <- function(file, finite_duration, dt) {
  nextnetR_empirical_tempornextnetR_(as.character(file), as.logical(finite_duration), as.double(dt))
}

#' @export
sirx_temporalnetwork <- function(graph, kappa0, kappa) {
  nextnetR_sirx_tempornextnetR_(graph, as.numeric(kappa0), as.numeric(kappa))
}

#' @export
empirical_network <- function(filename) {
  nextnetR_empirical_network(as.character(filename))
}

#' @export
adjacencylist_network <- function(adjacencylist, is_undirected = FALSE) {
  nextnetR_adjacencylist_network(lapply(adjacencylist, as.integer), as.logical(is_undirected))
}
