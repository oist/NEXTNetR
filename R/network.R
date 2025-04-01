#' @name network_types
#' @title Creating networks
#' 
#' @description 
#' NEXTNetR currently supports the following types of static, unweighted networks
#' 
#' * [empirical_network]
#' * [erdos_renyi_network]
#' * [fully_connected_network]
#' * [acyclic_network]
#' * [configmodel_network]
#' * [configmodel_clustered_network]
#' * [wattsstrogatz_network]
#' * [barabasialbert_network]
#' * [cubiclattice_network]
#' * [adjacencylist_network]
#'
#' the following static weighted networks
#'
#' * [weighted_adjacencylist_network]
#'
#' and following temporal networks
#'
#' * [empirical_temporalnetwork]
#' * [brownian_proximity_temporalnetwork]
#' * [sirx_temporalnetwork]
NULL

#' @title Create a Erdös-Rényi network 
#' 
#' @description
#' Creates a Erdös-Rényi network with the given size and average degree
#' 
#' @param size number of nodes
#' @param avg_degree average number of neighbour each node has
#' @returns a network object
#' 
#' @details
#' 
#' The name of this function was previously miss-spelled `erdos_reyni_network`,
#' and it is still available also under its old name.
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
erdos_renyi_network <- function(size, avg_degree) {
  nextnetR_erdos_renyi_network(as.integer(size), as.double(avg_degree))
}

#' @rdname erdos_renyi_network
#' @export
erdos_reyni_network <- erdos_renyi_network

#' @title Create a fully-connected network 
#' 
#' @description
#' Creates a fully-connected network with the given size
#' 
#' @param size number of nodes
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
fully_connected_network <- function(size) {
  nextnetR_fully_connected_network(as.integer(size))
}

#' @title Create an acyclic network 
#' 
#' @description
#' Creates an acyclic network with the given size and average degree. Node
#' degrees follow a Possonian distribution.
#' 
#' @param size number of nodes
#' @param avg_degree average number of neighbour each node has
#' @param reduced_root_degree if true, the degree of the root node is reduced by one
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
acyclic_network <- function(size, avg_degree, reduced_root_degree) {
  nextnetR_acyclic_network(as.integer(size), as.double(avg_degree), as.logical(reduced_root_degree))
}

#' @title Create an network with the specified node degrees
#' 
#' @description
#' Creates a network in nodes have the degrees specified
#' 
#' @param degrees a vector of length \eqn{N} listing the degrees of all nodes
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
configmodel_network <- function(degrees) {
  nextnetR_configmodel_network(as.integer(degrees))
}

#' @title Create an network with the specified node degrees and clustering
#' 
#' @description
#' Creates a network in nodes have the degrees specified and numbers of trianges specified
#' 
#' TODO: citation
#' 
#' @param degrees a vector of length \eqn{N} listing the degrees of all nodes
#' @param alpha_or_ck_or_triangles TODO
#' @param beta TODO
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
configmodel_clustered_network <- function(degrees, alpha_or_ck_or_triangles, beta) {
  if (is.numeric(alpha_or_ck_or_triangles) && (length(alpha_or_ck_or_triangles) == 1))
    nextnetR_configmodel_clustered_alpha_network(degrees, nextnetR_configmodel_clustered_alpha_network, beta)
  else if (is.function(alpha_or_ck_or_triangles))
    nextnetR_configmodel_clustered_ck_network(degrees, alpha_or_ck_or_triangles, beta)
  else
    nextnetR_configmodel_clustered_triangles_network(degrees, as.integers(alpha_or_ck_or_triangles), beta)
}

#' @title Create an Watts-Strogatz network
#' 
#' @description
#' Create an Watts-Strogatz  network with the given size and parameters \eqn{k}, \eqn{p}
#' 
#' @param size number of nodes in the network
#' @param k parameter k
#' @param p parameter p
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
wattsstrogatz_network <- function(size, k, p) {
  nextnetR_watts_strogatz_network(as.integer(size), as.integer(k), as.numeric(p))
}


#' @title Create an Barabasi-Albert network
#' 
#' @description
#' Create an Barabasi-Albert prefential attachment network with the given size and parameter \eqn{m}
#' 
#' @param size number of nodes in the network
#' @param m number of nodes each new node attaches to
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
barabasialbert_network <- function(size, m) {
  nextnetR_barabasialbert_network(as.integer(size), as.integer(m))
}

#' @name cubiclattice_network
#' @title Create an cubic lattice
#' 
#' @description
#' Creates a cubic lattice in the specified number of dimensions
#' 
#' @param length the number of nodes on each side of the \eqn{d}-dimensional cube
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}

#' @rdname cubiclattice_network
#' @export
cubiclattice2d_network <- function(length) {
  nextnetR_cubiclattice2d_network(length)
}

#' @rdname cubiclattice_network
#' @export
cubiclattice3d_network <- function(length) {
  nextnetR_cubiclattice3d_network(length)
}

#' @rdname cubiclattice_network
#' @export
cubiclattice4d_network <- function(length) {
  nextnetR_cubiclattice4d_network(length)
}

#' @rdname cubiclattice_network
#' @export
cubiclattice5d_network <- function(length) {
  nextnetR_cubiclattice5d_network(length)
}

#' @rdname cubiclattice_network
#' @export
cubiclattice6d_network <- function(length) {
  nextnetR_cubiclattice6d_network(length)
}

#' @rdname cubiclattice_network
#' @export
cubiclattice7d_network <- function(length) {
  nextnetR_cubiclattice7d_network(length)
}

#' @rdname cubiclattice_network
#' @export
cubiclattice8d_network <- function(length) {
  nextnetR_cubiclattice8d_network(length)
}

#' @title Creates a network for an adjacency list stored in a file
#' 
#' @description
#' The file must contain one line per node listing that node's neighbours
#' separated by commas. Nodes are referred to by their line numbers which
#' start with zero.
#' 
#' @param filename name of the file
#' @returns a network object
#' 
#' @export
empirical_network <- function(filename) {
  nextnetR_empirical_network(as.character(filename))
}

#' @title Create a network from an adjacency list
#' 
#' @description
#' Create an network object from an adjacencylist
#'
#' @param adjacencylist a list of vectors containing the neighbours of each node. Same format as the second entry in the return value of \code{\link{network_adjacencylist}}.
#' @param is_undirected true if the network is supposed to be undirected, i.e. contains a link from \eqn{i} to \eqn{j} exactly if it contains a link from \eqn{j} to \eqn{i}
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
adjacencylist_network <- function(adjacencylist, is_undirected = FALSE) {
  nextnetR_adjacencylist_network(lapply(adjacencylist, as.integer), as.logical(is_undirected))
}

#' @title Create a network from an adjacency list
#' 
#' @description
#' Create an network object from an adjacency list
#'
#' @param adjacencylist a list of vectors containing the neighbours of each node. Same format as the second entry in the return value of \code{\link{weighted_network_adjacencylist}}.
#' @param is_undirected true if the network is supposed to be undirected, i.e. contains a link from \eqn{i} to \eqn{j} exactly if it contains a link from \eqn{j} to \eqn{i}
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
weighted_adjacencylist_network <- function(adjacencylist, is_undirected = FALSE) {
  nextnetR_weighted_adjacencylist_network(lapply(adjacencylist, function(e) list(n=as.integer(e$n), w=as.numeric(e$w))),
                                          as.logical(is_undirected))
}

#' @title Create a Brownian proximity network
#' 
#' @description
#' Create a Brownian proximity network. `size` nodes are places randomly in two
#' dimensions, and connected by links if their distance does not exceed `radius`.
#' The size of the playing field is chosen such that the expected number of neighbours of
#' each node is `avg_degree`.
#' 
#' The resulting network is *temporal*, i.e. evolves over time. Non-infected nodes diffuse
#' with diffusivity \eqn{D_0}, infected nodes with diffusivity \eqn{D_1}. As the number of
#' infected nodes grows, diffusivities are additionally scaled by the factor \eqn{(1 - N_{inf} / N)^{\gamma}}.
#' 
#' The network evolves in descreet time steps of length `dt`. If left unspecified, a suitable
#' `dt` is chosen automatically based on the diffusivities and contact radius.
#'
#' @param size the number of nodes in the network
#' @param avg_degree the average number of neighbours of each node
#' @param radius the contact radius, i.e. the distance below which two nodes are connected by links
#' @param D0 the diffusivity of non-infected nodes
#' @param D1 the diffusivity of infected nodes
#' @param gamma the exponenent which with diffusivities are scaled as the epidemic grows
#' @param dt the time step used when evolving the network
#' @returns a network object
#' 
#' @seealso \code{\link{network_coordinates}}, \code{\link{network_bounds}}
#' 
#' @export
brownian_proximity_temporalnetwork <- function(size, avg_degree, radius, D0, D1=NULL,
                                        gamma=0.0, dt=NULL)
{
  if (is.null(D1))
    D1 <- D0  
  nextnetR_brownian_proximity_temporalnetwork(as.integer(size), as.double(avg_degree), as.double(radius),
                                              as.double(D0), as.double(D1), as.double(gamma),
                                              if (!is.null(dt)) dt else NULL)
}

#' TODO
#' 
#' @export
sirx_temporalnetwork <- function(graph, kappa0, kappa) {
  nextnetR_sirx_temporalnetwork(graph, as.numeric(kappa0), as.numeric(kappa))
}

#' TODO
#' 
#' @export
empirical_temporalnetwork <- function(file, finite_duration, dt) {
  nextnetR_empirical_temporalnetwork(as.character(file), as.logical(finite_duration), as.double(dt))
}

#' @name network_properties
#' @title Querying properties of networks
#' 
#' @description These functions allow properties and the topology of networks to be inspected
#'
#' @seealso \code{\link{network_types}}
#' 
#' @param nw a network object
#' @param nodes vector of node indices
#' @param indices vector of neighbour indices
#' 
#' @returns
#' * `network_is_undirected(nw)`
#'   returns true if the network is not directed, i.e. if there is a link from node \eqn{i}
#'   to \eqn{j} exactly if there is a link from node \eqn{j} to \eqn{i}.
#'   
#' * `network_size(nw)`
#'   returns the number of nodes in the network
#'   
#' * `network_outdegree(nw, node)`
#'   returns a vector of the same length as `nodes` containing  the out-degrees (i.e. number
#'   of outgoing links) of the nodes in `nodes`.
#'   
#' * `network_neighbour(nw, nodes, indices)`
#'   returns a vector of the same length as `nodes` and `indices` containing the neighbours
#'   with the given index of the given nodes.
#'   
#' * `network_adjacencylist(nw)`
#'   returns a named list which contains two entries, *nodes* and *neighbours*.
#'   *nodes* contains the indices of all nodes in the network, i.e. `1:network_size(nw)`.
#'   *neighbours* is a list of vectors, where *neighbours\[i\]* lists the neighbours of
#'   node \eqn{i}.
#'   
#' * `weighted_network_adjacencylist(nw)`
#'   returns a named list which contains two entries, *nodes* and *neighbours*.
#'   *nodes* contains the indices of all nodes in the network, i.e. `1:network_size(nw)`.
#'   *neighbours* is a list of named two-element lists containing vectores "n" and "w". The vector "n"
#'   in *neighbours\[i\]* contains the neighbours of node \eqn{i}, and the vector "w" contains the
#'   corresponding weights.
#'   
#' * `network_bounds(nw)`
#'   for networks embedded into \eqn{d}-dimensional space, this function returns a list containing
#'   two vectors of length \eqn{d}, \eqn{x} and \eqn{y}. \eqn{x} is a lower-bound and \eqn{y} the
#'   upper bound for the coordinates of the nodes in the network.
#'   
#' * `network_coordinates(nw, nodes)`
#'   for networks embedded into \eqn{d}-dimensional space, this function returns a \eqn{n\times d}{n x d}
#'   matrix containing the coordinates of the \eqn{n} nodes listed in `nodes`.
NULL

#' @rdname network_properties
#' @export
network_is_undirected <- function(nw) {
  nextnetR_network_is_undirected(nw)
}

#' @rdname network_properties
#' @export
network_size <- function(nw) {
  nextnetR_network_size(nw)
}

#' @rdname network_properties
#' @export
network_outdegree <- function(nw, nodes) {
  nextnetR_network_outdegree(nw, as.integer(nodes))
}

#' @rdname network_properties
#' @export
network_neighbour <- function(nw, nodes, indices) {
  nextnetR_network_neighbour(nw, as.integer(nodes), as.integer(indices))
}

#' @rdname network_properties
#' @export
network_adjacencylist <- function(nw) {
  nextnetR_network_adjacencylist(nw)
}

#' @rdname network_properties
#' @export
weighted_network_adjacencylist <- function(nw) {
  nextnetR_weighted_network_adjacencylist(nw)
}


#' @rdname network_properties
#' @export
network_bounds <- function(nw) {
  nextnetR_network_bounds(nw)
}

#' @rdname network_properties
#' @export
network_coordinates <- function(nw, nodes) {
  nextnetR_network_coordinates(nw, as.integer(nodes))
}

#' Computes the reproduction matrix \eqn{M_kk}
#' @export
network_reproduction_matrix <- function(nw) {
  nextnetR_reproduction_matrix(nw)
}
