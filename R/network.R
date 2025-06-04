#' @name network_types
#' @title Creating networks
#' 
#' @description 
#' 
#' In *NEXTNet*, networks are directed graphs, i.e. fully defined by a set of
#' nodes \eqn{V} and edges (or links) \eqn{E \subset V \times V} where self-edges
#' i.e. links of the form \eqn{(v,v)} are forbidden. *NEXTNet* provides various
#' different types of algorithms for creating synthetic networks, as well as
#' networks with are defined by specifying \eqn{V} and \eqn{E}.
#' 
#' In addition to the plain (unweighted and static) networks defined above, 
#' *NEXTNet* also supports *temporal* and *weighted* networks.
#' 
#' In *temporal* networks, the set of edges/links depends on the time \eqn{t},
#' i.e. \eqn{E=E(t)}. On such networks, epidemics can only spread from a node
#' \eqn{u} to a node \eqn{v} at times where \eqn{(u,v) \in E(t)}. At other times,
#' the infectiousness of \eqn{u} for the particular link \eqn{(u,v)} is effectively
#' zero. 
#' 
#' In *weighted* networks, each \eqn{e \in E} has an assigned weight
#' \eqn{w_e \geq 0}. These weights modulate the infectiousness
#' of nodes, see the discussion in [time_distributions]. Weighted networks can
#' be interpreted as a limit case of temporal networks in which edges fluctuate
#' with at a very high frequency. The weight then expresses the fraction of times
#' at which the link is present.
#' 
#' *NEXTNetR* supports the following types of static, unweighted networks:
#' 
#' * [empirical_network]: Network defined by an arbitrary adjacency list read from a file.
#' * [adjacencylist_network] Network defined by an arbitrary adjacency list.
#' * [erdos_renyi_network]: Erdős–Rényi network, i.e. edges are sampled i.i.d from a fully-connected network.
#' * [fully_connected_network]: Fully connected network, i.e. all possible edges exist.
#' * [acyclic_network]: Tree-shaped network.
#' * [configmodel_network]: Network with specified number of nodes of a certain degree.
#' * [configmodel_clustered_network]: Configuration model with clustering.
#' * [wattsstrogatz_network]: 
#' * [barabasialbert_network]: Barabási–Albert prefertial attachment network.
#' * [cubiclattice_network]: Cubic lattice in 2 up to 8 dimensions.
#'
#' the following static weighted networks
#'
#' * [empirical_weightednetwork]: Weighted network defined by an adjacency list read from a file
#' * [adjacencylist_weightednetwork]: Network defined by an arbitrary adjacency list with weighted edges.
#' * [erdos_renyi_weightednetwork]: Erdős–Rényi network with i.i.d edge weights.
#'
#' and the following temporal networks
#'
#' * [empirical_contact_temporalnetwork]: Network defined by contacts at pre-defined times read from a file.
#' * [erdos_renyi_temporalnetwork]: Erdős–Rényi network with temporally evolving edges.
#' * [brownian_proximity_temporalnetwork]: Proximity network for Brownian particles in two dimensions.
#' * [sirx_temporalnetwork]: Network version of the SIRX model proposed by Maier & Brockmann, 2020
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
#' Creates a network in which nodes have the degrees specified and numbers of triangles
#' specified using an implementation of the algorithm by Serrano & Boguñá,
#' 2005, Phys. Rev. E 72, 036133.
#' 
#' @param degrees a vector of length \eqn{N} listing the degrees of all nodes
#' @param ck arbitrary function \eqn{c(k)} for the probability that if \eqn{b}, \eqn{c}
#'           are neighbours of \eqn{a}, then \eqn{b} and \eqn{c} are neighbours.
#' @param triangles number of triangles overlapping nodes of degree 0, 1, 2, ...,
#' @param alpha set \eqn{c(k)=0.5*(k-1)^\alpha}
#' @param beta parameter that defines degree class probabilities \eqn{P(k)}
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
configmodel_clustered_network <- function(degrees, beta, alpha=NULL, ck=NULL, triangles=NULL) {
  if (!is.null(alpha))
    nextnetR_configmodel_clustered_alpha_network(as.integer(degrees), as.numeric(alpha), as.numeric(beta))
  else if (!is.null(ck))
    nextnetR_configmodel_clustered_ck_network(as.integer(degrees), as.function(ck), as.numeric(beta))
  else if (!is.null(triangles))
    nextnetR_configmodel_clustered_triangles_network(as.integer(degrees), as.integer(triangles), as.numeric(beta))
  else
    stop("either alpha, ck or triangles must be specified")
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
#'
#' Can read networks from a file, or download networks packages for
#' NEXT-Net directly from the NEXTNet-EmpiricalNetworks repository.
#'
#' The file must contain one line per node listing first a node and then
#' that node's neighbours, separated by whitespace (by default; other
#' separated can be specified). Node are identified by numbers starting
#' with `idxbase` (by default 1); the maximal node index that appears in
#' the file defines the size of the network. For undirected networks (i.e.
#' if `undirected=TRUE`), for every link (u,v) listed in the file the reverse
#' link (v,u) is added as well.
#'
#' Lines starting with the comment chracter '#' and skipped. If the first
#' non-comment line does not start with a numerical node index, it is
#' assumed to be a header line and skipped as well.
#' 
#' @param path name of the file
#' @param name name of a packaged empirical network, see \code{\link{packaged_empirical_network}}
#' @param group packaged empirical network group containing the network
#' @param undirected if `TRUE` the network is assumed to be undirected
#' @param simplify whether to remove self-edges and multi-edges
#' @param idxbase index of the first node (typically 1 or 0, default 1)
#' @param sep separator, by default whitespace
#' @param gzip whether the file is compressed
#' @param download.timeout for packaged networks the download timeout
#' @returns a network object
#' 
#' @export
empirical_network <- function(
  path, name=NULL, group="undirected", undirected=TRUE, simplify=FALSE,
  idxbase=1, sep=' ', gzip=grepl('\\.gz$', path), download.timeout=300
) {
  if (!is.null(name)) {
      if (!missing(path) || !(sep == ' '))
          stop("when loading packaged networks, neither path nor sep is supported")
      path <- packaged_empirical_network(as.character(name), as.character(group),
                                         timeout=download.timeout)
      idxbase <- 1
      sep <- ' '
      gzip <- TRUE
      if ((group == "undirected") && (!undirected))
          warn("loading undirected network as directed")
  }
  nextnetR_empirical_network(
    path.expand(as.character(path)), as.logical(undirected),
    as.logical(simplify), as.integer(idxbase),
    as.character(sep), as.logical(gzip))
}

#' @title Creates a weighted network for an adjacency list stored in a file
#' 
#' @description
#' The file must contain one line per node listing first a node and then
#' pairs of neighbours and their weights. Pairs are separated by whitespace,
#' neighbours and weights by a colon (':'); other separators can be specifed.
#' Nodes are identified by numbers starting with `idxbase` (by default 1);
#' the maximal node index that appears in the file defines the size of the
#' network. For undirected networks (i.e. if `undirected=TRUE`), for every
#' link (u,v) listed in the file the reverse link (v,u) is added as well.
#'
#' Lines starting with the comment chracter '#' and skipped. If the first
#' non-comment line does not start with a numerical node index, it is
#' assumed to be a header line and skipped as well.
#' 
#' @param path name of the file
#' @param undirected if `TRUE` the network is assumed to be undirected
#' @param simplify whether to remove self-edges and multi-edges
#' @param idxbase index of the first node (typically 1 or 0, default 1)
#' @param csep separator between neighbours/weight pairs, by default whitespace
#' @param wsep separator between neighbour and its weight, by default ':' 
#' @param gzip whether the file is compressed
#' @returns a network object
#' 
#' @export
empirical_weightednetwork <- function(
  path, undirected=TRUE, simplify=FALSE, idxbase=1, csep=' ', wsep=':',
  gzip=grepl('\\.gz$', path)
) {
  nextnetR_empirical_weightednetwork(
    path.expand(as.character(path)), as.logical(undirected),
    as.logical(simplify), as.integer(idxbase),
    as.character(csep), as.character(wsep), as.logical(gzip))
}

#' @title Create a network from an adjacency list
#' 
#' @description
#' Create an network object from an adjacencylist
#'
#' @param adjacencylist a list of vectors containing the neighbours of each node. Same format as return value of \code{\link{network_adjacencylist}}.
#' @param is_undirected `TRUE` if the network is supposed to be undirected, i.e. contains a link from \eqn{i} to \eqn{j} exactly if it contains a link from \eqn{j} to \eqn{i}.
#' @param above_diagonal set to `TRUE` if the network is undirected and the adjacencylist only contains edges \eqn{i,j} with \eqn{i \leq j}, i.e. represent the upper triangular submatrix of the adjacency matrix. Defaults to `TRUE` for undirected networks.
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
adjacencylist_network <- function(adjacencylist, is_undirected = TRUE, above_diagonal = is_undirected) {
  nextnetR_adjacencylist_network(lapply(adjacencylist, as.integer), as.logical(is_undirected), as.logical(above_diagonal))
}

#' @title Create a network from an adjacency list
#' 
#' @description
#' Create an network object from an adjacency list
#'
#' @param adjacencylist a list of vectors containing the neighbours of each node. Same format as the return value of \code{\link{weighted_network_adjacencylist}}.
#' @param is_undirected `TRUE` if the network is supposed to be undirected, i.e. contains a link from \eqn{i} to \eqn{j} exactly if it contains a link from \eqn{j} to \eqn{i}.
#' @param above_diagonal set to `TRUE` if the network is undirected and the adjacencylist only contains edges \eqn{i,j} with \eqn{i \leq j}, i.e. represent the upper triangular submatrix of the adjacency matrix. Defaults to `TRUE` for undirected networks.
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
adjacencylist_weightednetwork <- function(adjacencylist, is_undirected = FALSE, above_diagonal = is_undirected) {
  nextnetR_adjacencylist_weightednetwork(lapply(adjacencylist, function(e) list(n=as.integer(e$n), w=as.numeric(e$w))),
                                         as.logical(is_undirected), as.logical(above_diagonal))
}

#' @title Create a weighted Erdös-Rényi network
#' 
#' @description
#' Creates a weighted Erdös-Rényi network with the given size and average degree
#' and i.i.d. edge weights. The weight distribution is specified by a vector
#' of weights and a vector of corresponding probabilities.
#' 
#' @param size number of nodes
#' @param avg_degree average number of neighbour each node has
#' @param weights edge weights
#' @param probabilities edge weight probabilities, must be same length as `weights`
#' @returns a network object
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
erdos_renyi_weightednetwork <- function(size, avg_degree, weights, probabilities) {
  stopifnot(length(weights) == length(probabilities))
  nextnetR_erdos_renyi_weightednetwork(as.integer(size), as.numeric(avg_degree),
                                       as.numeric(weights), as.numeric(probabilities))
}

#' @title Create a temporal Erdös-Rényi network
#' 
#' @description
#' Creates a temporal Erdös-Rényi network with the given size and average degree.
#' 
#' @param size number of nodes (\eqn{n})
#' @param avg_degree average number of neighbour each node has (\eqn{k})
#' @param timescale time scale of network evolution (\eqn{\tau})
#' @returns a network object
#' 
#' @details
#' At any point in time, a temporal Erdös-Rényi network resembles a static
#' Erdös-Rényi network in that each possible edge exists in the network
#' with uniform probability \eqn{p_+ = k / (n - 1)} where \eqn{k} as the average
#' node degree and \eqn{n} the number of nodes. The state of each edge possible
#' evolves independently over time, appearing with rate \eqn{\lambda_+ = p_+ / \tau}
#' disappearing with rate \eqn{\lambda_- = p_- / \tau} where \eqn{\tau} is the
#' time scale.
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
erdos_renyi_temporalnetwork <- function(size, avg_degree, timescale) {
  nextnetR_erdos_renyi_temporalnetwork(as.integer(size), as.numeric(avg_degree),
                                       as.numeric(timescale))
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
brownian_proximity_temporalnetwork <- function(size, avg_degree, radius, D0,
                                               D1=NULL, gamma=0.0, dt=NULL)
{
  if (is.null(D1))
    D1 <- D0  
  nextnetR_brownian_proximity_temporalnetwork(
    as.integer(size), as.double(avg_degree), as.double(radius),
    as.double(D0), as.double(D1), as.double(gamma),
    if (!is.null(dt)) as.double(dt) else NULL)
}

#' @title Network SIRX model
#' 
#' @description
#' This temporal network implements a network version of the SIRX model which
#' was introduced for well-mixed populations by Maier & Brockmann.
#' Starting from an arbitrary network, nodes are removed with baseline rate
#' \eqn{\kappa_0} and infected nodes are removed with elevated rate
#' \eqn{\kappa_0 + \kappa}, see Maier & Brockmann 2020, Science 368 (6492), 742-74
#' for detais.
#' 
#' @param network base network. must be simple (i.e. no self-edges or multi-edges).
#' @param kappa0 \eqn{\kappa_0}
#' @param kappa \eqn{\kappa}
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
sirx_temporalnetwork <- function(network, kappa0, kappa) {
  nextnetR_sirx_temporalnetwork(network, as.numeric(kappa0), as.numeric(kappa))
}

#' @title Empirical temporal network comprising short contacts between nodes
#' 
#' @description
#' Can read networks from a file, or download networks packages for
#' NEXT-Net directly from the NEXTNet-EmpiricalNetworks repository.
#'
#' Each line of the file must list a single contact in the form
#' "<src> <dst> <time>" where <src> and <dst> are integers specifying the nodes
#' and <time> is a floating-point value. The highest node index appearing within
#' <src> and <dst> defines the size of the network. Finite-duration contacts
#' extend by dt/2 around the time specified in the file.
#' 
#' @param path name of the file
#' @param name name of a packaged empirical network, see \code{\link{packaged_empirical_network}}
#' @param group packaged empirical network group containing the network
#' @param finite_duration whether to treat contacts as having finite duration `dt`
#'        or to be instantaneous with effective weight `weight` * `dt`
#' @param dt duration of contacts
#' @param weight weight of contacts, instantaneous contacts have effective weight `weight` * `dt`.
#' @param gzip whether the file is compressed
#' @param download.timeout for packaged networks the download timeout
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
empirical_contact_temporalnetwork <- function(path, dt, weight, finite_duration=FALSE,
                                              name=NULL, group="contact",
                                              gzip=grepl('\\.gz$', path), download.timeout=300)
{
  if (!is.null(name)) {
      if (!missing(path))
          stop("when loading packaged networks path is not supported")
      path <- packaged_empirical_network(as.character(name), as.character(group),
                                         timeout=download.timeout)
      gzip <- TRUE
  }
  nextnetR_empirical_contact_temporalnetwork(
    path.expand(as.character(path)), as.logical(finite_duration),
                as.double(dt), as.double(weight), as.logical(gzip))
}

#' @title Activity-driven network model of Cai, Nie & Holme (2024).
#' 
#' @description
#' Here, nodes are initially inactive and have degree zero. Node \eqn{i} activates with
#' rate \eqn{a[i] * \eta} and upon activation connects to \eqn{m} other uniformly chosen nodes
#' (which are not necessarily active). Active nodes inactivate with constant rate
#' \eqn{b}. See Cai, Nie & Holme 2024, Phys. Rev. Research 6, L022017 for details.
#' Here, we implement a generalized version of the model in which the
#' activation and deacivation rates of infected nodes can differ from those of
#' non-infected node.
#' 
#' @param activities node-specific activities \eqn{a_1,a_2,\ldots}
#' @param m number of nodes an activated node connects to
#' @param eta activation rate 
#' @param b deactivation rate
#' @param eta_inf activation rate for infected nodes
#' @param b_inf deactivation rate of infected nodes
#' 
#' @seealso \code{\link{network_properties}}, \code{\link{network_types}}
#' 
#' @export
activity_driven_temporalnetwork <- function(activities, m, eta, b,
                                            eta_inf = eta, b_inf = b)
{
  nextnetR_activity_driven_temporalnetwork(
    as.numeric(activities), as.integer(m),
    as.numeric(eta), as.numeric(eta_inf),
    as.numeric(b), as.numeric(b_inf))
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
#' * `network_is_simple(nw)`
#'   returns true if the network is simple, i.e. if it neither contains
#'   self-edges nor multi-edges. Note that `false` does not imply that the
#'   network necessarily contains self- or multi-edges, only that is is not
#'   guaranteed to be simple by construction.
#'   
#' * `network_is_undirected(nw)`
#'   returns true if the network is not directed, i.e. if there is a link from node \eqn{i}
#'   to \eqn{j} exactly if there is a link from node \eqn{j} to \eqn{i}. Note that
#'   similar to `network_is_simple`, a return value of `true` guarantees that the
#'   network is undirected, but `false` does not imply the existence of an edge
#'   without a reversed counterpart.
#'
#' * `network_is_weighted(nw)`
#'   return true if the network is weighted, i.e. if links have associated weights
#'   
#' * `network_is_temporal(nw)`
#'   return true if the network is temporal, i.e. if links can appear and disappear
#'   during epidemic simulations
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
#' * `network_neighbour_weight(nw, nodes, indices)`
#'   returns a named list containing vectors "n" and "w". "n" is a vector of the same
#'   length as `nodes` and `indices` and contains the neighbours with the given index of the
#'   given nodes. "w" is a vector of the same length and contains the corresponding weights.
#'   
#' * `network_adjacencylist(nw, above_diagonal=network_is_undirected(nw))`
#'   returns a list of interger vectors where \eqn{i}-th vector contains the neighbours of
#'   node \eqn{i}. For undirected networks, of the edge pair \eqn{(i,j)} and \eqn{(j, i)}
#'   only the edge \eqn{(i,j)} with \eqn{i \leq j} is output, i.e. only the upper triangular
#'   submatrix of an symmetric adjacency matrix is considered. To include both edges,
#'   set `above_diagonal=FALSE`.
#'   
#' * `weighted_network_adjacencylist(nw, above_diagonal=network_is_undirected(nw))`
#'   returns a list of named lists containing vectors "n" and "w". The vector "n"
#'   in the \eqn{i}-th list contains neighbours of node \eqn{i}, and the vector "w" contains the
#'   corresponding weights. For undirected networks, of the edge pair \eqn{(i,j)} and \eqn{(j, i)}
#'   only the edge \eqn{(i,j)} with \eqn{i \leq j} is output, i.e. only the upper triangular
#'   submatrix of an symmetric adjacency matrix is considered. To include both edges,
#'   set `above_diagonal=FALSE`.
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
network_is_simple <- function(nw) {
  nextnetR_network_is_simple(nw)
}

#' @rdname network_properties
#' @export
network_is_weighted <- function(nw) {
  nextnetR_network_is_weighted(nw)
}

#' @rdname network_properties
#' @export
network_is_temporal <- function(nw) {
  nextnetR_network_is_temporal(nw)
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
  nodes <- as.integer(nodes)
  indices <- as.integer(indices)
  if (length(nodes) == 1)
    nodes <- rep(nodes, length(indices))
  else if (length(indices) == 1)
    indices <- rep(indices, length(nodes))
  else if (length(nodes) != length(indices))
    stop("length of either nodes or indices must be one or both lengths must agree")
  nextnetR_network_neighbour(nw, nodes, indices)
}

#' @rdname network_properties
#' @export
network_neighbour_weight <- function(nw, nodes, indices) {
  nodes <- as.integer(nodes)
  indices <- as.integer(indices)
  if (length(nodes) == 1)
    nodes <- rep(nodes, length(indices))
  else if (length(indices) == 1)
    indices <- rep(indices, length(nodes))
  else if (length(nodes) != length(indices))
    stop("length of either nodes or indices must be one or both lengths must agree")
  nextnetR_network_neighbour_weight(nw, nodes, indices)
}

#' @rdname network_properties
#' @export
network_adjacencylist <- function(nw, above_diagonal = network_is_undirected(nw)) {
  nextnetR_network_adjacencylist(nw, as.logical(above_diagonal))
}

#' @rdname network_properties
#' @export
weighted_network_adjacencylist <- function(nw, above_diagonal = network_is_undirected(nw)) {
  nextnetR_weighted_network_adjacencylist(nw, as.logical(above_diagonal))
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

#' Computes the reproduction matrix \eqn{M}, reproduction number \eqn{R},
#' and other statistics related to these quantites.
#'
#' @returns a named list containing entries *M*, *r*, *k1*, *k2*, *k3*,
#'          *m_bar*, *R0*, *R_r*, *R_pert*. See *Details* for a description of
#'          these quantites
#'
#' @details Computes the following quantities
#' * *M*: Reproduction matrix, \eqn{M_{ij}} is the number of susceptible individuals of \eqn{i}-th smallest degree connected to a node of \eqn{j}-th smallest degree.
#' * *r*: Degree correlation (assortativity)
#' * *k1*: first raw moment of the degree distribution.
#' * *k2*: second raw moment of the degree distribution.
#' * *k3*: third raw moment of the degree distribution.
#' * *m1*: average number of triangles a link is part of.
#' * *m2*: second moment related to \eqn{m_1}.
#' * *R0*: basic reproduction number \eqn{R_0 = (k_2 - k_1) / k_1}
#' * *R_r*: reproduction number taking assortativity but not clustering into accounts, like \eqn{R_{pert}} but with \eqn{m_1}, \eqn{m_2} set to zero.
#' * *R_pert*: reproduction number estimated from \eqn{R_0}, \eqn{r}, \eqn{m_1}, \eqn{m_2} and \eqn{k_1}, \eqn{k_2}, \eqn{k_3}.
#'
#' See Cure, Pflug & Pigolotti, Exponential rate of epidemic spreading on complex networks, *Physical Review E*, **111**, 044311 for a detailed discussion
#' of these quantities
#'
#' @export
network_reproduction_matrix <- function(nw) {
  nextnetR_reproduction_matrix(nw)
}

#' Downloads a packages empirical network and returns the file path
#'
#' See [NEXTNetR-EmpiricalNetworks](https://github.com/oist/NEXTNet-EmpiricalNetworks)
#' for a list of available networks.
#'
#' @param name name of the network
#' @param type type of network
#' @param format file format
#' @param timeout timeout for downloading
#' @returns path to the downloaded file
#'
#' Downloaded files are cached in the directory
#' `rappdirs::user_cache_dir(appname="EmpiricalNetworks", appauthor="NEXTNetR"'
#'
#' @export
packaged_empirical_network <- function(
  name, group="undirected", format="gz", timeout=300
) {
  # Return cached copy if available. Note that we currently don't check
  # whether the copy is still current; presumably packaged networks should
  # never change
  filename <- paste0(name, ".", format)
  cache <- rappdirs::user_cache_dir(appname="NEXTNetR-EmpiricalNetworks", appauthor="NEXTNetR")
  dir.create(file.path(cache, group), recursive=TRUE, showWarnings=FALSE)
  path <- file.path(cache, group, filename)
  if (file.exists(path))
    return(path)
  
  # Construct url and filenames
  prefix <- "https://github.com/oist/NEXTNet-EmpiricalNetworks/raw/refs/heads/master/"
  url <- paste0(prefix, "/", group, "/", filename)
  path.inprogress <- paste0(path, ".in-progress")
  if (file.exists(path.inprogress))
    file.remove(path.inprogress)
 
  # Make sure the timeout is set appropriately, reset afterwards
  op <-options(timeout = max(timeout, getOption("timeout")))
  on.exit(options(op))

  # Download
  message("Downloading ", filename, " to ", file.path(cache, group))
  if (download.file(url=url, destfile=path.inprogress, quiet=FALSE) == 0) {
    file.rename(path.inprogress, path)
    return(path)
  }
  
  # Handle error
  if (file.exists(path.inprogress))
    file.remove(path.inprogress)
  stop("failed to download " + filename + " from " + url)
}