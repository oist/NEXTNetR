% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network.R
\name{adjacencylist_weightednetwork}
\alias{adjacencylist_weightednetwork}
\title{Create a network from an adjacency list}
\usage{
adjacencylist_weightednetwork(adjacencylist, is_undirected = FALSE)
}
\arguments{
\item{adjacencylist}{a list of vectors containing the neighbours of each node. Same format as the second entry in the return value of \code{\link{weighted_network_adjacencylist}}.}

\item{is_undirected}{true if the network is supposed to be undirected, i.e. contains a link from \eqn{i} to \eqn{j} exactly if it contains a link from \eqn{j} to \eqn{i}}
}
\value{
a network object
}
\description{
Create an network object from an adjacency list
}
\seealso{
\code{\link{network_properties}}, \code{\link{network_types}}
}
