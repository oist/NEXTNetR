% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{simulation_functions}
\alias{simulation_functions}
\alias{simulation_transmissiontime}
\alias{simulation_resettime}
\alias{simulation_network}
\alias{simulation_graph}
\alias{simulation_options}
\alias{simulation_isinfected}
\alias{simulation_ninfected}
\alias{simulation_addinfections}
\title{Running simulations and querying their state and properties}
\usage{
simulation_transmissiontime(sim)

simulation_resettime(sim)

simulation_network(sim)

simulation_graph(sim)

simulation_options(sim)

simulation_isinfected(sim, nodes)

simulation_ninfected(sim)

simulation_addinfections(sim, nodes, times)
}
\arguments{
\item{sim}{a simulation object}
}
\value{
\itemize{
\item \code{simulation_transmissiontime(sim)}
return the \link[=time_properties]{transmission time} distribution
\item \code{simulation_resettime(sim)}
return the \link[=time_properties]{reset/recovery time} distribution
\item \code{simulation_network(sim)}
return the \link[=network_properties]{network}. Previusly called \code{simulation_graph}.
\item \code{simulation_options(sim)}
return the algorithm options specified when the simulation was created
\item \code{simulation_isinfected(sim, nodes)}
returns a boolen vector (of same length as \code{nodes}) containing true if the node is infected
\item \code{simulation_ninfected(sim)}
returns the current number of infected nodes
\item \code{simulation_addinfections(sim, nodes, times)}
markes the nodes in \code{nodes} as infected at the specific times in \code{times}
\item \code{simulation_run(sim, stop, opts)}
runs the simulation, see \code{\link{simulation_run}} for details
}
}
\description{
The functions allow simulations created with \code{\link{simulation}}
to be run and their state to be queried.
}
\seealso{
\code{\link{simulation}}
}
