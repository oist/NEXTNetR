% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{nmga_simulation}
\alias{nmga_simulation}
\title{Create a simulator using the non-Markovian Gillespie (nMGA) algorithm}
\usage{
nmga_simulation(nw, psi, rho = NULL, options = list())
}
\arguments{
\item{nw}{network to simulate on}

\item{psi}{infection \link[=time_distributions]{time distribution}}

\item{options}{named list of algorithm options, see below}

\item{rhi}{reset/recovery \link[=time_distributions]{time distribution}}
}
\value{
a simulation object. See \code{\link{simulation_run}} for how to
run the simulation and \link{simulation_functions} for other functions that
operate on \link[=simulation_types]{simulations}.
}
\description{
This is an alternative simulation algorithm that is generally
not recommended except for comparisons with \code{\link{nextreaction_simulation}}.
Only supports static, unweighted networks.
}
\details{
Possible options are
\itemize{
\item \emph{approx_threshold}: Threshold for infected nodes at which the approximate algorithm is used
\item \emph{max_dt}: Maximum timestep allowed for the approximate algorithm
\item \emph{tauprec}: Numerical precision used to invert the CDF in the exact algorithm
}
}
\seealso{
\code{\link{simulation_functions}}
}
