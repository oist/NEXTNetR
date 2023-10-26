#' Create a simulator using the NextReaction algorithm
#'
#' Possible options are
#'     shuffle_neighbours: Whether to shuffle the neighbours upon infecting a node. Default TRUE.
#'     edges_concurrent: Whether to activate all outgoing edges simultaenously or sequentially. If set to true, neighbours are implicitly shuffled and *shuffle_neighbours* thus has no effect. Default FALSE.
#' @export
nextreaction_simulation <- function(nw, psi, rho = NULL, options = list()) {
  episimR_nextreaction_simulation(nw, psi, rho, options)
}

#' @export
nextreaction_simulation_meanfield <- function(N, R0, psi, rho = NULL, options = list()) {
  episimR_nextreaction_simulation_meanfield(as.integer(N), as.double(R0), psi, rho, options)
}

#' Create a simulator using the non-Markovian Gillespie (nMGA) algorithm
#'
#' Possible options are
#'     approx_threshold: Threshold for infected nodes at which the approximate algorithm is used
#'     max_dt: Maximum timestep allowed for the approximate algorithm
#'     tauprec: Numerical precision used to invert the CDF in the exact algorithm
#' @export
nmga_simulation <- function(nw, psi, rho = NULL, options = list()) {
  episimR_nmga_simulation(nw, psi, rho, options)
}

#' @export
simulation_transmissiontime <- function(sim) {
  episimR_simulation_transmissiontime(sim)
}

#' @export
simulation_resettime <- function(sim) {
  episimR_simulation_resettime(sim)
}

#' @export
simulation_graph <- function(sim) {
  episimR_simulation_graph(sim)
}

#' @export
simulation_options <- function(sim) {
  episimR_simulation_options(sim)
}

#' @export
simulation_isinfected <- function(sim, nodes) {
  episimR_simulation_isinfected(sim, as.integer(nodes))
}

#' @export
simulation_ninfected <- function(sim) {
  episimR_simulation_ninfected(sim)
}

#' @export
simulation_addinfections <- function(sim, nodes, times) {
  episimR_simulation_addinfections(sim, as.integer(nodes), as.double(times))
}

#' @export
simulation_step <- function(sim_, steps, opts=list()) {
  episimR_simulation_run(sim_, list(epidemic_steps=as.integer(steps)), as.list(opts))
}

#' @export
simulation_run <- function(sim_, stop, opts=list()) {
  episimR_simulation_run(sim_, as.list(stop), as.list(opts))
}
