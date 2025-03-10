#' @name simulation_types
#' @title Creating simulations
#'
#' TODO: Writeme
#'
#' [nextreaction_simulation]
#' [nextreaction_simulation_meanfield]
#' [nmga_simulation]
NULL

#' Create a simulator using the NextReaction algorithm
#'
#' @param nw network to simulate on
#' @param psi infection [time distribution][time_distributions]
#' @param rhi reset/recivery [time distribution][time_distributions]
#' @param options named list of algorithm options, see below
#'
#' Possible options are
#'     shuffle_neighbours: Whether to shuffle the neighbours upon infecting a node. Default TRUE.
#'     edges_concurrent: Whether to activate all outgoing edges simultaenously or sequentially. If set to true, neighbours are implicitly shuffled and *shuffle_neighbours* thus has no effect. Default FALSE.
#'
#' @seealso simulation_properties
#' @export
nextreaction_simulation <- function(nw, psi, rho = NULL, options = list()) {
  nextnetR_nextreaction_simulation(nw, psi, rho, options)
}

#' @export
nextreaction_simulation_meanfield <- function(N, R0, psi, rho = NULL, options = list()) {
  nextnetR_nextreaction_simulation_meanfield(as.integer(N), as.double(R0), psi, rho, options)
}

#' Create a simulator using the non-Markovian Gillespie (nMGA) algorithm
#'
#' Possible options are
#'     approx_threshold: Threshold for infected nodes at which the approximate algorithm is used
#'     max_dt: Maximum timestep allowed for the approximate algorithm
#'     tauprec: Numerical precision used to invert the CDF in the exact algorithm
#' @export
nmga_simulation <- function(nw, psi, rho = NULL, options = list()) {
  nextnetR_nmga_simulation(nw, psi, rho, options)
}

#' @name simulation_functions
#' @title Running simulations and querying their state and properties
#' 
#' @description The functions allow [simulations][simulation_types] to be run and their state
#'   properties to be queried. See [simulations][simulation_types] for how to create a simulation.
#' @seealso nextreaction_simulation, nmga_simulation
#' 
#' @param sim a simulation object
#' 
#' @returns
#' * `simulation_transmissiontime(sim)`
#'   return the [transmission time][time_properties] distribution
#'   
#' * `simulation_resettime(sim)`
#'   return the [reset/recovery time][time_properties] distribution
#'   
#' * `simulation_network(sim)`
#'   return the [network][network_properties]
#'   
#' * `simulation_options(sim)`
#'   return the algorithm options specified when the simulation was created
#'   
#' * `simulation_isinfected(sim, nodes)`
#'   returns a boolen vector (of same length as `nodes`) containing true if the node is infected
#'   
#' * `simulation_ninfected(sim)`
#'   returns the current number of infected nodes
#'   
#' * `simulation_addinfections(sim, nodes, times)`
#'   markes the nodes in `nodes` as infected at the specific times in `times`
#'   
#' * `simulation_run(sim, stop, opts=list())`
#    runs 
NULL

#' @rdname simulation_functions
#' @export
simulation_transmissiontime <- function(sim) {
  nextnetR_simulation_transmissiontime(sim)
}

#' @rdname simulation_functions
#' @export
simulation_resettime <- function(sim) {
  nextnetR_simulation_resettime(sim)
}

#' @rdname simulation_functions
#' @export
simulation_graph <- function(sim) {
  nextnetR_simulation_graph(sim)
}

#' @rdname simulation_functions
#' @export
simulation_options <- function(sim) {
  nextnetR_simulation_options(sim)
}

#' @rdname simulation_functions
#' @export
simulation_isinfected <- function(sim, nodes) {
  nextnetR_simulation_isinfected(sim, as.integer(nodes))
}

#' @rdname simulation_functions
#' @export
simulation_ninfected <- function(sim) {
  nextnetR_simulation_ninfected(sim)
}

#' @rdname simulation_functions
#' @export
simulation_addinfections <- function(sim, nodes, times) {
  nextnetR_simulation_addinfections(sim, as.integer(nodes), as.double(times))
}

#' @rdname simulation_functions
#' @export
simulation_run <- function(sim_, stop, opts=list()) {
  nextnetR_simulation_run(sim_, as.list(stop), as.list(opts))
}

#' Outdated version of `simulation_run`
#
# @description This is a compatibility wrapper for `simulation_run`.
#
#' @rdname simulation_functions
#' @export
simulation_step <- function(sim_, steps, opts=list()) {
  nextnetR_simulation_run(sim_, list(epidemic_steps=as.integer(steps)), as.list(opts))
}
