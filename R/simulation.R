#' @title Create a simulation
#' 
#' @description Simulatiojn objects represent a single simulation on a
#' specific network and using specific tranmission and reset/recovery times.
#' The network can be *weighted* and/or *temporal*, see [network_types]. To
#' force a temporal network to be treated as a static network, set option
#' `static_network` to `TRUE`.
#'
#' @param nw network to simulate on
#' @param psi infection [time distribution][time_distributions]
#' @param rhi reset/recovery [time distribution][time_distributions]
#' @param options named list of algorithm options, see below
#' @param method the simulation method to use, currently only 'nextreaction'
#' @return a simulation object. See \code{\link{simulation_run}} for how to
#' run the simulation and [simulation_functions] for other functions that
#' operate on simulations.
#'
#' @details
#'
#' Possible options are
#' * *SIR*: Do not make recovered nodes susceptible again as in the SIR model. Default FALSE.
#' * *static_network*: Treat network as static, even if it is a temporal network. Default FALSE.
#' * *shuffle_neighbours*: Shuffle the neighbours upon infecting a node. Default TRUE.
#' * *edges_concurrent*: Activate all outgoing edges simultaenously or sequentially. If set to true, neighbours are implicitly shuffled and *shuffle_neighbours* thus has no effect. Default FALSE.
#'
#' @seealso \code{\link{simulation_functions}}
#' 
#' @example examples/basic_simulation.R
#' 
#' @md
#' @export
simulation <- function(nw, psi, rho = NULL, options = list(), method='nextreaction') {
  nextnetR_nextreaction_simulation(nw, psi, rho, options)
}

#' @name simulation_functions
#' @title Running simulations and querying their state and properties
#' 
#' @description The functions allow simulations created with \code{\link{simulation}}
#' to be run and their state to be queried.
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
#'   return the [network][network_properties]. Previusly called `simulation_graph`.
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
#' * `simulation_run(sim, stop, opts)`
#'   runs the simulation, see \code{\link{simulation_run}} for details
#'   
#' @seealso \code{\link{simulation}}
#'   
#' @md
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
simulation_network <- function(sim) {
  nextnetR_simulation_network(sim)
}

#' @rdname simulation_functions
#' @export
simulation_graph <- simulation_network

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

#' @title Runs a simulation until a stopping condition occurs
#' 
#' @description Runs the specified simulation until a stopping condition occurs
#' and returns a table listing all events that occurred. `simulation_run` can be
#' called repeatedly and will continue from the point where the last run stopped.
#' Note that for an epidemic to start, *initially infected* nodes must be defined
#' with \code{\link{simulation_addinfections}}.
#' 
#' @param sim a simulation object 
#' @param stop a list of stopping conditions
#' @param opts return value options
#' @return a `data.frame` containing the columns
#' * *time*: time since the start of the simulation
#' * *epidemic_step*: total number of epidemic events so far
#' * *network_step*: total number of network events so far
#' * *kind*: kind of event (*outside_infection*, *infection*, *reset*,
#'           *neighbour_added*, *neibhour_removed*, *instantenous_contact*)
#' * *node*: node affected by the event
#' * *neighbour*: neighbour (for *infection*, *neighbour_added*, *neibhour_removed*)
#' * *total_infected*: total number of infected nodes since simulations start
#' * *total_reset*: total number of recovered/reset nodes since simulation start
#' * *infected*: current number of infected nodes
#' 
#' @details 
#' 
#' Stopping conditions are specified as a named list of threshold values for
#' any combination of the following variables:
#' * *epidemic_steps*: Number of epidemic events (infection, outside_infection,
#'                     reset) since the simulation was started.
#' * *network_steps*:  Number of network events (edge added, edge removed,
#'                     instantaneous contact) since the simulation was started.
#'                     Only relevant for simulations on dynamic networks. 
#' * *time*:           time since the simulation was started
#' * *infected*:       Number of currently infected nodes
#' * *total_infected*: Total number of infections since the start of the simulation
#' * *total_reset*:    Total number of recoveries/resets since the start of the simulation
#'     
#' Possible return value options are:
#' * *network_events*: Include network events in the returned `data.frame`. Default *false*
#' * *epidemic_events*: Include epidemic events in the returned `data.frame`. Default *true*
#'     
#' @example examples/basic_simulation.R
#'
#' @seealso \code{\link{simulation_functions}}, \code{\link{simulation}},
#'
#' @md
#' @export
simulation_run <- function(sim, stop, opts=list()) {
  nextnetR_simulation_run(sim, as.list(stop), as.list(opts))
}

#' @title Outdated version of `simulation_run`
#'
#' @description This is a compatibility wrapper for `simulation_run`.
#'
#' @seealso \code{\link{simulation_run}}, \code{\link{simulation_functions}}
#'
#' @export
simulation_step <- function(sim_, steps, opts=list()) {
  nextnetR_simulation_run(sim_, list(epidemic_steps=as.integer(steps)), as.list(opts))
}
