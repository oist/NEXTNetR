#' @export
nextreaction_simulation <- function(nw, psi, rho = NULL, options = list()) {
  episimR_nextreaction_simulation(nw, psi, rho, options)
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
simulation_addinfections <- function(sim, nodes, times) {
  episimR_simulation_addinfections(sim, as.integer(nodes), as.double(times))
}

#' @export
simulation_step <- function(sim_, steps) {
  episimR_simulation_step(sim_, as.integer(steps))
}
