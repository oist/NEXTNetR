# Introduction

`{NEXTNetR}` (**N**ext-reaction-based **E**pidemics e**X**tended to **T**emporal **Net**works) is an R package for the efficient simulation of epidemics on complex networks (including weighted and temporal networks) with arbitrary transmission and recovery time distributions. *NEXTNetR* is an R wrapper around the C++ library [NEXTNet](https://github.com/oist/NEXTNet).

# Installation

The latest version of *NEXTNetR* can be installed directly from Github by executing the following in R:

    install.packages("remotes")
    remotes::install_github("oist/NEXTNetR")

# Synopsis

The following minimal example simulated an epidemic on an Erdős–Rényi network with lognormally distributed transmission time

	sim <- simulation(
		erdos_renyi_network(1e5, 5),
		lognormal_time(6, 30, 0.1))
	simulation_addinfections(sim, 1, 0.0)
	r <- simulation_run(sim, stop=list(total_infected=300e3))
	plot(r$time, r$infected, type='l')

See [Getting Started](articles/NEXTNetR.html) for stey-by-step instructions on how to use `{NEXTNetR}`.

# Supported Features

`{NEXTNetR}` offers a range of common types of artifical networks such as `erdos_renyi()`, see `help(network_types)` for a full list. `adjacencylist_network()` and `adjacencylist_weightednetwork()` allow arbitrary unweighted and weighted networks to be used for simulations. Transmission and recovery times can likewise be arbitrary probability distributions, see `help(time_distributions)`. 

`{NEXTNetR}` also allows simulations on *temporal* networks, i.e. networks which change over time, possibly in response to epidemic events. Amongst the temporal networks currently supported by `{NEXTNetR}` are `empirical_contact_temporalnetwork()` and `activity_driven_temporalnetwork()`, see `help(network_types)` for a full list.