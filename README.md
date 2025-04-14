# Introduction

*NEXTNetR* (*N*ext-reaction-based *E*pidemics e*X*tended to *T*emporal *Net*works) is an R package for the efficient simulation of epidemics on complex networks (including weighted and temporal networks) with arbitrary transmission and recovery time distributions. *NEXTNetR* is an R wrapper around the C++ library NEXTNet <https://github.com/oist/NEXTNet>.

See <https://oist.github.io/NEXTNetR/> for details.

# Installation

The latest version of *NEXTNetR* can be installed directly from Github by executing the following in R:

    install.packages(c("devtools", "git2r"))
    devtools::install_git("git@github.com:oist/NEXTNetR")

# Synopsis

The following R script sets up and runs a simulation on an Erdős–Rényi network with lognormally distributed transmission time (i.e. the generation time between subsequent infections) and lognormally distribured reset time (the time until an infected node becomes susecptible again). See [Basic Simulations of Simulations with NEXTNetR](https://oist.github.io/NEXTNetR/articles/basic_simulation.html) for a step-by step walk-through of this code.

	# Create contact network network
	g <- erdos_renyi_network(1e5, 5)
	# Create transmission and reset time distributions
	psi <- lognormal_time(3, 30)
	rho <- lognormal_time(50, 20)
	# Create simulation and specifiy initial set of infections
	sim <- nextreaction_simulation(g, psi, rho)
	simulation_addinfections(sim, nodes=c(1), times=c(0.0))
	# Run simulation until time t=100 or 200,000 infections have occured
	r <- simulation_run(sim, list(time=100, total_infected=2e5))
	# Plot the number of infected nodes against time
	plot(r$time, r$infected, type='l')

*NEXTNetR* offers a range of common types of artifical networks such as Erdős–Rényi, Barabási–Albert and Watts–Strogatz, and can run simulations on arbitrary empirical weighted networks defined by an adjacency list. *NEXTNetR* also allows simulations on *temporal* networks, i.e. networks which change over time, possibly in response to epidemic events. Amongst the temporal networks currently supported by *NEXTNetR* are empirical networks defined by contact times between nodes, activity-driven networks, and networks defined by the proximity of diffusing particles.

# Reference

See <https://oist.github.io/NEXTNetR/>.