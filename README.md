# Introduction

*NEXTNetR* (*N*ext-reaction-based *E*pidemics e*X*tended to *T*emporal *Net*works) is an R package for the efficient simulation of epidemics on complex networks with arbitrary transmission and recovery time distributions. *NEXTNetR* is an R wrapper around the C++ library NEXTNet <https://github.com/oist/NEXTNet>.

# Installation

The latest version of *NEXTNetR* can be installed directly from Github by executing the following in R:

    install.packages(c("devtools", "git2r"))
    devtools::install_git("git@github.com:oist/NEXTNetR")

# Synopsis

The following R script sets up and runs a simulation on an Erdös-Reyni network with lognormally distributed transmission time (i.e. the generation time between subsequent infections) and lognormally distribured reset time (the time until an infected node becomes susecptible again).

    # Load library
    library(NEXTNetR)

    # Parameters, which are network size (N), average number of contacts (K), 
    # mean infection transmission time (Mi) and variance (Vi),
    # mean reset time (Mr) and variance (Vr)
    N <- 10000
    K <- 5
    Mi <- 3
    Vi <- 30
    Mr <- 50
    Vr <- 20

    # Create contact network network
    g <- erdos_reyni_network(N, K)

    # Create transmission and reset time distributions
    psi <- lognormal_time(Mi, Vi)
    rho <- lognormal_time(Mr, Vr)

    # Create simulation and specifiy initial set of infections
    sim <- nextreaction_simulation(g, psi, rho)
    simulation_addinfections(sim, nodes=c(1), times=c(0.0))

    # Run simulation for 1e5 steps
    r <- simulation_run(sim, list(epidemic_steps=1e5))

    # Plot the number of infected nodes against time
    plot(r$time, r$infected, type='l')

# Reference

To run a simulation, two objects must be created first, a *network*, and a *transmission time* distribution.

## Creating and querying contact network networks

To create networks, the following functions are available

    erdos_reyni_network(size, avg_degree)
    fully_connected_network(size)
    acyclic_network(size, avg_degree, reduced_root_degree)
    configmodel_network(degrees)
    barabasialbert_network(size, m)
    configmodel_clustered_network(degrees, alpha_or_ck_or_triangles, beta)
    cubiclattice2d_network
    ...
    cubiclattice8d_network
    brownian_proximity_temporalnetwork(size, avg_degree, contact_radius, D, dt)
    empirical_temporalnetwork(file, finite_duration, dt)
    empirical_network(filename)
    adjacencylist_network(adjacencylist)
    
For user-defined networks, the topology is defined by the adjacency list, which must be a *list* whose length defines the number of nodes. The *i*-th element of the list must be an *integer vector* with elements from [*1*,*n*] listing the neighbours of the *i*-th node. For example, `list(c(2,3), c(1,3), c(1,2))` represents the complete network with 3 nodes.

The topology of an existing contact network network *g* returned by one of the functions above can be inspected with

    network_outdegree(network, nodes)
    network_neighbour(network, nodes, indices)
    network_adjacencylist(network)

For network embedded into d-dimensional space, the coordinates of the nodes can be queried with

    network_coordinates(network, nodes)

## Creating and querying time distributions

To create a t time distribution, the following functions are available

    exponential_time(lambda)
    lognormal_time(mean, var, p_infinity)
    gamma_time(mean, var, p_infinity)
    generic_time(density, survivalprobability, probability_is_trinary,
                 survivalquantile, quantile_is_trinary,
                 sample, p_infinity)

From a time distribution object *time* created by one of these function *n* independent random samples can be taken with `time_sample`. When sampling, the distribution can be conditioned on values *\>= t* for any *t \>= 0*, and can be modified to return the minimum of *m* independent samples instead.

    time_sample(n, time, t=0, m=1)

The density, hazardrate, cumulative survival function (i.e. 1 - CDF) and survival quantile function can be evaluated using the functions

    time_density(time, tau)
    time_hazardrate(time, tau)
    time_survivalprobability(time)
    time_survivalquantile(time, p, t = 0, m = 1)

## Creating, querying and running simulations

### Creating simulations using the NextReaction algorithm

Given a contact network *network*, transmission time distribution *psi* and reset time distribution *rho*, a simulation using the NextReaction algorithm is created with

    nextreaction_simulation(network, psi, rho = NULL, options = list())

The *options* parameter must be a named list specified the option names and their values. The supported options are:

*shuffle_neighbours*
: Shuffle neighbour order  upon infection of the node. Default *true*.

*edges_concurrent*
: Whether to activate all outgoing edges simultaenously or sequentially. If set to true, neighbours are implicitly shuffled and *shuffle_neighbours* thus has no effect. Default *false*.

*SIR*
: Whether to simulate using a SIR (*S*usceptiple, *I*nfected, *R*ecovered) or SIS (*S*usceptiple, *I*nfected, *S*usceptible) model. If set to true, the model is SIR, i.e. recovered individuals do not become susceptible again.

*static_network*
: Whether to treat a temporal network as static. By default, if a temporal network is specified, the network will evolve over time as the simulation commences. If *static_network* is set to true, the network is not evolved.

### Creating simulations using the nMGA algorithm

Given a contact network *network*, transmission time distribution *psi* and reset time distribution *rho*, a simulation using the nMGA algorithm is created with

    nmga_simulation(network, psi, rho = NULL, options = list())

The *options* parameter must be a named list which specifies the option names and their values. The supported options are:

*approx_threshold*
: Threshold for infected nodes at which the approximate algorithm is used

*max_dt*
: Maximum timestep allowed for the approximate algorithm

*tauprec*
: Numerical precision used to invert the CDF in the exact algorithm

### Creating simulations using the REGIR algorithm

The R wrapper currently does not support the REGIR algorithm

### Querying and running simulations, obtaining results

Before running a simulation, an initial set of infected nodes plus the times at which these become infected must be specified.

    simulation_addinfections(sim, nodes, times)

The simulation can then be run until one of a list of possible stopping conditions is reached with

    simulation_run(sim, stopping_condition)
    
The *stopping_condition* parameter must be a named list which specifies the stopping conditions to apply. The supported conditions are:

*epidemic_steps*
: Stop after the specified number of epidemic steps (i.e. infections or recoveries/resets)

*network_steps*
: Stop after the specified number of network topology changes

*time*
: Stop at the specified simulation time

*infected*
: Stop when the current number of infected nodes reaches the specified threshold. Recovered/reset nodes to not count towards infected nodes here.

*total_infected*
: Stop when the cumulative number of infection events reaches the specified threshold.

*total_reset*
: Stop when the cumulative number of recover/reset events reaches the specified threshold.

Running a simulation returns a `data.frame` which lists the event that occurred in each step in ascending order of time of occurrence. The result contains the 6 columns *time*, *kind* (infection or reset), *node* (1-based index of affected node), *total_infected* (total number of infection so far), *total_reset* (total number of reset so far), *infected* (number of currently infected nodes). *simulation_step* returns after the specified number of steps, but can be continued by simply calling *simulation_step* again.

The constituent parts of a simulation (contact network network and time distribution) and the current state of the simulation
(number of infected nodes and whether a node is infected or not) can be queried with

    simulation_transmissiontime(sim)
    simulation_resettime(sim)
    simulation_network(sim)
    simulation_options(sim)
    simulation_ninfected(sim)
    simulation_isinfected(sim, nodes)
