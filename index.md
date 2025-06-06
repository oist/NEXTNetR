# Introduction

`{NEXTNetR}` (**N**ext-reaction-based **E**pidemics e**X**tended to **T**emporal **Net**works) is an R package for the efficient simulation of epidemics on complex networks (including weighted and temporal networks) with arbitrary transmission and recovery time distributions. *NEXTNetR* is an R wrapper around the C++ library [NEXTNet](https://github.com/oist/NEXTNet).

# Installation

If [Git](https://git-scm.com/downloads) is available, the [latest released](https://github.com/oist/NEXTNetR/releases) version of *NEXTNetR* can be installed directly from Github by executing the following in R:

    install.packages("remotes")
    remotes::install_github("oist/NEXTNetR", ref="latest-release")

Alternatively, download the [latest released](https://github.com/oist/NEXTNetR/releases) version of *NEXTNetR-v\<version\>-pkg.tar.gz*. Then make sure all required dependencies are installed with `install.packages(c("BH", "cpp11", "rappdirs"))` and install *NEXTNetR* on the command line (*not* within R) with

    R CMD INSTALL NEXTNetR-v<version>-pkg.tar.gz
   
Since *NEXT-Net* is implemented in C++, a C++ compiler is required to install *NEXTNetR*. On Linux a compiler should typically be already available, on Mac OS R a suitable compiler is provided by [XCode](https://developer.apple.com/xcode/) or the [XCode Command Line Tools](https://mac.install.guide/commandlinetools/), and on Windows by [RTools](https://cran.rstudio.com/bin/windows/Rtools/).

# Synopsis

The following minimal example simulated an epidemic on an Erdős–Rényi network with lognormally distributed transmission time

	library(NEXTNetR)
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
