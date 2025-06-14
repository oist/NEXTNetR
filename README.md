# Introduction

*NEXTNetR* (**N**ext-reaction-based **E**pidemics e**X**tended to **T**emporal **Net**works) is an R package for the efficient simulation of epidemics on complex networks (including weighted and temporal networks) with arbitrary transmission and recovery time distributions. *NEXTNetR* is an R wrapper around the C++ library [*NEXTNet*](https://github.com/oist/NEXTNet).

See the [*NEXTNetR* website](https://oist.github.io/NEXTNetR/) for a reference and usage examples.

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

See [Getting Started](https://oist.github.io/NEXTNetR/articles/NEXTNetR.html) for a step-by-step walkthrough of NEXTNetR's features.

*NEXTNetR* offers a range of common types of artifical networks such as Erdős–Rényi, Barabási–Albert and Watts–Strogatz, and can run simulations on arbitrary empirical weighted networks defined by an adjacency list. *NEXTNetR* also allows simulations on *temporal* networks, i.e. networks which change over time, possibly in response to epidemic events. Amongst the temporal networks currently supported by *NEXTNetR* are empirical networks defined by contact times between nodes, activity-driven networks, and networks defined by the proximity of diffusing particles.

# Documentation

See <https://oist.github.io/NEXTNetR/>.
