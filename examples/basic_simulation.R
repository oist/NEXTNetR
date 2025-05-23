# Create contact network network
nw <- erdos_renyi_network(1e5, 5)
# Create transmission and reset time distributions
psi <- lognormal_time(6, 30, 0.1)
rho <- weibull_time(shape=5, scale=50)
# Create simulation and specifiy initial set of infections
sim <- simulation(nw, psi, rho)
simulation_addinfections(sim, nodes=c(1), times=c(0.0))
# Run simulation until time t=100 or 200,000 infections have occured
r <- simulation_run(sim, stop=list(time=300, total_infected=300e3))
# Plot the number of infected nodes against time
plot(r$time, r$infected, type='l')
