# Create contact network network
g <- erdos_renyi_network(1e5, 5)
# Create transmission and reset time distributions
psi <- lognormal_time(3, 30)
rho <- weibull_time(shape=5, scale=50)
# Create simulation and specifiy initial set of infections
sim <- nextreaction_simulation(g, psi, rho)
simulation_addinfections(sim, nodes=c(1), times=c(0.0))
# Run simulation until time t=100 or 200,000 infections have occured
r <- simulation_run(sim, stop=list(time=300, total_infected=300e3))
# Plot the number of infected nodes against time
plot(r$time, r$infected, type='l')
