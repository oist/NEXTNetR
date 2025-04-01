# Create contact network network
g <- erdos_reyni_network(1e5, 5)
# Create transmission and reset time distributions
psi <- lognormal_time(3, 30)
rho <- lognormal_time(50, 20)
# Create simulation and specifiy initial set of infections
sim <- nextreaction_simulation(g, psi, rho)
simulation_addinfections(sim, nodes=c(1), times=c(0.0))
# Run simulation until 1000 individuals are infected, at for at most 1e5 steps
r <- simulation_run(sim, list(epidemic_steps=1e5, infected=1000))
