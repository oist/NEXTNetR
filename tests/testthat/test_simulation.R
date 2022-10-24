test_that("next_reaction", {
    # Create network
    N <- 10
    g <- episimR_fully_connected_graph(N)
    
    # Create transmission time
    psi <- episimR_gamma_time(1, 2, 0.0);
    
    # Create simulator
    s <- episimR_nextreaction_simulation(g, psi, NULL);
    
    # Run
    episimR_simulation_addinfections(s, 1, 0);
    r <- episimR_simulation_step(s, N-1);
    
    print(r)
})
