test_that("type checking", {
    g <- episimR_fully_connected_graph(1)

    expect_error(episimR_nextreaction_simulation(g, g, g));
})

test_that("next_reaction", {
    # Create network
    N <- 1
    g <- fully_connected_graph(N)

    # Create transmission time
    psi <- gamma_time(1, 2, 0.0);
    
    # Create simulator
    s <- nextreaction_simulation(g, psi, NULL);
    
    # Run
    simulation_addinfections(s, 1L, 0.0);
    r <- simulation_step(s, N);

    # Rudimentary check result
    expect_equal(nrow(r), N);
    expect_equal(ncol(r), 3);

    # Further runs should exit immediately since all nodes are infected
    r2 <- simulation_step(s, 1);
    expect_equal(nrow(r2), 0);
    expect_equal(ncol(r2), 3);        
})
