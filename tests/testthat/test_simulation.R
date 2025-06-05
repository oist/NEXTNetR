test_that("type checking", {
    g <- fully_connected_network(1)

    expect_error(simulation(g, g, g));
})

test_that("options", {
  g <- fully_connected_network(1)
  psi <- gamma_time(1, 2, 0.0);
  
  # Implicitly set default value
  s <- simulation(g, psi, NULL)
  expect_equal(simulation_options(s)$shuffle_neighbours, TRUE)
  
  # Explicitly set default value
  s <- simulation(g, psi, NULL, list(shuffle_neighbours=TRUE))
  expect_equal(simulation_options(s)$shuffle_neighbours, TRUE)
  
  # Explicitly set non-default value
  s <- simulation(g, psi, NULL, list(shuffle_neighbours=FALSE))
  expect_equal(simulation_options(s)$shuffle_neighbours, FALSE)
  
  # Warn on unexpected options
  expect_warning(simulation(g, psi, NULL, list(no_such_option=NULL)))
})

test_that("basic simulation", {
    # Create network
    N <- 10
    g <- fully_connected_network(N)

    # Create transmission time
    psi <- gamma_time(1, 2, 0.0);
    
    # Create simulator
    s <- simulation(g, psi, NULL);
    
    # No node should be infected
    expect_equal(simulation_isinfected(s, 1:N), rep(FALSE, N))
    
    # Setup initial infection
    simulation_addinfections(s, 1L, 0.0);
    
    # Run first step
    expect_equal(simulation_isinfected(s, 1:N), rep(FALSE, N))
    r1 <- simulation_run(s, stop=list(epidemic_steps=1));
    expect_equal(simulation_ninfected(s)$infected, 1)
    expect_equal(simulation_ninfected(s)$total_infected, 1)
    expect_equal(simulation_ninfected(s)$total_reset, 0)
    expect_equal(simulation_isinfected(s, 1:N), c(TRUE, rep(FALSE, N-1)))
    expect_equal(nrow(r1), 1);
    expect_equal(ncol(r1), 10);
    expect_equal(r1$infected, 1);
    expect_equal(r1$total_infected, 1);
    expect_equal(r1$total_reset, 0);
    
    # Run further steps
    r2 <-simulation_run(s, stop=list(epidemic_steps=N));
    expect_equal(simulation_ninfected(s)$infected, N)
    expect_equal(simulation_ninfected(s)$total_infected, N)
    expect_equal(simulation_ninfected(s)$total_reset, 0)
    expect_equal(simulation_isinfected(s, 1:N), rep(TRUE, N))
    r <- rbind(r1, r2)

    # Rudimentary check result
    expect_equal(nrow(r), N);
    expect_equal(ncol(r), 10);
    expect_equal(r$infected, 1:N);
    expect_equal(r$total_infected, 1:N);
    expect_equal(r$total_reset, rep(0, N));

    # Further runs should exit immediately since all nodes are infected
    r3 <- simulation_run(s, stop=list(epidemic_steps=1));
    expect_equal(simulation_ninfected(s)$infected, N)
    expect_equal(simulation_ninfected(s)$total_infected, N)
    expect_equal(simulation_ninfected(s)$total_reset, 0)
    expect_equal(simulation_isinfected(s, 1:N), rep(TRUE, N))
    expect_equal(nrow(r3), 0);
    expect_equal(ncol(r3), 10);
    expect_equal(r3$infected, r3$total_infected - r3$total_reset)
})
