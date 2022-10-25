test_that("fully_connected", {
    # Create graph
    N <- 5
    g <- fully_connected_graph(N)
    
    # Check that all nodes have outdegree N-1
    d <- graph_outdegree(g, 1:N)
    expect_equal(d, rep(N-1, N))
    
    # Check that all nodes have all other nodes except themselves as neighbours
    n <- graph_neighbour(g, rep(1:N, each=N-1), rep(1:(N-1), times=N))
    expect_equal(length(n), N*(N-1))
    m <- matrix(n, nrow=N, ncol=N-1, byrow=TRUE)
    for(i in 1:N)
        expect_equal(sort(m[i,]), setdiff(1:N, i))
    
    # Check adjacencylist generation
    a <- graph_adjacencylist(g)
    for(i in 1:N)
        expect_equal(sort(a$neighbours[[i]]), setdiff(1:N, i))
})
