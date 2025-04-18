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

test_that("userdefined", {
    al <- list()
    al[[1]] <- c(2,3,4)
    al[[2]] <- c(1,3)
    al[[3]] <- c(1,2,4)
    al[[4]] <- c(1,3)

    g <- userdefined_graph(al)
    expect_equal(graph_outdegree(g, 1:4), c(3,2,3,2))
    expect_equal(graph_neighbour(g, rep(1:4, each=4), rep(1:4, times=4)),
                 c(2,3,4,NA, 1,3,NA,NA, 1,2,4,NA, 1,3,NA,NA))
    expect_equal(graph_adjacencylist(g)$neighbours, al)
})
