test_that("fully_connected", {
    # Create network
    N <- 5
    g <- fully_connected_network(N)
    
    # Check that all nodes have outdegree N-1
    d <- network_outdegree(g, 1:N)
    expect_equal(d, rep(N-1, N))
    
    # Check that all nodes have all other nodes except themselves as neighbours
    n <- network_neighbour(g, rep(1:N, each=N-1), rep(1:(N-1), times=N))
    expect_equal(length(n), N*(N-1))
    m <- matrix(n, nrow=N, ncol=N-1, byrow=TRUE)
    for(i in 1:N)
        expect_equal(sort(m[i,]), setdiff(1:N, i))
    
    # Check adjacencylist generation
    a <- network_adjacencylist(g)
    for(i in 1:N)
        expect_equal(sort(a$neighbours[[i]]), if (i < N) ((i+1):N) else integer())
})

test_that("adjacencylist", {
    al <- list()
    al[[1]] <- c(2,3,4)
    al[[2]] <- c(1,3)
    al[[3]] <- c(1,2,4)
    al[[4]] <- c(1,3)

    g <- adjacencylist_network(al, is_undirected=FALSE)
    expect_equal(network_outdegree(g, 1:4), c(3,2,3,2))
    expect_equal(network_neighbour(g, rep(1:4, each=4), rep(1:4, times=4)),
                 c(2,3,4,NA, 1,3,NA,NA, 1,2,4,NA, 1,3,NA,NA))
    expect_equal(network_adjacencylist(g)$neighbours, al)
})
