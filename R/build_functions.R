## Libraries

require(MASS)
require(truncnorm)
require(igraph)

## ===============================================================
## Functions to build different types of matrices and growth rates,
## The matrices are required to be negative definite A + A^T < 0
## ===============================================================


## ===================================================
## Growth rates
## ========================================================
build.normal.gr <- function(n, params){
    u <- params$u
    sd <- params$sd
    return (rnorm(n, mean = u, sd = sd))
}

build.truncated.normal.gr <- function(n, params){
    u <- params$u
    sd <- params$sd
    return (rtruncnorm(n, a = 0, mean = u, sd = sd))
}


## ===========================================================================
## Matrices with different structure and parameterization upon that structure
## ===========================================================================

### general functions
get.build.structure <- function(n){
    switch(n,
        "1" = build.no.structure,
        "2" = build.random.structure,
        "3" = build.powerlaw,
        "4" = build.2.blocks,
        "5" = build.full.bipartite,
        "6" = build.full.inblocks)
}

get.build.matrix <- function(n){
    switch(n,
           "1" = build.fixed.matrix,
           "2" = build.normal,
           "3" = build.cor.normal)
}


build.random.structure.dstable <- function(n, params){
    build.structure <- params$build.structure
    S <- build.structure(n, params)
    params$S <- S
    return (build.structure.dstable(n, params))
}

build.structure.dstable <- function(n, params){
    build.matrix <- params$build.matrix
    A <- build.matrix(n, params)
    S <- params$S
    A <- S * A
    D <- A + t(A)
    e <- max(eigen(D, symmetric = TRUE, only.values = TRUE)$values)
    perturb <- 1e-6
    if (e >= 0){
        A <- A - diag(n) * (e+perturb) / 2
    }
    return (A)
}
### Types of matrices

build.fixed.matrix <- function(n, params){
    u <- params$u
    d <- params$d
    A <- matrix(u, nrow = n, ncol = n)
    diag(A) <- d
    return (A)
}

build.normal <- function(n, params){
    u <- params$u
    sd <- params$sd
    d <- params$d
    A <- matrix(rnorm(n^2, mean = u, sd = sd), nrow = n ) + d * diag(n)
    return (A)
}

build.cor.normal <- function(n, params){
    u <- params$u
    sd <- params$sd
    d <- params$d
    p <- params$p
    mu <- c(u, u)
    Sigma <- matrix(c(sd^2, p * sd^2, p * sd^2, sd^2), nrow = 2)
    A <- diag(rnorm(n, mean = u, sd = sd) + d)
    for(i in seq(1, n - 1)) {
        for (j in seq(i + 1, n)){
            Z <- mvrnorm(1, mu, Sigma)
            A[i,j] <- Z[1]
            A[j,i] <- Z[2]
        }
    }
    return (A)
}

### Types of Structures

#### ER random graph
build.random.structure <- function(n, params){
    connectance <- params$c
    S <- matrix(0, nrow = n, ncol = n)
    S[upper.tri(S)] <- rbinom(n * (n - 1)/2,size = 1, p = connectance)
    S <- S + t(S)
    diag(S) <- 1
    return(S)    
}


#### Mean field matrix
build.fixed.matrix.structure <- function(n, params){
    A <- build.fixed.matrix(n, params)
    build.structure <- params$build.structure
        
    while(TRUE){
        S <- build.structure(n, params)
        B <- A * S
        e <- max(eigen(B, symmetric = TRUE, only.values = TRUE)$values)
        if (e < -1e-15)
            break
    }
    return (B)
}

#### Power law degree matrix
build.powerlaw <- function(n, params){
    gamma <- params$exponent
    c <- params$c
    edges <- c * n * (n-1)/2
    PG <- sample_fitness_pl(n, edges, gamma)
    M <- matrix(get.adjacency(PG), nrow = n)
    diag(M) <- 1
    return (M)
}


### Completely filled matrix, i.e ER with C = 1
build.no.structure <- function(n, params){
    return (1)
}


### Matrix with two blocks
build.2.blocks <- function(n, params){
    M <- matrix(0, n, n)
    c1 <- params$c_within
    c2 <- params$c_between
    ratio <- params$block_ratio
    m <- round(ratio * n)
    ## build connection within blocks
    params$c <- c1
    A <- build.random.structure(m, params)
    B <- build.random.structure(n-m, params)
    ## build connection among blocks
    C <- matrix(rbinom(m * (n -m), size = 1 , p = c2), m , n - m)
    M[1:m , 1:m] <- A
    M[1:m, (m+1):n] <- C
    M[(m+1):n, 1:m] <- t(C)
    M[(m+1):n, (m+1):n] <- B
    
    return (M)
}



### For a given connectance, get the maximally bipartite matrix
### assuming within block and between block ER graphs
build.full.bipartite <- function(n, params){
    c <- params$c
    m <- round(params$block_ratio * n)
    v <- m * (m - 1) + (n - m) * (n - m - 1)
    b <- 2 * (n - m) * m
    c_ratio <- params$c_ratio
    w <- c_ratio * v + b
    c_within <- (c * n * (n - 1))/w
    params$c_within <- c_within
    params$c_between <- c_within * c_ratio
    return (build.2.blocks(n, params))
}

### For a given connectance, get the maximally modular matrix
### assuming within block and between block ER graphs

build.full.inblocks <- function(n, params){
    c <- params$c
    m <- round(params$block_ratio * n)
    v <- m * (m - 1) + (n - m) * (n - m - 1)
    b <- 2 * (n - m) * m
    c_ratio <- params$c_ratio
    w <- c_ratio * b + v
    c_between <- c * n * (n - 1) / w
    params$c_between <- c_between
    params$c_within <- c_between * c_ratio
    return (build.2.blocks(n, params))    
}




## =================================================
## Translate from LV to Replicator dynamics
## =================================================

build.rep <- function(A, r, n){
    M <- matrix(0, n+1, n+1)
    M[1:n, n+1] <- r
    M[1:n, 1:n] <- A
    return (M)
}

