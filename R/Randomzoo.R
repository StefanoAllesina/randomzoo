## Load required packages and files
require(tidyverse)
source("lemke-howson.R")
source("build_functions.R")

## ============================================================
## Example of dynamical pruning under lotka volterra dynamics
## ============================================================

## Set parameter values
n <- 1000 ## number of species in the original community
params_interactions <- list(build.structure = build.no.structure, ## completely filled matrix
               build.matrix = build.normal,          ## normal matrix
               u = 0,                                ## mean of interactions
               d = -1,                               ## original mean diagonal (this may change to ensure stability)
               sd = 1                                ## standard deviation of interactions
               )
params_growth_rates <- list(u = 0,                   ## mean of growth rates
                            sd = 1                   ## sd of growth rates   
                            )
A <- build.random.structure.dstable(n, params_interactions) ## set interaction matrix
r <- build.normal.gr(n, params_growth_rates)                ## set growth rate
M <- build.rep(A, r, n)                                     ## set payoff matrix for replicator dynamics

### Run dynamics
y <- lemkeHowson_symmetric(M)  ## get nash equilibrium of the game
S <- y$subset                  ## Subset of coexisting species
x <- toLV(y$eq, n)[S]          ## equilibrium of LV model
B <- M[S, S]                   ## Final interaction Matrix
v <- r[S]                      ## Final vector of growth rates

### Plotting

## Compare eigenvalue distribution

evals_pruned <- eigen(B, only.values = TRUE)$values  
evals_original  <- eigen(A, only.values = TRUE)$values

Evals_pruned <- data_frame(x = Re(evals_pruned),
                           y = Im(evals_pruned), type = "pruned")
Evals_original <- data_frame(x = Re(evals_original),
                             y = Im(evals_original), type = "original")
Evals <- Evals_pruned %>% bind_rows(Evals_original)

Evals %>% ggplot(aes(x = x, y = y)) + geom_point()  +
    facet_grid(.~type) + theme_bw() + coord_fixed()

## Compare distribution of coefficients

Dist_coeff_o <- data_frame(a = A[upper.tri(A)], type = "original")
Dist_coeff <- data_frame(a = B[upper.tri(B)], type = "pruned")
Dist_coeff <- Dist_coeff %>% bind_rows(Dist_coeff_o)

Dist_coeff %>% ggplot(aes(x = a, fill = type, colour = type)) +
    geom_density(alpha = 0.5) + facet_grid(.~type) + theme_bw()


## Compared distribution of growth rates

Dist_growth_rates <- data_frame(r = v, type = "pruned") %>% bind_rows(data_frame(r = r, type = "original"))

Dist_growth_rates %>% ggplot(aes(x = r, fill = type, colour = type)) +
    geom_density(alpha = 0.5) + facet_grid(.~type) + theme_bw()



## ==============================================================
## Numerical computation of distribution of attractor's sizes
## ==============================================================

get_distribution_sizes_attractors <- function(n, nsim,
                                              params_interactions,
                                              params_growth_rates){
    sizes <- map_dbl(1:nsim, function(i){
        A <- build.random.structure.dstable(n, params_interactions)
        r <- build.normal.gr(n, params_growth_rates)
        M <- build.rep(A, r, n)
        y <- lemkeHowson_symmetric(M)
        return (length(y$subset))
    })

    D <- data_frame(n = n, k = sizes, mu = params_interactions$u, sd = params_interactions$sd, gamma = params_growth_rates$u)
    
    ## convert to frequencies, if needed the grouping can be done for all different parameters
    D <- D %>% group_by(k) %>% summarise(count = n()) %>% mutate(freq = count / sum(count)) %>%
        ungroup()

    return (D)
}

n <- 15
nsim <- 50000
set.seed(171) ## set random seed for simulations
attractors <- get_distribution_sizes_attractors(n, nsim,
                                                params_interactions,
                                                params_growth_rates )

## In this case we can get a binomal expectation
attractors <- attractors %>% mutate(expectation = choose(n, k) *  0.5^n)

## Plotting
attractors %>% ggplot(aes(x = k, y = freq)) +
    geom_bar(stat = "identity", alpha = 0.5, colour = "red", fill = "red") +
    geom_point(aes(y = expectation), shape = 4, size = 2) + theme_bw()












               
              
               




