
# Libraries ---------------------------------------------------------------

library(tidyverse)
library(rethinking)

# Binomial likelihoods ----------------------------------------------------

# Given that the probability is 50%, the probability of getting 6 positives out of 9 samples is 16.4%
dbinom(6, size = 9, prob = 0.5)


# Grid Approximation ------------------------------------------------------

# define the grid
p_grid <- seq(from = 0, to = 1, length.out = 20)
# define different priors to try
prior <- rep(1, 20)
prior <- ifelse(p_grid < 0.5, 0, 1)
prior <- exp(-5 * abs(p_grid - 0.5))
# compute likelihood for each p
likelihood <- dbinom(6, size = 9, prob = p_grid)
# compute the unstandardized posterior
unstd.posterior <- prior * likelihood
# standardized the posterior
posterior <- unstd.posterior / sum(unstd.posterior)
# plot the posterior
plot(p_grid, posterior, type = "b", 
    xlab = "probability of water",
    ylab = "posterior probability")
mtext("20 points")

#The Grid approximation scales poorly when there are more than a parameters


# Quadratic Approximation --------------------------------------------------

# The region near the peak of the posterior distribution will be nearly Gaussian, so we can approximate it with a Gaussian distribution.

# need to use the quap tool
library(rethinking)

# quap finds mode of posterior distribution for arbitrary fixed effect models and then produce an approximation of the full posterior using the quadratic curvature at the mode.

globe.qa <- quap(
    alist(
        W ~ dbinom( W+L ,p) ,  # binomial likelihood
        p ~ dunif(0,1)     # uniform prior
), data=list(W=6,L=3) )

# display summary of quadratic approximation
precis( globe.qa )


# MCMC ----------------------------------------

n_samples <- 1000
p <- rep( NA , n_samples )
p[1] <- 0.5
W <- 6
L <- 3
for ( i in 2:n_samples ) {
    p_new <- rnorm( 1 , p[i-1] , 0.1 )
    if ( p_new < 0 ) p_new <- abs( p_new )
    if ( p_new > 1 ) p_new <- 2 - p_new
    q0 <- dbinom( W , W+L , p[i-1] )
    q1 <- dbinom( W , W+L , p_new )
    p[i] <- ifelse( runif(1) < q1/q0 , p_new , p[i-1] )
}

# plot approximation using MCMC
dens(p, xlim = c(0, 1))
# plot real density
curve(dbeta(x, W+1, L+1), lty = 2, add = TRUE)


p_grid <- seq( from=0 , to=1 , length.out=1000 )
prob_p <- rep( 1 , 1000 )
prob_data <- dbinom( 6 , size=9 , prob=p_grid )
posterior <- prob_data * prob_p
posterior <- posterior / sum(posterior)

samples <- sample( p_grid , prob=posterior , size=1e5 , replace=TRUE )
w <- rbinom(1e4, size = 9, prob = samples)

simplehist(w)


# 95% crediable interval
quantile(samples, c(0.025, 0.975))
# 50% of the sample from the central
PI(samples, prob = 0.5)
# narrowest interval containing 50% of the samples
HPDI(samples, prob = 0.5)
