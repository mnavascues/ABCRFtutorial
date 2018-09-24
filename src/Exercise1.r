##########################################################################
## Exercise 1:

# data from 100 tosses
toss <- sample(c('H','T'),size=100,replace=T)
total_tosses <- length(toss)
heads_count  <- length(which(toss=="H"))
likelihood_profile <- flip.coin.likelihood(n = total_tosses,
                                           k = heads_count,
                                           p = seq(0,1,0.001))
plot(x    = likelihood_profile$p,
     y    = likelihood_profile$likelihood, type="l")

# data from 1000 tosses
toss <- sample(c('H','T'),size=1000,replace=T)
total_tosses <- length(toss)
heads_count  <- length(which(toss=="H"))
likelihood_profile <- flip.coin.likelihood(n = total_tosses,
                                           k = heads_count,
                                           p = seq(0,1,0.001))
lines(x    = likelihood_profile$p,
      y    = likelihood_profile$likelihood)
