##########################################################################
## Exercise 3:

# 100 simulations
likelihood_profile_approx <-
  flip.coin.likelihood.approx(n   = total_tosses,
                              k   = heads_count,
                              p   = seq(0,1,0.01),
                              rep = 10)
plot(x    = likelihood_profile$p,
     y    = likelihood_profile$likelihood,
     xlab = expression(italic(p)), ylab = "Likelihood",
     lwd  = 2, type = "l")
lines(x   = likelihood_profile_approx$p ,
      y   = likelihood_profile_approx$likelihood,
      col = 7, lwd = 2)

# 10000 simulations
likelihood_profile_approx <-
  flip.coin.likelihood.approx(n   = total_tosses,
                              k   = heads_count,
                              p   = seq(0,1,0.01),
                              rep = 10000)
lines(x   = likelihood_profile_approx$p ,
      y   = likelihood_profile_approx$likelihood,
      col = 4, lwd = 2)
