##########################################################################
## Exercise 8:

# plot theta and tajima's D from simulations
plot(log10(sim_theta),sim_SS,xlab=expression(log(theta*"'")),ylab="SS'",ylim=c(-4,4))
abline(h=target_SS,col=7)

# repeat ABC analysis with pi as SS
tolerance <- 0.1
target_SS <- target1_pi
sim_SS    <- as.matrix(sim1_pi); colnames(sim_SS) <- "pi"
abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = as.matrix(sim_SS),
                  tol     = tolerance,
                  transf  = "log",
                  method  = "loclinear")
sim_theta_adjusted_pi <- as.vector(abc_result$adj.values)
sim_weights_pi        <- abc_result$weights

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted_pi), breaks=seq(-2,2,0.05), col=Blue_transparency, freq=F, add=T, weight=sim_weights_pi)
box()
wtd.hist(x      = log10(sim_theta_adjusted_D),
         breaks = seq(-2,3,0.05),
         col    = Vermillion_transparency,
         freq   = FALSE, add = TRUE,
         weight = sim_weights_D); box()



