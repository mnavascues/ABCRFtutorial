##########################################################################
## Exercise 7:

# Let's invcrease the tolerance
tolerance <- 0.5

abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = as.matrix(sim_SS),
                  tol     = tolerance,
                  method  = "rejection")
sim_theta_kept_tol2 <- as.vector(abc_result$unadj.values)
sim_SS_kept_tol2    <- abc_result$ss

plot(log10(sim_theta),sim_SS,xlab=expression(log[10](theta*"'")),ylab="SS'",ylim=c(0,50))
abline(h=target_SS,col=6)
{abline(h=max(abc_result$ss),col=6,lty=2)
  abline(h=min(abc_result$ss),col=6,lty=2)}
points(log10(sim_theta_kept_tol2),sim_SS_kept_tol2,col=6)



# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log(theta)),
     ylab="probability density")
hist(log10(sim_theta_kept_tol2), breaks=seq(-1,2,0.05), 
     col=Blue_transparency, freq=F, add=T)
hist(log10(sim_theta_kept_tol1), breaks=seq(-1,2,0.05), 
     col=Vermillion_transparency, freq=F, add=T)
box()
