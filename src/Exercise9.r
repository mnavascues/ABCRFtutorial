##########################################################################
## Exercise 9:

target1_SS <- cbind(target1_pi, 
                    target1_TajimasD)
target2_SS <- cbind(target2_pi, 
                    target2_TajimasD)
sim_SS     <- cbind(sim1_pi,
                    sim1_TajimasD)
colnames(target1_SS) <- colnames(target2_SS) <-
  colnames(sim_SS) <- c("pi","TD")


abc_result1 <- abc(target=target1_SS,
                   param=sim_theta,
                   sumstat=as.matrix(sim_SS),
                   tol=tolerance,
                   transf = "log",
                   method="loclinear")
abc_result2 <- abc(target=target2_SS,
                   param=sim_theta,
                   sumstat=as.matrix(sim_SS),
                   tol=tolerance,
                   transf = "log",
                   method="loclinear")

sim_theta_adjusted_1 <- as.vector(abc_result1$adj.values)
sim_weights_1        <- abc_result1$weights
sim_theta_adjusted_2 <- as.vector(abc_result2$adj.values)
sim_weights_2        <- abc_result2$weights


# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted_1), breaks=seq(-2,2,0.05),
         col=Vermillion_transparency, freq=F, add=T, weight=sim_weights_1)
wtd.hist(log10(sim_theta_adjusted_2), breaks=seq(-2,2,0.05),
         col=Blue_transparency, freq=F, add=T, weight=sim_weights_2)
box()

