##########################################################################
## Exercise 4:

# Prior believe that the coin is fair
sh1 <- sh2 <- 5
sim_parameter_p2 <- rbeta(100000,sh1,sh2)
sim_data_p2      <- rbinom(n    = length(sim_parameter_p2),
                           size = 10,
                           prob = sim_parameter_p2)

sim_parameter_kept_p2 <- sim_parameter_p2[which(sim_data_p2==
                                                  heads_count)]
sim_data_kept_p2      <- sim_data_p2[which(sim_data_p2==
                                             heads_count)]

hist(x      = sim_parameter_p2,
     breaks = seq(0,1,0.02),
     col    = "grey",
     freq   = F,
     ylim   = c(0,5),
     main   = "",
     xlab   = expression(italic(p)),
     ylab   = "probability density")
hist(sim_parameter_kept_p2, breaks=seq(0,1,0.02), col=Blue_transparency, freq=F, add=T)
box()
lines(x   = seq(0,1,0.001),
      y   = dbeta(x=seq(0,1,0.001),sh1+heads_count,sh2+tails_count),
      col = 6,
      lwd = 3)
lines(x   = seq(0,1,0.001),
      y   = dbeta(x=seq(0,1,0.001),1+heads_count,1+tails_count),
      col = 7,
      lwd = 3)
p_hat  <- median(sim_parameter_kept_p2)
p_95CI <- quantile(sim_parameter_kept_p2,probs=c(0.025,0.975))
cat(paste(c("Point estimate of p:",round(p_hat,2),"; 95%CI:",
            round(as.vector(p_95CI),2),"\n"),sep=" ",collapse = " "))

