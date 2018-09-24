##########################################################################
## Exercise 6:

# ABC rejection using pi

target_SS <- target1_pi
sim_SS    <- sim1_pi

plot(log10(sim_theta),
     sim_SS,
     xlab=expression(log[10](theta*"'")),ylab="SS'",ylim=c(0,50))

sim_theta_kept <- sim_theta[which(sim_SS==target_SS)]
sim_SS_kept    <- sim_SS[which(sim_SS==target_SS)]
abline(h=target_SS,col=7)
points(log10(sim_theta_kept),sim_SS_kept,col=7)

cat(paste("The number of simulations kept is",
          length(sim_theta_kept), "out of",  length(sim_SS), "\n"))

if (length(sim_SS_kept)>0){
  hist(log10(sim_theta), breaks=seq(-1,2,0.05), col="grey",
       freq=F, ylim=c(0,6), main="", xlab=expression(log[10](theta)),
       ylab="probability density")
  hist(log10(sim_theta_kept), breaks=seq(-1,2,0.05), col=Vermillion_transparency, freq=F, add=T)
  box()
}

# ABC rejection using S, pi and NH

sim_kept <- which(sim1_S==target1_S)
cat(paste("For the moment,", length(sim_kept),
          "simulations kept out of",  length(sim_SS), "\n"))
sim_kept <- intersect(sim_kept, which(sim1_pi==target1_pi))
cat(paste("For the moment,", length(sim_kept),
          "simulations kept out of",  length(sim_SS), "\n"))
sim_kept <- intersect(sim_kept, which(sim1_NH==target1_NH))
cat(paste("For the moment,", length(sim_kept),
          "simulations kept out of",  length(sim_SS), "\n"))
# etc
