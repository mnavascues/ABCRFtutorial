##########################################################################
## Exercise 10:

posterior_theta <- sample(sim_theta_adjusted_2,size=1000,replace=T,prob=sim_weights_2)

if (file.exists("data/abc_sims_posterior.txt")){
  file_removed <- file.remove("data/abc_sims_posterior.txt")
}
ms(nsam       = sample_size, 
   opt        = "-t tbs",
   tbs.matrix = cbind(posterior_theta),
   temp.file  = "data/abc_sims_posterior.txt")

msout <- ms.inp.multi(sample_size, 1000, ms.output.file="data/abc_sims_posterior.txt")
# simulated summary statistics from posterior
posterior_S        <- S(msout)
posterior_pi       <- thetaPi(msout)
posterior_NH       <- NH(msout)
posterior_TajimasD <- tajimaD(msout, thetaW(msout), posterior_pi)
posterior_FayWuH   <- fayWuH(msout)
posterior_FuLiD    <- fuliD(msout, thetaS1(msout) )

hist(posterior_NH, breaks=seq(0,40,1), col="grey", freq=F, xlab="Number of haplotypes", main="")
abline(v=target2_NH,col=7,lwd=3)
box()
hist(posterior_FuLiD, breaks=seq(-5,5,0.2), col="grey", freq=F, xlab="Fu and Li's D", main="")
abline(v=target2_FuLiD,col=7,lwd=3)
box()
