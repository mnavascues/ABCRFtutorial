################################################
# Workshop on Approximate Bayesian Computation 
#     from rejection to random forests
################################################

# Miguel Navascués
# miguel.navascues@inra.fr
# INRA

## ---- rcode_0_A
sessionInfo()$running
sessionInfo()$platform
R.version.string
.Platform$GUI
## ---- rcode_0_B
# install.packages(c("phyclust","abc","weights","abcrf","tree"))
library("phyclust") # for a version of Hudson ms (coalescent simulator)
library("abc")      # for approximate Bayesian computation
library("weights")  # for weighted histograms
library("abcrf")    # for ABC with random forests
library("tree")     # for CART
## ---- rcode_0end

sessionInfo()$otherPkgs$phyclust$Version
sessionInfo()$otherPkgs$abc$Version
sessionInfo()$otherPkgs$weights$Version
sessionInfo()$otherPkgs$abcrf$Version
sessionInfo()$otherPkgs$tree$Version

recordSessionInfo <- sessionInfo()
save(recordSessionInfo,file="results/sessionInfo.RData")


#################################
# Introduction to the likelihood 
#################################

## ---- rcode_1
toss         <- c("T","H","H","T","H","H","T","H","H","H")

total_tosses <- length(toss)
heads_count  <- length(which(toss=="H"))
tails_count  <- total_tosses - heads_count
## ---- rcode_1end

cat(paste("number of tosses:",total_tosses,"\n"))
cat(paste("number of heads:",heads_count,"\n"))
cat(paste("number of tails:",tails_count,"\n"))

## ---- rcode_2
p_heads <- 0.5
p_tails <- 1-p_heads
p_combo <- p_heads^heads_count * p_tails^tails_count
## ---- rcode_2end
message <- c("Probability of combination (",toss,"):",p_combo,"\n")
cat(paste(message, sep=" ",collapse=" "))


## ---- rcode_3
combinations <- choose(total_tosses, heads_count)
## ---- rcode_3end
cat(paste("There are",combinations,"combination of", heads_count, "head and", tails_count, "tails\n"))


## ---- rcode_4
likelihood <- combinations *
              p_heads^heads_count * p_tails^tails_count
## ---- rcode_4end
cat(paste("Likelihood of p=",p_heads,"given",heads_count,"heads in",total_tosses,"tosses:",likelihood,"\n"))
## ---- rcode_4end

## ---- rcode_5
source("src/flip_coin_likelihood.r")
likelihood_profile <- flip_coin_likelihood(n = total_tosses,
                                           k = heads_count,
                                           p = seq(0,1,0.001))
## ---- rcode_5b
maxL <- likelihood_profile$p[which(likelihood_profile$likelihood==
        max(likelihood_profile$likelihood))]
CI95 <- likelihood_profile$p[which(likelihood_profile$likelihood>
        exp(log(max(likelihood_profile$likelihood))-1.92))]
CI95 <- CI95[c(1,length(CI95))]
## ---- rcode_5end

plot(x    = likelihood_profile$p,
     y    = likelihood_profile$likelihood,
     xlab = "p",
     ylab = "Likelihood",
     type = "l")
abline(v   = maxL,
       col = "red")
{abline(v   = CI95[1],
       lty = 2,
       col = "red")
abline(v   = CI95[2],
       lty = 2,
       col = "red")}


##############################
# APPROXIMATION 1: SIMULATION
##############################

## ---- rcode_6
toss_simulation <- rbinom(n = 10, size = 1, prob = p_heads)
toss_simulation[which(toss_simulation==1)] <- "H"
toss_simulation[which(toss_simulation==0)] <- "T"  
## ---- rcode_6end
cat(paste("Toss simulation:\n"));(toss_simulation)


## ---- rcode_7
# Note that it prints on screen the first 10 replicates only
likelihood_approx <- flip_coin_likelihood_approx(n     = total_tosses,
                                                 k     = heads_count,
                                                 p     = p_heads,
                                                 rep   = 1000,
                                                 trace = T)
## ---- rcode_7end
cat(paste("Estimated likelihood of p=",p_heads,"given",heads_count,"heads in",
          total_tosses,"tosses:",likelihood_approx$likelihood,"\n"))


## ---- rcode_8
likelihood_profile_approx <-
              flip_coin_likelihood_approx(n     = total_tosses,
                                          k     = heads_count,
                                          p     = seq(0,1,0.01),
                                          rep   = 10000,
                                          trace = F)
## ---- rcode_8end

plot(x    = likelihood_profile_approx$p,
     y    = likelihood_profile_approx$likelihood,
     xlab = "p",
     ylab = "Likelihood",
     lwd  = 2,
     type = "l")
lines(x   = likelihood_profile$p,
      y   = likelihood_profile$likelihood,
      col = "blue")

####################
# BAYESIAN APPROACH            
####################

# EXACT REJECTION ALGORITHM
#---------------------------
# 1. generate p' from prior
# 2. generate D' from simulation with p'
# 3. accept p' if D'=D; return to 1

## ---- rcode_9
sim_parameter_pi1 <- runif(10000,0,1)  # Uninformative prior = rbeta(100000,1,1)
sim_data_pi1      <- rbinom(n    = length(sim_parameter_pi1),
                            size = total_tosses,
                            prob = sim_parameter_pi1)
## ---- rcode_9end



plot(x    = sim_parameter_pi1,
     y    = sim_data_pi1,
     xlab = expression(italic(p)*"'"),
     ylab = expression(italic(D)*"'"))
sim_parameter_kept_pi1 <- sim_parameter_pi1[which(sim_data_pi1==heads_count)]
sim_data_kept_pi1      <- sim_data_pi1[which(sim_data_pi1==heads_count)]
abline(h=heads_count,col="red")
points(sim_parameter_kept_pi1,sim_data_kept_pi1,col="red")

# Estimation of posterior probability
hist(x      = sim_parameter_pi1,
     breaks = seq(0,1,0.02),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,5),
     main   = "",
     xlab   = expression(italic(p)),
     ylab   = "probability density")
hist(x      = sim_parameter_kept_pi1,
     breaks = seq(0,1,0.02),
     col    = rgb(1,0,0,0.5),
     freq   = FALSE,
     add    = TRUE)
box()
lines(x   = seq(0,1,0.001),
      y   = dbeta(x=seq(0,1,0.001),1+heads_count,1+tails_count),
      col = "red",
      lwd = 2)

## ---- rcode_10
p_hat  <- median(sim_parameter_kept_pi1)
p_95CI <- quantile(sim_parameter_kept_pi1,probs=c(0.025,0.975))
## ---- rcode_10end
cat(paste(c("Point estimate of p:",round(p_hat,2),"; 95%CI:",
          round(as.vector(p_95CI),2),"\n"),sep=" ",collapse = " "))



# Let's try another prior

## ---- rcode_11
# Prior believe that the coin is fair
sh1 <- sh2 <- 2
sim_parameter_pi2 <- rbeta(100000,sh1,sh2)
sim_data_pi2      <- rbinom(n    = length(sim_parameter_pi2),
                            size = 10,
                            prob = sim_parameter_pi2)

sim_parameter_kept_pi2 <- sim_parameter_pi2[which(sim_data_pi2==
                                                  heads_count)]
sim_data_kept_pi2      <- sim_data_pi2[which(sim_data_pi2==
                                             heads_count)]
## ---- rcode_11end


hist(x      = sim_parameter_pi2,
     breaks = seq(0,1,0.02),
     col    = "grey",
     freq   = F,
     ylim   = c(0,5),
     main   = "",
     xlab   = expression(italic(p)),
     ylab   = "probability density")
hist(sim_parameter_kept_pi2, breaks=seq(0,1,0.02), col=rgb(0,0,1,0.5), freq=F, add=T)
box()
lines(x   = seq(0,1,0.001),
      y   = dbeta(x=seq(0,1,0.001),sh1+heads_count,sh2+tails_count),
      col = "blue",
      lwd = 2)
lines(x   = seq(0,1,0.001),
      y   = dbeta(x=seq(0,1,0.001),1+heads_count,1+tails_count),
      col = "red",
      lwd = 2)

## ---- rcode_12
p_hat  <- median(sim_parameter_kept_pi2)
p_95CI <- quantile(sim_parameter_kept_pi2,probs=c(0.025,0.975))
## ---- rcode_12end
cat(paste(c("Point estimate of p:",round(p_hat,2),"; 95%CI:",
            round(as.vector(p_95CI),2),"\n"),sep=" ",collapse = " "))




######################################
# APPROXIMATION 2: SUMMARY STATISTICS
######################################

# EXACT REJECTION ALGORITHM ON SUMMARY STATISTICS
#-------------------------------------------------
# 1. generate theta' from prior
# 2. generate D' from simulation with theta'
# 3. calculate SS' from D' 
# 4. accept theta' if SS'=SS; return to 1

# OBSERVED DATA (actually simulated)
## ---- rcode_13
sample_size <- 50
true_theta  <- 10
if (file.exists("data/dataset1.txt")){
  file_removed <- file.remove("data/dataset1.txt")
}
ms(nsam      = sample_size,
   opt       = paste("-t",true_theta),
   temp.file = "data/dataset1.txt")
true_theta0  <- 50
true_theta1  <- 0.5
true_tau     <- 0.5
if (file.exists("data/dataset2.txt")){
  file_removed <- file.remove("data/dataset2.txt")
}
ms(nsam      = sample_size,
   opt       = paste("-t",true_theta0,"-eN",true_tau/true_theta0,true_theta1/true_theta0),
   temp.file = "data/dataset2.txt")
## ---- rcode_13end


## ---- rcode_14_A
source("src/ms.r")
# read file
seq_data         <- ms.inp.multi(sample_size, 1, ms.output.file="data/dataset1.txt")
# number of segregating sites
target1_S        <- S(seq_data)
# nucleotide diversity (mean number of pairwise differences)
target1_pi       <- thetaPi(seq_data) # 
# number of haplotypes
target1_NH       <- NH(seq_data)
# site frequency spectrum
target1_SFS      <- SFS(seq_data)     
# Tajima's D
target1_TajimasD <- tajimaD(seq_data, thetaW(seq_data), target1_pi)
# Fay and Wu H
target1_FayWuH   <- fayWuH(seq_data)
# Fu and Li D
target1_FuLiD    <- fuliD(seq_data, thetaS1(seq_data) )
## ---- rcode_14_Aend
(seq_data)
cat(paste("Data set 1:\n",
          "Sample size:",sample_size,"\n",
          "Number of polymorphic sites:",target1_S,"\n",
          "Number of haplotypes:",target1_NH,"\n"))
colnames(target1_SFS)<-1:(sample_size-1)
barplot(height = target1_SFS/target1_S,
        main   = "Unfolded Site Frequency Spectrum",
        xlab   = "derived allele count in sample",
        ylab   = "Proportion of sites")
box()


## ---- rcode_14_B
seq_data         <- ms.inp.multi(sample_size, 1, ms.output.file="data/dataset2.txt")
target2_S        <- S(seq_data) 
target2_pi       <- thetaPi(seq_data) # 
target2_NH       <- NH(seq_data)
target2_SFS      <- SFS(seq_data)     
target2_TajimasD <- tajimaD(seq_data, thetaW(seq_data), target1_pi)
target2_FayWuH   <- fayWuH(seq_data)
target2_FuLiD    <- fuliD(seq_data, thetaS1(seq_data) )
## ---- rcode_14_Bend
rm(seq_data)
cat(paste("Data set 2:\n",
          "Sample size:",sample_size,"\n",
          "Number of polymorphic sites:",target2_S,"\n",
          "Number of haplotypes:",target2_NH,"\n"))
colnames(target2_SFS)<-1:(sample_size-1)
barplot(height = target2_SFS/target2_S,
        main   = "Unfolded Site Frequency Spectrum",
        xlab   = "derived allele count in sample",
        ylab   = "Proportion of sites")
box()




## ---- rcode_15
# set number of simulations
nsim <- 10000
# Sample parameter theta=4*N*mu from prior 
sim_theta <- 10^runif(nsim,min=-1,max=2) # uniform prior on log10(theta)

if (file.exists("data/abc_sims1.txt")){
  file_removed <- file.remove("data/abc_sims1.txt")
}
ms(nsam       = sample_size, 
   opt        = "-t tbs",
   tbs.matrix = cbind(sim_theta),
   temp.file  = "data/abc_sims1.txt")

msout <- ms.inp.multi(sample_size, ndraws=nsim, ms.output.file="data/abc_sims1.txt")
sim_S        <- S(msout)
sim_pi       <- thetaPi(msout)
sim_NH       <- NH(msout)
sim_SFS      <- SFS(msout)
sim_TajimasD <- tajimaD(msout, thetaW(msout), sim_pi)
sim_FayWuH   <- fayWuH(msout)
sim_FuLiD    <- fuliD(msout, thetaS1(msout) )
## ---- rcode_15end
sim_pi       <- sim_pi[which(sim_S!=0)]
sim_NH       <- sim_NH[which(sim_S!=0)]
sim_SFS      <- sim_SFS[which(sim_S!=0)]
sim_TajimasD <- sim_TajimasD[which(sim_S!=0)]
sim_FayWuH   <- sim_FayWuH[which(sim_S!=0)]
sim_FuLiD    <- sim_FuLiD[which(sim_S!=0)]
sim_theta    <- as.matrix(sim_theta[which(sim_S!=0)])
sim_S        <- sim_S[which(sim_S!=0)]
colnames(sim_theta) <- "theta"




# Set summary statistic to be used
target_SS <- target1_S
sim_SS    <- sim_S

plot(log10(sim_theta),sim_SS,xlab=expression(log[10](theta*"'")),ylab="SS'",ylim=c(0,100))

sim_theta_kept <- sim_theta[which(sim_SS==target_SS)]
sim_SS_kept    <- sim_SS[which(sim_SS==target_SS)]
abline(h=target_SS,col="red")
length(sim_SS_kept)
points(log10(sim_theta_kept),sim_SS_kept,col="red")


# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
hist(log10(sim_theta_kept), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T)
box()


target_SS <- target1_pi
sim_SS    <- sim_pi


plot(log10(sim_theta),sim_SS,xlab=expression(log[10](theta*"'")),ylab="SS'",ylim=c(0,100))

sim_theta_kept <- sim_theta[which(sim_SS==target_SS)]
sim_SS_kept    <- sim_SS[which(sim_SS==target_SS)]
abline(h=target_SS,col="red")
length(sim_SS_kept)
points(log10(sim_theta_kept),sim_SS_kept,col="red")

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,6),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
hist(log10(sim_theta_kept), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T)
box()













########################################
# APPROXIMATION 3: TOLERANCE (DISTANCE)
########################################

# REJECTION ALGORITHM (ABC)
#-------------------------------------------------
# 1. generate theta' from prior
# 2. generate D' from simulation with theta'
# 3. calculate SS' from D' 
# 4. accept theta' if distance(SS',SS) < delta; return to 1

# set tolerance
tolerance <- 0.1

sim_SS <- as.matrix(sim_SS)
colnames(sim_SS) <- "pi"

abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = sim_SS,
                  tol     = tolerance,
                  method  = "rejection")
sim_theta_kept_tol1 <- as.vector(abc_result$unadj.values)
sim_SS_kept_tol1    <- abc_result$ss


plot(log10(sim_theta),sim_SS,xlab=expression(log[10](theta*"'")),ylab="SS'",ylim=c(0,100))
abline(h=target_SS,col="red")
{abline(h=max(abc_result$ss),col="red",lty=2)
 abline(h=min(abc_result$ss),col="red",lty=2)}
points(log10(sim_theta_kept_tol1),sim_SS_kept_tol1,col="red")



# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
hist(log10(sim_theta_kept_tol1), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T)
box()



# Let's invcrease the tolerance
tolerance <- 0.5

abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = as.matrix(sim_SS),
                  tol     = tolerance,
                  method  = "rejection")
sim_theta_kept_tol2 <- as.vector(abc_result$unadj.values)
sim_SS_kept_tol2    <- abc_result$ss


plot(log10(sim_theta),sim_SS,xlab=expression(log[10](theta*"'")),ylab="SS'",ylim=c(0,100))
abline(h=target_SS,col="red")
{abline(h=max(abc_result$ss),col="red",lty=2)
 abline(h=min(abc_result$ss),col="red",lty=2)}
points(log10(sim_theta_kept_tol2),sim_SS_kept_tol2,col="red")



# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log(theta)),
     ylab="probability density")
hist(log10(sim_theta_kept_tol2), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T)
hist(log10(sim_theta_kept_tol1), breaks=seq(-1,2,0.05), col=rgb(0,0,1,0.5), freq=F, add=T)
box()




















########################
# REGRESSION ADJUSTMENT
########################

# set tolerance
tolerance <- 0.5

abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=as.matrix(sim_SS),
                  tol=tolerance,
                  transf = "log",
                  method="loclinear")

sim_theta_kept     <- as.vector(abc_result$unadj.values)
sim_theta_adjusted <- as.vector(abc_result$adj.values)
sim_SS_kept        <- abc_result$ss
local_regression <- lm(log10(sim_theta_kept)~sim_SS_kept, weights=abc_result$weights)

plot(sim_SS_kept,
     abc_result$weights,
     xlab="SS'",
     ylab="weight")
abline(v=target_SS,col="red")

plot(log10(sim_theta),sim_SS,xlab=expression(log[10](theta*"'")),ylab="SS'",ylim=c(0,100))
abline(h=target_SS,col="red")
{abline(h=max(abc_result$ss),col="red",lty=2)
 abline(h=min(abc_result$ss),col="red",lty=2)}
points(log10(sim_theta_kept),sim_SS_kept,col="red")
abline(a=-local_regression$coefficients[1]/local_regression$coefficients[2],
       b=1/local_regression$coefficients[2],
       col="blue",lwd=3)
sim_kept <- 1:10
points(log10(sim_theta_kept)[sim_kept],sim_SS_kept[sim_kept],col="blue")
arrows(x0  = log10(sim_theta_kept)[sim_kept],
       y0  = sim_SS_kept[sim_kept],
       x1  = log10(sim_theta_adjusted)[sim_kept],
       y1  = target_SS,
       col = "blue",
       length = 0.1)
points(log10(sim_theta_adjusted),array(target_SS,length(sim_theta_adjusted)),col="blue")


# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log(theta)),
     ylab="probability density")
hist(log10(sim_theta_kept), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T )
wtd.hist(log10(sim_theta_adjusted), breaks=seq(-1,2,0.05), col=rgb(0,1,0,0.5), freq=F, add=T, weight=abc_result$weights)
hist(log10(sim_theta_kept_tol1), breaks=seq(-1,2,0.05), col=rgb(0,0,1,0.5), freq=F, add=T)
box()



















####################################
# INFORMATION ON SUMMARY STATISTICS
####################################

# set tolerance
tolerance <- 0.1

target_SS <- target1_pi
sim_SS <- as.matrix(sim_pi)
colnames(sim_SS) <- "pi"


abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=as.matrix(sim_SS),
                  tol=tolerance,
                  transf = "log",
                  method="loclinear")

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
wtd.hist(log10(sim_theta_adjusted_pi), breaks=seq(-2,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, weight=sim_weights_pi)
box()

# Tajima's D
target_SS <- target1_TajimasD
sim_SS    <- cbind(sim_TajimasD)
colnames(sim_SS) <- "Tajima's D"

abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=as.matrix(sim_SS),
                  tol=tolerance,
                  transf = "log",
                  method="loclinear")

sim_theta_adjusted_D <- as.vector(abc_result$adj.values)
sim_weights_D        <- abc_result$weights

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted_D), breaks=seq(-2,3,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, weight=sim_weights_D)
box()
wtd.hist(log10(sim_theta_adjusted_pi), breaks=seq(-2,3,0.05), col=rgb(0,0,1,0.5), freq=F, add=T, weight=sim_weights_pi)


sim_theta_kept     <- as.vector(abc_result$unadj.values)
sim_theta_adjusted <- as.vector(abc_result$adj.values)
sim_SS_kept        <- abc_result$ss
local_regression <- lm(log10(sim_theta_kept)~sim_SS_kept, weights=sim_weights_D)

plot(log10(sim_theta),sim_SS,xlab=expression(log(theta*"'")),ylab="SS'",ylim=c(-4,4))
abline(h=target_SS,col="red")
{abline(h=max(abc_result$ss),col="red",lty=2)
 abline(h=min(abc_result$ss),col="red",lty=2)}
points(log10(sim_theta_kept),sim_SS_kept,col="red")
abline(a=-local_regression$coefficients[1]/local_regression$coefficients[2],
       b=1/local_regression$coefficients[2],
       col="blue",lwd=3)
points(log10(sim_theta_adjusted),array(target_SS,length(sim_theta_adjusted)),col="blue")

cor(sim_TajimasD,as.vector(sim_theta))
cor(sim_pi,as.vector(sim_theta))

target_SS <- cbind(target1_pi, 
                    target1_TajimasD)
sim_SS     <- cbind(sim_pi, 
                    sim_TajimasD)
colnames(target_SS) <- colnames(sim_SS) <- c("pi","TD")

abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=sim_SS,
                  tol=tolerance,
                  transf = "log",
                  method="loclinear")

sim_theta_adjusted_piD <- as.vector(abc_result$adj.values)
sim_weights_piD        <- abc_result$weights

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted_piD), breaks=seq(-2,3,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, weight=sim_weights_piD)
box()
wtd.hist(log10(sim_theta_adjusted_pi), breaks=seq(-2,2,0.05), col=rgb(0,0,1,0.5), freq=F, add=T, weight=sim_weights_pi)



target1_noise <- 76
sim_noise     <- 10^runif(length(sim_pi),-2,3)

target_SS <- cbind(target1_noise, 
                   target1_pi)
sim_SS     <- cbind(sim_noise, 
                    sim_pi)
colnames(target_SS) <- colnames(sim_SS) <- c("noise","pi")


abc_result <- abc(target=target_SS,
                   param=sim_theta,
                   sumstat=as.matrix(sim_SS),
                   tol=tolerance,
                   transf = "log",
                   method="loclinear")

sim_theta_adjusted_pi_noise <- as.vector(abc_result$adj.values)
sim_weights_pi_noise        <- abc_result$weights

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted_pi_noise), breaks=seq(-2,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, weight=sim_weights_pi_noise)
wtd.hist(log10(sim_theta_adjusted_pi), breaks=seq(-1,2,0.05), col=rgb(0,1,0,0.5), freq=F, add=T, weight=sim_weights_pi)
#wtd.hist(log10(sim_theta_adjusted_piD), breaks=seq(-2,2,0.05), col=rgb(0,0,1,0.5), freq=F, add=T, weight=sim_weights_piD)
box()

















######################
# BEWARE OF THE PRIOR
######################

# prior on mutation rate (3e-8 - 1e-6) and population size (1-30 million)

sim_mu <- 10^runif(nsim,min=-7.5,max=-6)
sim_Ne <- 10^runif(nsim,min=6,max=7.5)
sim_theta_composite <- 4*sim_Ne*sim_mu

hist(log10(sim_theta_composite),
     breaks=seq(-2,3,0.05),
     col="grey",
     freq=F,
     ylim=c(0,1),
     xlim=c(-1,2),
     main="",
     xlab=expression(log[10](theta)),
     ylab="probability density")
box()








target1_SS <- cbind(target1_pi, 
                    target1_TajimasD)
target2_SS <- cbind(target2_pi, 
                    target2_TajimasD)
sim_SS     <- cbind(sim_pi,
                    sim_TajimasD)
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
wtd.hist(log10(sim_theta_adjusted_1), breaks=seq(-2,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T, weight=sim_weights_1)
wtd.hist(log10(sim_theta_adjusted_2), breaks=seq(-2,2,0.05), col=rgb(0,1,0,0.5), freq=F, add=T, weight=sim_weights_2)
box()








##################
# QUALITY CONTROL
##################

# 1. Assessment of the method per se (validation through simulted data)

# cross validation

cross_validation_result <- cv4abc(param   = sim_theta,
                                  sumstat = as.matrix(sim_SS),
                                  abc.out = abc_result1,
                                  nval    = 100,
                                  tols    = 0.1)

plot(x    = log10(cross_validation_result$true$theta),
     y    = log10(cross_validation_result$estim$tol0.1),
     xlab = expression("true "*log(theta)),
     ylab = expression(log(hat(theta))),
     xlim = c(-1,2), ylim = c(-1,2))
abline(a=0,b=1)


# 2. Assessment of the applicability of the method to the observed data

# Are the simulations producing pseudo-data similar to target data?
# Summary statistics used in the analysis
plot(sim_S,sim_pi,xlab="S",ylab=expression(pi),log="xy")
points(target1_S,target1_pi,col="red",cex=3,pch="*")
points(target2_S,target2_pi,col="green",cex=3,pch="*")

plot(sim_NH,sim_pi,xlab="Number of Haplotypes",ylab=expression(pi),log="xy")
points(target1_NH,target1_pi,col="red",cex=3,pch="*")
points(target2_NH,target2_pi,col="green",cex=3,pch="*")

plot(sim_FuLiD,sim_pi,xlab="Fu and Li's D",ylab=expression(pi),log="y",xlim=c(-6,2))
points(target1_FuLiD,target1_pi,col="red",cex=3,pch="*")
points(target2_FuLiD,target2_pi,col="green",cex=3,pch="*")

sim_SS <- cbind(sim_S,
                sim_pi, 
                sim_NH,
                sim_TajimasD,
                sim_FayWuH,
                sim_FuLiD)
target1_SS <- cbind(target1_S,
                    target1_pi, 
                    target1_NH,
                    target1_TajimasD,
                    target1_FayWuH,
                    target1_FuLiD)
target2_SS <- cbind(target2_S,
                    target2_pi, 
                    target2_NH,
                    target2_TajimasD,
                    target2_FayWuH,
                    target2_FuLiD)
colnames(sim_SS) <- colnames(target1_SS) <-
                    colnames(target2_SS) <- c("S","pi","NH","TD","FWH","FLD")

PCA_stats   <- princomp(~S+pi+NH+TD+FWH+FLD,as.data.frame(sim_SS)) 
PCA_target1 <- predict(PCA_stats,as.data.frame(target1_SS))
PCA_target2 <- predict(PCA_stats,as.data.frame(target2_SS))

plot(PCA_stats$scores[,1], PCA_stats$scores[,2],
     xlab="PC1",ylab="PC2")
points(PCA_target1[1],PCA_target1[2],col="red",cex=3,pch="*")
points(PCA_target2[1],PCA_target2[2],col="green",cex=3,pch="*")
plot(PCA_stats$scores[,3], PCA_stats$scores[,4],
     xlab="PC3",ylab="PC4")
points(PCA_target1[3],PCA_target1[4],col="red",cex=3,pch="*")
points(PCA_target2[3],PCA_target2[4],col="green",cex=3,pch="*")

plot(PCA_stats$scores[,5], PCA_stats$scores[,6],
     xlab="PC5",ylab="PC6",ylim=c(-5,3))
points(PCA_target1[5],PCA_target1[6],col="red",cex=3,pch="*")
points(PCA_target2[5],PCA_target2[6],col="green",cex=3,pch="*")


# Simulations from the posterior probability of the parameter
# Sample parameter theta=4*N*mu from posterior 
posterior_theta <- sample(sim_theta_adjusted1,size=1000,replace=T,prob=sim_weights1)

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
posterior_TajimasD <- tajimaD(msout, thetaW(msout), sim_pi)
posterior_FayWuH   <- fayWuH(msout)
posterior_FuLiD    <- fuliD(msout, thetaS1(msout) )

hist(posterior_NH, breaks=seq(0,40,1), col="grey", freq=F, xlab="Number of haplotypes", main="")
abline(v=target1_NH,col="red",lwd=3)
box()
hist(posterior_FuLiD, breaks=seq(-5,4,0.2), col="grey", freq=F, xlab="Fu and Li's D", main="")
abline(v=target1_FuLiD,col="red",lwd=3)
box()



posterior_theta <- sample(sim_theta_adjusted2,size=1000,replace=T,prob=sim_weights2)

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
posterior_TajimasD <- tajimaD(msout, thetaW(msout), sim_pi)
posterior_FayWuH   <- fayWuH(msout)
posterior_FuLiD    <- fuliD(msout, thetaS1(msout) )

hist(posterior_NH, breaks=seq(0,40,1), col="grey", freq=F, xlab="Number of haplotypes", main="")
abline(v=target2_NH,col="red",lwd=3)
box()
hist(posterior_FuLiD, breaks=seq(-5,5,0.2), col="grey", freq=F, xlab="Fu and Li's D", main="")
abline(v=target2_FuLiD,col="red",lwd=3)
box()























######################
# MODEL CHOICE IN ABC
######################

# set number of simulations
nsim <- 20000

# Consider two alternative models: (1) constant size and (2) population expansion
# Sample models from prior
sim_model <- sample(c("C","V"),size=nsim,replace=T)

# Sample parameter theta=4*N*mu from prior for model (1) constant size
sim_theta0 <- array(NA,nsim)
sim_theta0[sim_model=="C"] <- 10^runif(length(which(sim_model=="C")),min=-1,max=2)

# Sample parameter theta0=4*N0*mu, theta1=4*N1*mu, tau=t*mu from prior for model (2) population expansion
sim_theta0[sim_model=="V"] <- 10^runif(length(which(sim_model=="V")),min=-1,max=2) 
sim_theta1 <- array(NA,nsim)
sim_theta1[sim_model=="V"] <- 10^runif(length(which(sim_model=="V")),min=-1,max=2) 
sim_theta1[sim_model=="C"] <- sim_theta0[sim_model=="C"] 
sim_tau <- array(NA,nsim)
sim_tau[sim_model=="V"] <- runif(length(which(sim_model=="V")),min=0,max=2) 
sim_tau[sim_model=="C"] <- 0 

ref_table_params <- cbind(sim_theta0,sim_theta1,sim_tau)
colnames(ref_table_params) <- c("theta0","theta1","tau")

tbs_params <-  cbind(ref_table_params[,"theta0"],
                     ref_table_params[,"tau"]/ref_table_params[,"theta0"],
                     ref_table_params[,"theta1"]/ref_table_params[,"theta0"])

if (file.exists("data/abc_sims2.txt")){
  file_removed <- file.remove("data/abc_sims2.txt")
}
ms(nsam       = sample_size, 
   opt        = "-t tbs -eN tbs tbs",
   tbs.matrix = tbs_params,
   temp.file  = "data/abc_sims2.txt")

msout <- ms.inp.multi(sample_size, nsim, ms.output.file="data/abc_sims2.txt")
# simulated summary statistics from posterior
sim_S        <- S(msout)
sim_pi       <- thetaPi(msout)
sim_NH       <- NH(msout)
sim_TajimasD <- tajimaD(msout, thetaW(msout), sim_pi)
sim_FayWuH   <- fayWuH(msout)
sim_FuLiD    <- fuliD(msout, thetaS1(msout) )

sim_pi       <- sim_pi[which(sim_S!=0)]
sim_NH       <- sim_NH[which(sim_S!=0)]
sim_SFS      <- sim_SFS[which(sim_S!=0)]
sim_TajimasD <- sim_TajimasD[which(sim_S!=0)]
sim_FayWuH   <- sim_FayWuH[which(sim_S!=0)]
sim_FuLiD    <- sim_FuLiD[which(sim_S!=0)]
sim_model    <- sim_model[which(sim_S!=0)]
ref_table_params <- ref_table_params[which(sim_S!=0),]
sim_S        <- sim_S[which(sim_S!=0)]

modelC <- which(sim_model=="C")
modelV <- which(sim_model=="V")

# Model (and prior) quality control
plot(sim_S[modelV],sim_pi[modelV],xlab="S",ylab=expression(pi),log="xy")
points(sim_S[modelC],sim_pi[modelC],col="grey")
points(target1_S,target1_pi,col="red",cex=3,pch="*")
points(target2_S,target2_pi,col="green",cex=3,pch="*")
plot(sim_FuLiD[modelV],sim_pi[modelV],xlab="Fu and Li's D",ylab=expression(pi),log="y",xlim=c(-6,2))
points(sim_FuLiD[modelC],sim_pi[modelC],col="grey")
points(target1_FuLiD,target1_pi,col="red",cex=3,pch="*")
points(target2_FuLiD,target2_pi,col="green",cex=3,pch="*")

# Model choice
tolerance <- 0.1
target1_SS <- cbind(target1_S,target1_pi,target1_NH,target1_TajimasD,target1_FayWuH,target1_FuLiD)
target2_SS <- cbind(target2_S,target2_pi,target2_NH,target2_TajimasD,target2_FayWuH,target2_FuLiD)
sim_SS    <- cbind(sim_S,sim_pi,sim_NH,sim_TajimasD,sim_FayWuH,sim_FuLiD)
colnames(target1_SS) <- colnames(target2_SS) <- 
                        colnames(sim_SS) <- c("S","pi","NH","TD","FWH","FLD")

abc_model_choice1 <- postpr(target=target1_SS,
                            index=sim_model,
                            sumstat=sim_SS,
                            tol=tolerance,
                            method="mnlogistic")
                           
abc_model_choice1$pred                           

abc_model_choice2 <- postpr(target=target2_SS,
                            index=sim_model,
                            sumstat=sim_SS,
                            tol=tolerance,
                            method="mnlogistic")

abc_model_choice2$pred                           






################
# ABCRF
################

# CART
#------

theta  <- ref_table_params[modelV,"theta1"]
TajD   <- sim_SS[modelV,"TD"]
pi     <- sim_SS[modelV,"pi"]
ref_table <- data.frame(theta,TajD,pi)

regression_tree <- tree(theta ~ TajD + pi,
                        data=ref_table)

plot(regression_tree)
text(regression_tree,cex=0.75)

plot(ref_table$TajD,
     ref_table$pi,
     xlab="Tajima's D",
     ylab=expression(pi),
     col=terrain.colors(5)[theta],
     pch=20)
partition.tree(regression_tree,
               ordvars=c("TajD","pi"),
               add=T,cex=1.5,col="blue")


TajD   <- sim_SS[,"TD"]
FuLiD  <- sim_SS[,"FLD"]
ref_table <- data.frame(sim_model,TajD,FuLiD)
classification_tree <- tree(sim_model ~ TajD + FuLiD,
                       data=ref_table)

plot(classification_tree)
text(classification_tree,cex=0.75)

plot(ref_table$TajD[which(sim_model=="V")],
     ref_table$FuLiD[which(sim_model=="V")],
     xlab="Tajima's D",
     ylab="Fu and Li's D",
     pch=20)
points(ref_table$TajD[which(sim_model=="C")],
       ref_table$FuLiD[which(sim_model=="C")],
       col="grey",pch=20)
partition.tree(classification_tree,
               ordvars=c("TajD","FuLiD"),
               add=T,cex=3,col="red")



# Random Forest
#---------------

# Explain Bagging (note out-of-bag)

theta     <- ref_table_params[modelC,"theta0"]
pi        <- sim_SS[modelC,"pi"]

ref_table <- data.frame(theta,pi)
regression_tree <- tree(log10(theta) ~ pi,
                        data=ref_table)
plot(regression_tree)
text(regression_tree,cex=0.75)

plot(ref_table$pi,
     log10(ref_table$theta),
     xlab=expression(pi),
     ylab=expression(log[10]*theta),
     pch=20,log="x")
partition.tree(regression_tree,
               ordvars=c("pi","theta"),
               add=T,cex=1.5,col="blue",lwd=2)

for (i in 1:100){
  random_sample <- sample(length(theta),size=100,replace=T)
  ref_table_random_sample <- ref_table[random_sample,]
  regression_tree_random_sample <- tree(log10(theta) ~ pi,
                                        data=ref_table_random_sample)
  partition.tree(regression_tree_random_sample,
                 ordvars=c("pi","theta"),
                 add=T,cex=1.5,col="red",lwd=1)
  
}




# DO NOT USE FWH statistics, abcrf crashes, don't know why

ref_table <- data.frame(sim_model,sim_SS[,-5])
model_RF <- abcrf(formula = sim_model~.,
                  data    = ref_table,
                  lda     = F,
                  ntree   = 500,
                  paral   = T)
plot(model_RF,
     training=ref_table)
model_RF$prior.err

# error and number of trees, it does not work for this dataset, don't know why
#err.abcrf(model_RF,
#          training=ref_table,
#          paral=T)



predict(object=model_RF$model.rf,
        data=as.data.frame(rbind(target1_SS,target2_SS)))$predictions
rowSums(predict(object=model_RF$model.rf,predict.all = TRUE,
        data=as.data.frame(rbind(target1_SS,target2_SS)))$predictions==1)


model_selection_result_RF <- predict(object         = model_RF,
                                     obs            = as.data.frame(rbind(target1_SS,target2_SS)),
                                     training       = ref_table,
                                     ntree          = 500,
                                     paral          = T,
                                     paral.predict  = T)
(model_selection_result_RF)





log10theta <- log10(ref_table_params[modelV,"theta1"])
ref_table <- data.frame(log10theta,sim_SS[modelV,])
RFmodel_theta <- regAbcrf(formula = log10theta~.,
                           data    = ref_table,
                           ntree   = 1000,
                           paral   = T)


plot(RFmodel_theta)
#err.regAbcrf(object   = RFmodel_theta,
#             training = ref_table,
#             paral    = T)
posterior_theta_RF <- predict(object    = RFmodel_theta,
                              obs       = as.data.frame(target1_SS),
                              training  = ref_table,
                              paral     = T)
(posterior_theta_RF)
densityPlot(object    = RFmodel_theta,
            obs       = as.data.frame(target1_SS),
            training  = ref_table,
            main      = expression(log[10]*theta),
            paral     = T)
lines(density(log10theta), col="grey")



