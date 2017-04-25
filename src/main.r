################################################
# Workshop on Approximate Bayesian Computation 
#     from rejection to random forests
################################################

# Miguel Navascués
# miguel.navascues@inra.fr
# INRA

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

p_combo <- p_tails * p_heads * p_heads * 
           p_tails * p_heads * p_heads * 
           p_tails * p_heads * p_heads * p_heads
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
abline(v   = CI95[1],
       lty = 2,
       col = "red")
abline(v   = CI95[2],
       lty = 2,
       col = "red")


##############################
# APPROXIMATION 1: SIMULATION
##############################

## ---- rcode_6
toss_simulation <- rbinom(n = 10, size = 1, prob = 0.5)
toss_simulation[which(toss_simulation==1)] <- "H"
toss_simulation[which(toss_simulation==0)] <- "T"  
## ---- rcode_6end
cat(paste("Toss simulation:\n"));(toss_simulation)


## ---- rcode_7
# Note that it prints on screen the first 10 replicates only
likelihood_approx <- flip_coin_likelihood_approx(n     = total_tosses,
                                                 k     = heads_count,
                                                 p     = 0.5,
                                                 rep   = 1000,
                                                 trace = T)
## ---- rcode_7end
cat(paste("Estimated likelihood of p=0.5 given",heads_count,"heads in",
          total_tosses,"tosses:",likelihood_approx$likelihood,"\n"))


## ---- rcode_8
likelihood_profile_approx <-
              flip_coin_likelihood_approx(n     = total_tosses,
                                          k     = heads_count,
                                          p     = seq(0,1,0.01),
                                          rep   = 1000,
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
sim_parameter_pi1 <- runif(10000,0,1)  # Uninformative prior
sim_data_pi1      <- rbinom(n    = length(sim_parameter_pi1),
                            size = 10,
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
     freq   = F,
     ylim   = c(0,5),
     main   = "",
     xlab   = expression(italic(p)),
     ylab   = "probability density")
hist(sim_parameter_kept_pi1, breaks=seq(0,1,0.02), col=rgb(1,0,0,0.5), freq=F, add=T)
box()

## ---- rcode_10
p_hat  <- median(sim_parameter_kept_pi1)
p_95CI <- quantile(sim_parameter_kept_pi1,probs=c(0.025,0.975))
## ---- rcode_10end

# Let's try another prior

## ---- rcode_11
# Prior believe that the coin is fair
sim_parameter_pi2 <- rbeta(10000,6,6)
sim_data_pi2      <- rbinom(n    = length(sim_parameter_pi2),
                            size = 10,
                            prob = sim_parameter_pi2)

sim_parameter_kept_pi2 <- sim_parameter_pi2[which(sim_data_pi2==
                                                  heads_count)]
sim_data_kept_pi2      <- sim_data_pi2[which(sim_data_pi2==
                                             heads_count)]
## ---- rcode_11end


hist(x      = sim_parameter_pi1,
     breaks = seq(0,1,0.02),
     col    = "grey",
     freq   = F,
     ylim   = c(0,5),
     main   = "",
     xlab   = expression(italic(p)),
     ylab   = "probability density")
hist(sim_parameter_kept_pi1, breaks=seq(0,1,0.02), col=rgb(1,0,0,0.5), freq=F, add=T)
hist(sim_parameter_pi2, breaks=seq(0,1,0.02), col=rgb(0,0,0,0.3), freq=F, add=T)
hist(sim_parameter_kept_pi2, breaks=seq(0,1,0.02), col=rgb(0,0,1,0.3), freq=F, add=T)
box()

## ---- rcode_12
p_hat  <- median(sim_parameter_kept_pi2)
p_95CI <- quantile(sim_parameter_kept_pi2,probs=c(0.025,0.975))
## ---- rcode_12end




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
require("phyclust", quietly = TRUE)       # for a version of Hudson ms (coalescent simulator)
sample_size <- 50
true_theta  <- 10
if (file.exists("data/dataset1.txt")){
  file_removed <- file.remove("data/dataset1.txt")
}
ms(nsam      = sample_size,
   opt       = paste("-t",true_theta),
   temp.file = "data/dataset1.txt")
detach("package:phyclust", unload=TRUE)
detach("package:ape", unload=TRUE)
## ---- rcode_13end


## ---- rcode_14
source("src/ms.r")
# read file
seq_data         <- ms.inp.multi(sample_size, 1, ms.output.file="data/dataset1.txt")
# number of segregating sites
target1_S        <- seq_data[[1]]$nss 
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
## ---- rcode_14end

## ---- rcode_15
require("phyclust", quietly = TRUE)       

# set number of simulations
nsim <- 10000
# Sample parameter theta=4*N*mu from prior 
sim_theta <- 10^runif(nsim,min=-1,max=2) # uniform prior on log10(theta)
# simulated summary statistics
sim_S             <- array(NA,nsim)
sim_pi            <- array(NA,nsim)
sim_NH            <- array(NA,nsim)
sim_SFS           <- matrix(NA,nsim,sample_size-1)
sim_TajimasD      <- array(NA,nsim)
sim_FayWuH        <- array(NA,nsim)
sim_FuLiD         <- array(NA,nsim)
for (sim in 1:nsim){
  ms(nsam=sample_size, nreps=1, opt=paste("-t",sim_theta[sim]), temp.file ="temp_ms_file.txt")
  msout <- ms.inp.multi(sample_size, 1, ms.output.file="temp_ms_file.txt")
  system("rm temp_ms_file.txt")
  sim_S[sim]             <- msout[[1]]$nss
  sim_pi[sim]            <- thetaPi(msout)
  sim_NH[sim]            <- NH(msout)
  sim_SFS[sim,]          <- SFS(msout)
  sim_TajimasD[sim]      <- tajimaD(msout, thetaW(msout), sim_pi[sim])
  sim_FayWuH[sim]        <- fayWuH(msout)
  sim_FuLiD[sim]         <- fuliD(msout, thetaS1(msout) )
}
## ---- rcode_15end



# Set summary statistic to be used
target_SS <- target_S
sim_SS    <- sim_S

# target_SS <- target_pi
# sim_SS    <- sim_pi


plot(log10(sim_theta),sim_SS,xlab=expression(log(theta*"'")),ylab="SS'",ylim=c(0,100))

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
     xlab=expression(log(theta)),
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





# SETTINGS

# install.packages(c("phyclust","abc","weights"))
# load required packages
library("abc")            # for approximate Bayesian computation
library("abcrf")            # for approximate Bayesian computation
library("weights")        # for weighted histograms



# set tolerance
tolerance <- 0.1
# tolerance <- 0.5
# sim_theta_kept_previous <- sim_theta_kept
# tolerance <- 0.99
# sim_theta_kept_previous2 <- sim_theta_kept


abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=as.matrix(sim_SS),
                  tol=tolerance,
                  transf="log",
                  method="rejection")

sim_theta_kept <- as.vector(abc_result$unadj.values)
sim_SS_kept    <- abc_result$ss


plot(log10(sim_theta),sim_SS,xlab=expression(log(theta*"'")),ylab="SS'",ylim=c(0,100))
abline(h=target_SS,col="red")
abline(h=max(abc_result$ss),col="red",lty=2)
abline(h=min(abc_result$ss),col="red",lty=2)
points(log10(sim_theta_kept),sim_SS_kept,col="red")



# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log(theta)),
     ylab="probability density")
hist(log10(sim_theta_kept), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.5), freq=F, add=T)
# hist(log10(sim_theta_kept_previous), breaks=seq(-1,2,0.05), col=rgb(0,0,1,0.5), freq=F, add=T)
# hist(log10(sim_theta_kept_previous2), breaks=seq(-1,2,0.05), col=rgb(0,1,0,0.5), freq=F, add=T)
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

plot(log10(sim_theta),sim_SS,xlab=expression(log(theta*"'")),ylab="SS'",ylim=c(0,100))
abline(h=target_SS,col="red")
abline(h=max(abc_result$ss),col="red",lty=2)
abline(h=min(abc_result$ss),col="red",lty=2)
points(log10(sim_theta_kept),sim_SS_kept,col="red")
abline(a=-local_regression$coefficients[1]/local_regression$coefficients[2],
       b=1/local_regression$coefficients[2],
       col="orange",lwd=3)
points(log10(sim_theta_adjusted),array(target_SS,length(sim_theta_adjusted)),col="orange")


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
wtd.hist(log10(sim_theta_adjusted), breaks=seq(-1,2,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=abc_result$weights)
# hist(log10(sim_theta_kept_previous), breaks=seq(-1,2,0.05), col=rgb(0,0,1,0.5), freq=F, add=T)
box()



















#################################################################
# SUFFICIENT SUMMARY STATISTICS (INFORMATIVE SUMMARY STATISTICS)
#################################################################

# set tolerance
tolerance <- 0.2

# Set summary statistic to be used: pi
target_SS <- target_pi
sim_SS    <- sim_pi
# Tajima's D
# target_SS <- target_TajimasD
# sim_SS    <- sim_TajimasD

abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=as.matrix(sim_SS),
                  tol=tolerance,
                  transf = "log",
                  method="loclinear")

sim_theta_adjusted <- as.vector(abc_result$adj.values)
sim_weights        <- abc_result$weights

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log(theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
#wtd.hist(log10(sim_theta_adjusted), breaks=seq(-2,3,0.05), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
box()

# sim_theta_kept     <- as.vector(abc_result$unadj.values)
# sim_theta_adjusted <- as.vector(abc_result$adj.values)
# sim_SS_kept        <- abc_result$ss
# local_regression <- lm(log10(sim_theta_kept)~sim_SS_kept, weights=sim_weights)

# plot(log10(sim_theta),sim_SS,xlab=expression(log(theta*"'")),ylab="SS'",ylim=c(-4,4))
# abline(h=target_SS,col="red")
# abline(h=max(abc_result$ss),col="red",lty=2)
# abline(h=min(abc_result$ss),col="red",lty=2)
# points(log10(sim_theta_kept),sim_SS_kept,col="red")
# abline(a=-local_regression$coefficients[1]/local_regression$coefficients[2],
#        b=1/local_regression$coefficients[2],
#        col="orange",lwd=3)
# points(log10(sim_theta_adjusted),array(target_SS,length(sim_theta_adjusted)),col="orange")














##########################
# CURSE OF DIMENSIONALITY
##########################
# set tolerance
tolerance <- 0.2

# Set summary statistic to be used: Site Frequency Spectrum
target_SS <- target_SFS
sim_SS    <- sim_SFS
# Number of polymorphic sites
# target_SS <- target_S
# sim_SS    <- sim_S
# sim_theta_adjusted_previous <- sim_theta_adjusted
# sim_weights_previous        <- sim_weights

abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=as.matrix(sim_SS),
                  tol=tolerance,
                  transf = "log",
                  method="loclinear")

sim_theta_adjusted <- as.vector(abc_result$adj.values)
sim_weights        <- abc_result$weights

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log(theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
# wtd.hist(log10(sim_theta_adjusted), breaks=seq(-2,3,0.05), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
# wtd.hist(log10(sim_theta_adjusted_previous), breaks=seq(-2,3,0.05), col=rgb(1,0.5,0,0.7), freq=F, add=T, weight=sim_weights_previous)
box()














##################
# QUALITY CONTROL
##################

# 1. Assessment of the method per se (validation through simulted data)

# 2. Assessment of the applicability of the method to the observed data

# set tolerance
tolerance <- 0.1

# Set summary statistic to be used: Site Frequency Spectrum
target_SS <- cbind(target_S,target_pi)
sim_SS    <- cbind(sim_S,sim_pi)
dimnames(target_SS)[[2]]<-dimnames(sim_SS)[[2]]<-c("S","pi")


abc_result <- abc(target=target_SS,
                  param=sim_theta,
                  sumstat=as.matrix(sim_SS),
                  tol=tolerance,
                  transf = "log",
                  method="loclinear")

sim_theta_adjusted <- as.vector(abc_result$adj.values)
sim_weights        <- abc_result$weights

# Estimation of posterior probability
hist(log10(sim_theta),
     breaks=seq(-1,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,4),
     main="",
     xlab=expression(log(theta)),
     ylab="probability density")
wtd.hist(log10(sim_theta_adjusted), breaks=seq(-1,2,0.05), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
box()


# Are the simulations producing pseudo-data similar to target data?
# Summary statistics used in the analysis
plot(sim_S,sim_pi,xlab="S",ylab=expression(pi),log="xy")
points(target_S,target_pi,col="red",cex=3,pch="*")


# Simulations from point estimate of the parameter
theta_point_estimate <- summary(abc_result)[3]

ms(nsam=sample_size, nreps=1000, opt=paste("-t",theta_point_estimate), temp.file ="data.txt")
pred_seq_data <- ms.inp.multi(sample_size, 1000, ms.output.file="data.txt")
#predicted_S             <- S(pred_seq_data)              # number of segregating sites
predicted_pi            <- thetaPi(pred_seq_data)              # nucleotide diversity (mean number of pairwise differences)
predicted_NH            <- NH(pred_seq_data)
#predicted_SFS           <- SFS(pred_seq_data)                  # site frequency spectrum
predicted_TajimasD      <- tajimaD(pred_seq_data, thetaW(pred_seq_data), predicted_pi)
#predicted_FayWuH        <- fayWuH(pred_seq_data)
#predicted_FuLiD         <- fuliD(pred_seq_data, thetaS1(pred_seq_data) )

hist(predicted_NH, breaks=seq(0,40,1), col="grey", freq=F, xlab="Number of haplotypes", main="")
abline(v=target_NH,col="red",lwd=3)
box()
hist(predicted_TajimasD, breaks=seq(-4,4,0.2), col="grey", freq=F, xlab="Tajima's D", main="")
abline(v=target_TajimasD,col="red",lwd=3)
box()

# Simulations from the posterior probability of the parameter
# Sample parameter theta=4*N*mu from posterior 
posterior_theta <- sample(sim_theta_adjusted,size=1000,replace=T,prob=sim_weights)

# simulated summary statistics
#predicted_S             <- array(NA,1000)
predicted_pi            <- array(NA,1000)
predicted_NH            <- array(NA,1000)
#predicted_SFS           <- matrix(NA,1000,sample_size-1)
predicted_TajimasD      <- array(NA,1000)
#predicted_FayWuH        <- array(NA,1000)
#predicted_FuLiD         <- array(NA,1000)
for (sim in 1:1000){
  ms(nsam=sample_size, nreps=1, opt=paste("-t",posterior_theta[sim]), temp.file ="temp_ms_file.txt")
  pred_seq_data <- ms.inp.multi(sample_size, 1, ms.output.file="temp_ms_file.txt")
  system("rm temp_ms_file.txt")
  #predicted_S[sim]             <- pred_seq_data[[1]]$nss
  predicted_pi[sim]            <- thetaPi(pred_seq_data)
  predicted_NH[sim]            <- NH(pred_seq_data)
  #predicted_SFS[sim,]          <- SFS(pred_seq_data)
  predicted_TajimasD[sim]      <- tajimaD(pred_seq_data, thetaW(pred_seq_data), predicted_pi[sim])
  #predicted_FayWuH[sim]        <- fayWuH(pred_seq_data)
  #predicted_FuLiD[sim]         <- fuliD(pred_seq_data, thetaS1(pred_seq_data) )
}

hist(predicted_NH, breaks=seq(0,40,1), col="grey", freq=F, xlab="Number of haplotypes", main="")
abline(v=target_NH,col="red",lwd=3)
box()
hist(predicted_TajimasD, breaks=seq(-4,4,0.2), col="grey", freq=F, xlab="Tajima's D", main="")
abline(v=target_TajimasD,col="red",lwd=3)
box()

# NEW OBSERVED DATA (actually simulated)
# theta0=4*N0*mu, theta1=4*N1*mu, tau=t*mu
sample_size  <- 50
true_theta0  <- 50
true_theta1  <- 0.5
true_tau     <- 0.5
ms(nsam=sample_size,nreps=1,
   opt=paste("-t",true_theta0,"-eN",true_tau/true_theta0,true_theta1/true_theta0),
   temp.file ="data.txt")
seq_data <- ms.inp.multi(sample_size, 1, ms.output.file="data.txt")
target_S             <- seq_data[[1]]$nss              # number of segregating sites
target_pi            <- thetaPi(seq_data)              # nucleotide diversity (mean number of pairwise differences)
target_NH            <- NH(seq_data)
target_SFS           <- SFS(seq_data)                  # site frequency spectrum
target_TajimasD      <- tajimaD(seq_data, thetaW(seq_data), target_pi)
target_FayWuH        <- fayWuH(seq_data)
target_FuLiD         <- fuliD(seq_data, thetaS1(seq_data) )

























######################
# MODEL CHOICE IN ABC
######################

# set number of simulations
nsim <- 20000

# Consider two alternative models: (1) constant size and (2) population expansion
# Sample models from prior
sim_model <- sample(2,size=nsim,replace=T)

# Sample parameter theta=4*N*mu from prior for model (1) constant size
sim_theta <- array(NA,nsim)
sim_theta[sim_model==1] <- 10^runif(length(which(sim_model==1)),min=-1,max=2)

# Sample parameter theta0=4*N0*mu, theta1=4*N1*mu, tau=t*mu from prior for model (2) population expansion
sim_theta0 <- array(NA,nsim)
sim_theta0[sim_model==2] <- 10^runif(length(which(sim_model==2)),min=-1,max=2) 
sim_theta1 <- array(NA,nsim)
sim_theta1[sim_model==2] <- 10^runif(length(which(sim_model==2)),min=-1,max=2) 
sim_tau <- array(NA,nsim)
sim_tau[sim_model==2] <- runif(length(which(sim_model==2)),min=0,max=2) 

# simulated summary statistics
sim_S             <- array(NA,nsim)
sim_pi            <- array(NA,nsim)
sim_NH            <- array(NA,nsim)
sim_SFS           <- matrix(NA,nsim,sample_size-1)
sim_TajimasD      <- array(NA,nsim)
sim_FayWuH        <- array(NA,nsim)
sim_FuLiD         <- array(NA,nsim)
for (sim in 1:nsim){
  system("rm temp_ms_file.txt")
  if (sim_model[sim]==1){
    ms(nsam=sample_size, nreps=1, opt=paste("-t",sim_theta[sim]), temp.file ="temp_ms_file.txt")
  }else if (sim_model[sim]==2){
    ms(nsam=sample_size,
       nreps=1,
       opt=paste("-t",sim_theta0[sim],"-eN",sim_tau[sim]/sim_theta0[sim],sim_theta1[sim]/sim_theta0[sim]),
       temp.file ="temp_ms_file.txt")
  }
  msout <- ms.inp.multi(sample_size, 1, ms.output.file="temp_ms_file.txt")
  system("rm temp_ms_file.txt")
  sim_S[sim]             <- msout[[1]]$nss
  sim_pi[sim]            <- thetaPi(msout)
  sim_NH[sim]            <- NH(msout)
  sim_SFS[sim,]          <- SFS(msout)
  sim_TajimasD[sim]      <- tajimaD(msout, thetaW(msout), sim_pi[sim])
  sim_FayWuH[sim]        <- fayWuH(msout)
  sim_FuLiD[sim]         <- fuliD(msout, thetaS1(msout) )
}

# Model (and prior) quality control
plot(sim_S,sim_pi,xlab="S",ylab=expression(pi),log="xy")
points(target_S,target_pi,col="red",cex=3,pch="*")
plot(sim_TajimasD,sim_FayWuH,xlab="Tajima'S D",transfylab="Fay and Wu H")
points(target_TajimasD,target_FayWuH,col="red",cex=3,pch="*")

# Model choice
tolerance <- 0.1
target_SS <- cbind(target_S,target_pi,target_TajimasD,target_FayWuH)
sim_SS    <- cbind(sim_S,sim_pi,sim_TajimasD,sim_FayWuH)
dimnames(target_SS)[[2]]<-dimnames(sim_SS)[[2]]<-c("S","pi","TajimasD","FayWuH")

abc_model_choice <- postpr(target=target_SS,
                           index=sim_model,
                           sumstat=sim_SS,
                           tol=tolerance,
                           method="mnlogistic")
                           
summary(abc_model_choice)                           

sim_params <- cbind(sim_theta0,sim_theta1,sim_tau)
dimnames(sim_params)[[2]] <- c("theta0","theta1","tau")

abc_result <- abc(target=target_SS,
                  param=sim_params[sim_model==2,],
                  sumstat=sim_SS[sim_model==2,],
                  tol=tolerance,
                  transf = c("log","log","none"),
                  method="loclinear")



sim_theta0_adjusted <- as.vector(abc_result$adj.values[,"theta0"])
sim_theta1_adjusted <- as.vector(abc_result$adj.values[,"theta1"])
sim_tau_adjusted <- as.vector(abc_result$adj.values[,"tau"])
sim_weights        <- abc_result$weights

par(mfrow=c(1,3))
# Estimation of posterior probability
hist(log10(sim_theta0),
     breaks=seq(-1,2,0.1),
     col="grey",
     freq=F,
     ylim=c(0,2),
     main="",
     xlab=expression(log(theta[0])),
     ylab="probability density")
wtd.hist(log10(sim_theta0_adjusted), breaks=seq(-1,3,0.1), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
box()
hist(log10(sim_theta1),
     breaks=seq(-1,2,0.1),
     col="grey",
     freq=F,
     ylim=c(0,2),
     main="",
     xlab=expression(log(theta[1])),
     ylab="probability density")
wtd.hist(log10(sim_theta1_adjusted), breaks=seq(-2,2,0.1), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
box()
hist(sim_tau,
     breaks=seq(0,2,0.05),
     col="grey",
     freq=F,
     ylim=c(0,2),
     main="",
     xlab=expression(tau),
     ylab="probability density")
wtd.hist(sim_tau_adjusted, breaks=seq(0,3,0.05), col=rgb(1,0,0,0.7), freq=F, add=T, weight=sim_weights)
box()

# model checking from point estimates
theta0_point_estimate <- summary(abc_result)[3,1]
theta1_point_estimate <- summary(abc_result)[3,2]
tau_point_estimate    <- summary(abc_result)[3,3]

ms(nsam=sample_size, nreps=1000,
   opt=paste("-t",theta0_point_estimate,"-eN",tau_point_estimate/theta0_point_estimate,theta1_point_estimate/theta0_point_estimate),
   temp.file ="data.txt")
pred_seq_data <- ms.inp.multi(sample_size, 1000, ms.output.file="data.txt")
#predicted_S             <- S(pred_seq_data)              # number of segregating sites
#predicted_pi            <- thetaPi(pred_seq_data)              # nucleotide diversity (mean number of pairwise differences)
predicted_NH            <- NH(pred_seq_data)
#predicted_SFS           <- SFS(pred_seq_data)                  # site frequency spectrum
#predicted_TajimasD      <- tajimaD(pred_seq_data, thetaW(pred_seq_data), predicted_pi)
#predicted_FayWuH        <- fayWuH(pred_seq_data)
predicted_FuLiD         <- fuliD(pred_seq_data, thetaS1(pred_seq_data) )

par(mfrow=c(1,2))
hist(predicted_NH, breaks=seq(0,40,1), col="grey", freq=F, xlab="Number of haplotypes", main="")
abline(v=target_NH,col="red",lwd=3)
box()
hist(predicted_FuLiD, col="grey", freq=F, xlab="Fu and Li D", main="")
abline(v=target_FuLiD,col="red",lwd=3)
box()















