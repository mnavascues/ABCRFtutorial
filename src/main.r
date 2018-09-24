##############################################################################
# A Biologist Guide to Approximate Bayesian Computation in Population Genetics
##############################################################################

# Miguel Navascués
# miguel.navascues@inra.fr
# INRA


# install.packages(c("phyclust","abc","weights","abcrf","tree"))

# LOAD PACKAGES, SAVE DATA ON SESSION
#=====================================
## ---- rcode_0

suppressMessages(library("phyclust")) # for a version of Hudson ms (coalescent simulator)
suppressMessages(library("abc"))      # for approximate Bayesian computation
suppressMessages(library("weights"))  # for weighted histograms
suppressMessages(library("abcrf"))    # for ABC with random forests
suppressMessages(library("tree"))     # for CART
suppressMessages(library("highr"))    # for R highligth in pdf handouts (not necessary for tutorial)

source("src/color_blind_palette.R")

recordSessionInfo <- sessionInfo()
R_environment     <- .Platform$GUI
R_version         <- R.version.string
#dir.create("results")
save(recordSessionInfo, R_version, R_environment, file="results/sessionInfo.RData")

# Get information on session
#===========================
# recordSessionInfo$platform
# recordSessionInfo$otherPkgs$abc$Version


#################################
# Introduction to the likelihood 
#################################

# The experiment: A coin is tossed 10 times and the results (heads/tails) are recorded:
#set.seed(29031978)
toss <- sample(c('H','T'),size=10,replace=T,prob=c(0.6,0.4))

## ---- rcode_0end
## ---- rcode_1

toss

total_tosses <- length(toss)
heads_count  <- length(which(toss=="H"))
tails_count  <- total_tosses - heads_count

## ---- rcode_1end

{cat(paste("number of tosses:",total_tosses,"\n"))
 cat(paste("number of heads:",heads_count,"\n"))
 cat(paste("number of tails:",tails_count,"\n"))}

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

## ---- rcode_5
source("src/flip_coin_likelihood.r")
likelihood_profile <- flip.coin.likelihood(n = total_tosses,
                                           k = heads_count,
                                           p = seq(0,1,0.001))
## ---- rcode_5b
maxL  <- max(likelihood_profile$likelihood)
p_hat <- likelihood_profile$p[which(likelihood_profile$likelihood
                                                          ==maxL)]
CI95  <- likelihood_profile$p[which(likelihood_profile$likelihood>
                                             exp(log(maxL)-1.92))]
CI95  <- CI95[c(1,length(CI95))]
## ---- rcode_5end

cat(paste("Maximum likelihood estimate:",p_hat,"\n",
          "with 95% confidence interval:",CI95[1],",",CI95[2],"\n"))

## ---- Exercise1
# source("src/Exercise1.r")
## ---- Exercise1end

plot(x    = likelihood_profile$p,
     y    = likelihood_profile$likelihood,
     xlab = "p",
     ylab = "Likelihood",
     type = "l")
abline(v = p_hat, col = 2, lwd = 2)
abline(h = exp(log(maxL)-1.92), col = 6, lwd = 2)
{abline(v = CI95[1], lty = 2, col = 2, lwd = 2)
 abline(v = CI95[2], lty = 2, col = 2, lwd = 2)}


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
likelihood_approx <- flip.coin.likelihood.approx(n  =total_tosses,
                                                 k  =heads_count,
                                                 p  =p_heads,
                                                 rep=1000)
## ---- rcode_7end
cat(paste("Estimated likelihood of p=",p_heads,"given",heads_count,"heads in",
          total_tosses,"tosses:",likelihood_approx$likelihood,
          "(true value:",likelihood,")\n"))

## ---- Exercise2
# source("src/Exercise2.r")
## ---- Exercise2end

## ---- rcode_8
likelihood_profile_approx <-
                 flip.coin.likelihood.approx(n   = total_tosses,
                                             k   = heads_count,
                                             p   = seq(0,1,0.01),
                                             rep = 1000)
## ---- rcode_8end

plot(x    = likelihood_profile$p,
     y    = likelihood_profile$likelihood,
     xlab = expression(italic(p)), ylab = "Likelihood",
     lwd  = 2, type = "l")
lines(x   = likelihood_profile_approx$p ,
      y   = likelihood_profile_approx$likelihood,
      col = 2, lwd = 2)

## ---- Exercise3
# source("src/Exercise3.r")
## ---- Exercise3end


####################
# BAYESIAN APPROACH            
####################

## ---- rcode_9a
total_draws     <- 10000
sim_parameter_p <- runif(total_draws,0,1) 
prior_p_approx  <- sum(sim_parameter_p==p_heads)/total_draws
## ---- rcode_9a_end
cat(paste("Prior probability of p=",p_heads,":",prior_p_approx,"\n"));
## ---- rcode_9end


# EXACT REJECTION ALGORITHM
#---------------------------
# 1. generate p' from prior
# 2. generate D' from simulation with p'
# 3. accept p' if D'=D; return to 1


## ---- rcode_9b
sim_parameter_p1 <- runif(10000,0,1)
sim_data_p1      <- rbinom(n    = length(sim_parameter_p1),
                           size = total_tosses,
                           prob = sim_parameter_p1)
sim_parameter_kept_p1 <- sim_parameter_p1[which(sim_data_p1==heads_count)]
sim_data_kept_p1      <- sim_data_p1[which(sim_data_p1==heads_count)]
## ---- rcode_9b_end



plot(x    = sim_parameter_p1,
     y    = sim_data_p1,
     xlab = expression(italic(p)*"'"),
     ylab = expression(italic(D)*"'"))
abline(h=heads_count,col=7)
points(sim_parameter_kept_p1,sim_data_kept_p1,col=7)

# Estimation of posterior probability
hist(x      = sim_parameter_p1,
     breaks = seq(0,1,0.02),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,5),
     main   = "",
     xlab   = expression(italic(p)),
     ylab   = "probability density")
hist(x      = sim_parameter_kept_p1,
     breaks = seq(0,1,0.02),
     col    = Vermillion_transparency,
     freq   = FALSE,
     add    = TRUE)
box()
lines(x   = seq(0,1,0.001),
      y   = dbeta(x=seq(0,1,0.001),1+heads_count,1+tails_count),
      col = 7,
      lwd = 3)

## ---- rcode_10
p_hat  <- median(sim_parameter_kept_p1)
p_95CI <- quantile(sim_parameter_kept_p1,probs=c(0.025,0.975))
## ---- rcode_10end
cat(paste(c("Point estimate of p:",round(p_hat,2),"; 95%CI:",
          round(as.vector(p_95CI),2),"\n"),sep=" ",collapse = " "))



## ---- Exercise4
# source("src/Exercise4.r")
## ---- Exercise4end



# OBSERVED DATA (actually simulated)
## ---- rcode_13
sample_size <- 50
true_theta  <- 10
if (!file.exists("data/dataset1.txt")){
  ms(nsam      = sample_size,
     opt       = paste("-t",true_theta),
     temp.file = "data/dataset1.txt")
}
true_theta0  <- 50
true_theta1  <- 0.5
true_tau     <- 0.5
if (!file.exists("data/dataset2.txt")){
  ms(nsam      = sample_size,
     opt       = paste("-t", true_theta0,
                       "-eN", true_tau/true_theta0,
                       true_theta1/true_theta0),
     temp.file = "data/dataset2.txt")
}
## ---- rcode_13end


## ---- rcode_14_A
source("src/ms.r")
seq_data <- ms.inp.multi(sample_size, 1,
                         ms.output.file="data/dataset1.txt")

target1_S        <- S(seq_data)
target1_pi       <- thetaPi(seq_data) 
target1_NH       <- NH(seq_data)
target1_SFS      <- SFS(seq_data)     
target1_TajimasD <- tajimaD(seq_data, thetaW(seq_data), target1_pi)
target1_FayWuH   <- fayWuH(seq_data)
target1_FuLiD    <- fuliD(seq_data, thetaS1(seq_data) )
## ---- rcode_14_Aend

seq_data[[1]]$nss
seq_data[[1]]$mutations
seq_data[[1]]$sample[1:5,1:10]

cat(paste("Data set 2:\n",
          "Sample size:",sample_size,"\n",
          "Number of polymorphic sites:",target1_S,"\n",
          "Number of haplotypes:",target1_NH,"\n"))

colnames(target1_SFS)<-1:(sample_size-1)
barplot(height = target1_SFS/target1_S,
        main   = "Unfolded Site Frequency Spectrum",
        xlab   = "derived allele count in sample",
        ylab   = "Proportion of sites")
box()

## ---- Exercise5
source("src/Exercise5.r")
## ---- Exercise5end



######################################
# APPROXIMATION 2: SUMMARY STATISTICS
######################################

# EXACT REJECTION ALGORITHM ON SUMMARY STATISTICS
#-------------------------------------------------
# 1. generate theta' from prior
# 2. generate D' from simulation with theta'
# 3. calculate SS' from D' 
# 4. accept theta' if SS'=SS; return to 1


## ---- rcode_15
nsim <- 10000
sim_theta <- 10^runif(nsim,min=-1,max=2) # log uniform prior

if (file.exists("data/abc_sims1.txt")){
  file_removed <- file.remove("data/abc_sims1.txt")
}
ms(nsam       = sample_size, 
   opt        = "-t tbs",
   tbs.matrix = cbind(sim_theta),
   temp.file  = "data/abc_sims1.txt")

msout <- ms.inp.multi(sample_size, 
                      ndraws=nsim, 
                      ms.output.file="data/abc_sims1.txt")
sim1_S <- S(msout)

## ---- rcode_15_b
sim1_pi       <- thetaPi(msout)
sim1_NH       <- NH(msout)
sim1_SFS      <- SFS(msout)
sim1_TajimasD <- tajimaD(msout, thetaW(msout), sim1_pi)
sim1_FayWuH   <- fayWuH(msout)
sim1_FuLiD    <- fuliD(msout, thetaS1(msout) )

sim1_pi       <- sim1_pi[which(sim1_S!=0)]
sim1_NH       <- sim1_NH[which(sim1_S!=0)]
sim1_SFS      <- sim1_SFS[which(sim1_S!=0)]
sim1_TajimasD <- sim1_TajimasD[which(sim1_S!=0)]
sim1_FayWuH   <- sim1_FayWuH[which(sim1_S!=0)]
sim1_FuLiD    <- sim1_FuLiD[which(sim1_S!=0)]
sim_theta    <- as.matrix(sim_theta[which(sim1_S!=0)])
sim1_S        <- sim1_S[which(sim1_S!=0)]
colnames(sim_theta) <- "theta"
## ---- rcode_15end




## ---- rcode_16
target_SS <- target1_S
sim_SS    <- sim1_S
sim_theta_kept <- sim_theta[which(sim_SS==target_SS)]
sim_SS_kept    <- sim_SS[which(sim_SS==target_SS)]
## ---- rcode_16end

plot( log10(sim_theta), sim_SS,
      xlab = expression(log[10](theta*"'")),
      ylab = "SS'",
      ylim = c(0,100))

abline(h=target_SS,col=7)
points(log10(sim_theta_kept),sim_SS_kept,col=7)

cat(paste("The number of simulations kept is",
          length(sim_SS_kept), "out of",  length(sim_SS)))

# Estimation of posterior probability
hist(x      = log10(sim_theta),
     breaks = seq(-1,2,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,6),
     main   = "",
     xlab   = expression(log[10](theta)),
     ylab   = "probability density")
hist(x      = log10(sim_theta_kept),
     breaks = seq(-1,2,0.05), 
     col    = Vermillion_transparency,
     freq   = FALSE,
     add    = TRUE); box()


## ---- Exercise6
source("src/Exercise6.r")
## ---- Exercise6end


########################################
# APPROXIMATION 3: TOLERANCE (DISTANCE)
########################################

# REJECTION ALGORITHM (ABC)
#-------------------------------------------------
# 1. generate theta' from prior
# 2. generate D' from simulation with theta'
# 3. calculate SS' from D' 
# 4. accept theta' if distance(SS',SS) < delta; return to 1

## ---- rcode_17
tolerance  <- 0.1
target_SS  <- target1_pi
sim_SS     <- as.matrix(sim1_pi); colnames(sim_SS) <- "pi"
abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = sim_SS,
                  tol     = tolerance,
                  method  = "rejection")
sim_theta_kept_tol1 <- as.vector(abc_result$unadj.values)
sim_SS_kept_tol1    <- abc_result$ss
## ---- rcode_17end


plot(x    = log10(sim_theta),
     y    = sim_SS,
     xlab = expression(log[10](theta*"'")),
     ylab = "SS'",
     ylim = c(0,50))
abline(h = target_SS, col = 7)
{abline(h = max(abc_result$ss), col = 7, lty = 2)
 abline(h = min(abc_result$ss), col = 7, lty = 2)}
points(x   = log10(sim_theta_kept_tol1),
       y   = sim_SS_kept_tol1,
       col = 7)

cat(paste("The number of simulations kept is",
          length(sim_SS_kept_tol1), "out of",  dim(sim_SS)[1],
          ";\n determined by tolerance of",tolerance))


# Estimation of posterior probability
hist(x      = log10(sim_theta),
     breaks = seq(-1,2,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,4),
     main   = "",
     xlab   = expression(log[10](theta)),
     ylab   = "probability density")
hist(x      = log10(sim_theta_kept_tol1),
     breaks = seq(-1,2,0.05),
     col    = Vermillion_transparency, 
     freq   = FALSE,
     add    = TRUE); box()


## ---- Exercise7
source("src/Exercise7.r")
## ---- Exercise7end



########################
# REGRESSION ADJUSTMENT
########################




## ---- rcode_18
tolerance  <- 0.5
target_SS  <- target1_pi
sim_SS     <- as.matrix(sim1_pi); colnames(sim_SS) <- "pi"
abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = as.matrix(sim_SS),
                  tol     = tolerance,
                  transf  = "log",
                  method  = "loclinear")

sim_theta_kept     <- as.vector(abc_result$unadj.values)
sim_theta_adjusted <- as.vector(abc_result$adj.values)
sim_SS_kept        <- abc_result$ss
local_regression   <- lm(log10(sim_theta_kept)~sim_SS_kept,
                         weights=abc_result$weights)
## ---- rcode_18end

plot(sim_SS_kept,
     abc_result$weights,
     xlab="SS'",
     ylab="weight")
abline(v=target_SS,col=7)

plot(x    = log10(sim_theta),
     y    = sim_SS,
     xlab = expression(log[10](theta*"'")),
     ylab = "SS'",
     ylim = c(0,30))
abline(h = target_SS, col = 6)
{abline(h = max(abc_result$ss), col = 6, lty = 2)
 abline(h = min(abc_result$ss), col = 6, lty = 2)}
points(log10(sim_theta_kept), sim_SS_kept, col = 6)
abline(a   = -local_regression$coefficients[1]/local_regression$coefficients[2],
       b   = 1/local_regression$coefficients[2],
       col = 5,
       lwd = 3)
sim_kept <- 1:10
points(x   = log10(sim_theta_kept)[sim_kept],
       y   = sim_SS_kept[sim_kept],
       col = 5)
arrows(x0     = log10(sim_theta_kept)[sim_kept],
       y0     = sim_SS_kept[sim_kept],
       x1     = log10(sim_theta_adjusted)[sim_kept],
       y1     = target_SS,
       col    = 5,
       length = 0.1)
points(x   = log10(sim_theta_adjusted),
       y   = array(target_SS,length(sim_theta_adjusted)),
       col = 5)


# Estimation of posterior probability
hist(x      = log10(sim_theta),
     breaks = seq(-1,2,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,4),
     main   = "",
     xlab   = expression(log(theta)),
     ylab   = "probability density")
hist(x      = log10(sim_theta_kept),
     breaks = seq(-1,2,0.05),
     col    = Blue_transparency, freq=FALSE, add=TRUE )
wtd.hist(x      = log10(sim_theta_adjusted),
         breaks = seq(-1.5,2,0.05),
         col    = Yellow_transparency, freq=FALSE, add=TRUE,
         weight = abc_result$weights)
hist(x      = log10(sim_theta_kept_tol1),
     breaks = seq(-1,2,0.05),
     col    = Vermillion_transparency, freq=FALSE, add=TRUE)
box()



















####################################
# INFORMATION ON SUMMARY STATISTICS
####################################


## ---- rcode_19
target_SS <- target1_TajimasD
sim_SS    <- cbind(sim1_TajimasD); colnames(sim_SS) <- "Tajima's D"
abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = as.matrix(sim_SS),
                  tol     = tolerance,
                  transf  = "log",
                  method  = "loclinear")
sim_theta_adjusted_D <- as.vector(abc_result$adj.values)
sim_weights_D        <- abc_result$weights
## ---- rcode_19end

hist(x      = log10(sim_theta),
     breaks = seq(-1,2,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,4),
     main   = "",
     xlab   = expression(log[10](theta)),
     ylab   = "probability density")
wtd.hist(x      = log10(sim_theta_adjusted_D),
         breaks = seq(-2,3,0.05),
         col    = BluishGreen_transparency,
         freq   = FALSE, add = TRUE,
         weight = sim_weights_D); box()

## ---- Exercise8
source("src/Exercise8.r")
## ---- Exercise8end



## ---- rcode_20
cor(sim1_TajimasD,as.vector(sim_theta))
cor(sim1_pi,as.vector(sim_theta))

## ---- rcode_21
target_SS <- cbind(target1_pi, 
                   target1_TajimasD)
sim_SS    <- cbind(sim1_pi, 
                   sim1_TajimasD)
colnames(target_SS) <- colnames(sim_SS) <- c("pi","TD")
abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = sim_SS,
                  tol     = tolerance,
                  transf  = "log",
                  method  = "loclinear")
sim_theta_adjusted_piD <- as.vector(abc_result$adj.values)
sim_weights_piD        <- abc_result$weights
## ---- rcode_21end

hist(x      = log10(sim_theta),
     breaks = seq(-1,2,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,4),
     main   = "",
     xlab   = expression(log[10](theta)),
     ylab   = "probability density")
wtd.hist(x      = log10(sim_theta_adjusted_piD),
         breaks = seq(-2,3,0.05),
         col    = BluishGreen_transparency,
         freq   = FALSE, add = TRUE,
         weight = sim_weights_piD); box()
wtd.hist(x      = log10(sim_theta_adjusted_pi),
         breaks = seq(-2,2,0.05), 
         col    = Yellow_transparency, 
         freq   = FALSE, add = TRUE,
         weight = sim_weights_pi)

## ---- rcode_22
target1_noise <- 76
sim_noise     <- 10^runif(length(sim1_pi),-2,3)
target_SS <- cbind(target1_noise, 
                   target1_pi)
sim_SS     <- cbind(sim_noise, 
                    sim1_pi)
colnames(target_SS) <- colnames(sim_SS) <- c("noise","pi")
abc_result <- abc(target  = target_SS,
                  param   = sim_theta,
                  sumstat = as.matrix(sim_SS),
                  tol     = tolerance,
                  transf  = "log",
                  method  = "loclinear")
sim_theta_adjusted_pi_noise <- as.vector(abc_result$adj.values)
sim_weights_pi_noise        <- abc_result$weights
## ---- rcode_22end

hist(x      = log10(sim_theta),
     breaks = seq(-1,2,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,4),
     main   = "",
     xlab   = expression(log[10](theta)),
     ylab   = "probability density")
wtd.hist(x      = log10(sim_theta_adjusted_pi_noise), 
         breaks = seq(-2,2,0.05), 
         col    = ReddishPurple_transparency, 
         freq   = FALSE, add = TRUE, 
         weight = sim_weights_pi_noise)
wtd.hist(log10(sim_theta_adjusted_pi), breaks=seq(-1,2,0.05),
         col=Yellow_transparency, freq=F, add=T, weight=sim_weights_pi)
wtd.hist(log10(sim_theta_adjusted_piD), breaks=seq(-2,2,0.05), 
         col=BluishGreen_transparency, freq=F, add=T, weight=sim_weights_piD)
box()

















######################
# BEWARE OF THE PRIOR
######################

# prior on mutation rate (3e-8 - 1e-6) and population size (1-30 million)

## ---- rcode_23
sim_mu <- 10^runif(nsim,min=-7.5,max=-6)
sim_Ne <- 10^runif(nsim,min=6,max=7.5)
sim_theta_composite <- 4*sim_Ne*sim_mu
## ---- rcode_23end

hist(x      = log10(sim_theta_composite),
     breaks = seq(-2,3,0.05),
     col    = "grey",
     freq   = FALSE,
     ylim   = c(0,1),
     xlim   = c(-1,2),
     main   = "",
     xlab   = expression(log[10](theta)),
     ylab   = "probability density"); box()

## ---- Exercise9
source("src/Exercise9.r")
## ---- Exercise9end

##################
# QUALITY CONTROL
##################

# 1. Assessment of the method per se (validation through simulted data)

# cross validation

## ---- rcode_24
cross_validation_result <- cv4abc(param   = sim_theta,
                                  sumstat = as.matrix(sim_SS),
                                  abc.out = abc_result1,
                                  nval    = 100,
                                  tols    = 0.1)
## ---- rcode_24end

plot(x    = log10(cross_validation_result$true$theta),
     y    = log10(cross_validation_result$estim$tol0.1),
     xlab = expression(log(theta)),
     ylab = expression(log(hat(theta))),
     xlim = c(-1,2), ylim = c(-1,2)); abline(a = 0, b = 1)


# 2. Assessment of the applicability of the method to the observed data

# Are the simulations producing pseudo-data similar to target data?
# Summary statistics used in the analysis
par(mfrow=c(3,1),mar=c(4.2,4.2,1,1))
plot(sim1_S,sim1_pi,xlab="S",ylab=expression(pi),log="xy")
points(target1_S,target1_pi,col=7,cex=4,pch="*")
points(target2_S,target2_pi,col=6,cex=4,pch="*")

plot(sim1_NH,sim1_pi,xlab="Number of Haplotypes",ylab=expression(pi),log="xy")
points(target1_NH,target1_pi,col=7,cex=4,pch="*")
points(target2_NH,target2_pi,col=6,cex=4,pch="*")

plot(sim1_FuLiD,sim1_pi,xlab="Fu and Li's D",ylab=expression(pi),log="y",xlim=c(-6,2))
points(target1_FuLiD,target1_pi,col=7,cex=4,pch="*")
points(target2_FuLiD,target2_pi,col=6,cex=4,pch="*")


# Simulations from the posterior probability of the parameter
# Sample parameter theta=4*N*mu from posterior 

## ---- rcode_25
posterior_theta <- sample(sim_theta_adjusted_1,
                          size    = 1000,
                          replace = TRUE,
                          prob    = sim_weights_1)
if (file.exists("data/abc_sims_posterior.txt")){
  file_removed <- file.remove("data/abc_sims_posterior.txt")
}
ms(nsam       = sample_size, 
   opt        = "-t tbs",
   tbs.matrix = cbind(posterior_theta),
   temp.file  = "data/abc_sims_posterior.txt")
msout <- ms.inp.multi(sample_size, 1000, 
                      ms.output.file="data/abc_sims_posterior.txt")
posterior_S        <- S(msout)
posterior_pi       <- thetaPi(msout)
posterior_NH       <- NH(msout)
posterior_TajimasD <- tajimaD(msout, thetaW(msout), posterior_pi)
posterior_FayWuH   <- fayWuH(msout)
posterior_FuLiD    <- fuliD(msout, thetaS1(msout) )
## ---- rcode_25end


par(mfrow=c(2,1))
hist(x      = posterior_NH, 
     breaks = seq(0,40,1),
     col    = "grey",
     freq   = FALSE,
     xlab   = "Number of haplotypes", main = "")
abline(v = target1_NH, col = 7, lwd = 2); box()
hist(x      = posterior_FuLiD,
     breaks = seq(-5,4,0.2),
     col    = "grey",
     freq   = FALSE, 
     xlab   = "Fu and Li's D", main = "")
abline(v = target1_FuLiD, col = 7, lwd = 2); box()


## ---- Exercise10
source("src/Exercise10.r")
## ---- Exercise10end
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))






















######################
# MODEL CHOICE IN ABC
######################

# Consider two alternative models: (1) constant size and (2) population expansion



## ---- rcode_26
nsim <- 10000
sim_theta0 <- 10^runif(nsim,min=-1,max=2) 
sim_theta1 <- 10^runif(nsim,min=-1,max=2) 
sim_tau    <- runif(nsim,min=0,max=2) 

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
sim2_S        <- S(msout)
sim2_pi       <- thetaPi(msout)
sim2_NH       <- NH(msout)
sim2_TajimasD <- tajimaD(msout, thetaW(msout), sim2_pi)
sim2_FayWuH   <- fayWuH(msout)
sim2_FuLiD    <- fuliD(msout, thetaS1(msout) )
## ---- rcode_26_b

sim2_pi       <- sim2_pi[which(sim2_S!=0)]
sim2_NH       <- sim2_NH[which(sim2_S!=0)]
sim2_TajimasD <- sim2_TajimasD[which(sim2_S!=0)]
sim2_FayWuH   <- sim2_FayWuH[which(sim2_S!=0)]
sim2_FuLiD    <- sim2_FuLiD[which(sim2_S!=0)]
ref_table_params <- ref_table_params[which(sim2_S!=0),]
sim2_S        <- sim2_S[which(sim2_S!=0)]

## ---- rcode_26_end

plot(sim2_S,sim2_pi,xlab="S",ylab=expression(pi),log="xy")
points(sim1_S,sim1_pi,col="grey")
points(target1_S,target1_pi,col=7,cex=3,pch="*")
points(target2_S,target2_pi,col=6,cex=3,pch="*")
plot(sim2_FuLiD,sim2_pi,xlab="Fu and Li's D",ylab=expression(pi),log="y",xlim=c(-6,2))
points(sim1_FuLiD,sim1_pi,col="grey")
points(target1_FuLiD,target1_pi,col=7,cex=3,pch="*")
points(target2_FuLiD,target2_pi,col=6,cex=3,pch="*")



# Model choice


sim_model    <- c(array("C",9000),array("V",9000))
sim_S        <- c(sim1_S[1:9000],sim2_S[1:9000])
sim_pi       <- c(sim1_pi[1:9000],sim2_pi[1:9000])
sim_NH       <- c(sim1_NH[1:9000],sim2_NH[1:9000])
sim_TajimasD <- c(sim1_TajimasD[1:9000],sim2_TajimasD[1:9000])
sim_FayWuH   <- c(sim1_FayWuH[1:9000],sim2_FayWuH[1:9000])
sim_FuLiD    <- c(sim1_FuLiD[1:9000],sim2_FuLiD[1:9000])




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

theta  <- ref_table_params[,"theta1"]
TajD   <- sim2_TajimasD
pi     <- sim2_pi
ref_table <- data.frame(theta,TajD,pi)

regression_tree <- tree(theta ~ TajD + pi,
                        data=ref_table)

plot(regression_tree)
text(regression_tree,cex=0.75)

plot(ref_table$TajD,
     ref_table$pi,
     xlab="Tajima's D",
     ylab=expression(pi),
     col=cbPalette1[round(log10(theta)+3)],
     #col=grey(1-theta/max(theta)),
     pch=20)
partition.tree(regression_tree,
               ordvars=c("TajD","pi"),
               add=T,cex=1)


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

theta     <- ref_table_params[,"theta0"]
pi        <- sim2_pi

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

err.abcrf(model_RF,
          training=ref_table,
          paral=T)

###################################################
# HERE BE DRAGONS



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





log10theta <- log10(sim_theta)
sim_SS <- cbind(sim1_S,sim1_pi,sim1_NH,
                sim1_TajimasD,sim1_FayWuH,sim1_FuLiD)
ref_table <- data.frame(log10theta,sim_SS)
colnames(ref_table) <- c("log10theta","S","pi","NH","TD","FWH","FLD")
RFmodel_theta <- regAbcrf(formula = log10theta~.,
                           data    = ref_table,
                           ntree   = 1000,
                           paral   = T)


plot(RFmodel_theta)
err.regAbcrf(object   = RFmodel_theta,
             training = ref_table,
             paral    = T)
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

