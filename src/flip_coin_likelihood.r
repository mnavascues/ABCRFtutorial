flip_coin_likelihood <- function(n,k,p){
  q <- 1-p
  l <- n - k
  likelihood <- choose(n,k) * p^k * q^l
  return (list(p=p,likelihood=likelihood))
}

flip_coin_likelihood_approx <- function(n,k,p,rep,trace=F){
  likelihood_approx <- array(NA,length(p))
  for (i in 1:length(p)){
    counts <- array(NA,rep)
    for (j in 1:rep){
      toss_simulation<-rbinom(n,1,p[i])
      counts[j] <- length(which(toss_simulation==1))
      if (j<11 && trace){
        toss_simulation[which(toss_simulation==1)] <- "H"  
        toss_simulation[which(toss_simulation==0)] <- "T"  
        print(toss_simulation)
      }
    }
    likelihood_approx[i] <- length(which(counts==k))/rep
  }
  return (list(p=p,likelihood=likelihood_approx))
}
