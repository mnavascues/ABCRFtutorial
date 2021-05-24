# function to calculate likelihood for the coin toss experiment
flip_coin_likelihood <- function(p,h,t){
  likelihood <- choose(h+t,h) * p^h * (1-p)^t
  return (list(p=p,l=likelihood))
}

flip_coin_likelihood_approx <- function(p,h,t,num_of_sims){
  likelihood <- array(dim=length(p))
  for (i in seq_along(p)){
    likelihood[i] <- sum(rbinom(n=num_of_sims,size=h+t,prob=p[i])==h)/num_of_sims  
  }
  return (list(p=p,l=likelihood))
}

flip_coin_posterior_approx <- function(h,t,a,b,num_of_sims){
  p_prime <- rbeta(num_of_sims,a,b)
  d_prime <- rbinom(num_of_sims,h+t,p_prime)
  p_prime_accept <- p_prime[which(d_prime==h)]
  d_prime_accept <- d_prime[which(d_prime==h)]
  return(list(p_prime=p_prime,
              d_prime=d_prime,
              p_prime_accept=p_prime_accept,
              d_prime_accept=d_prime_accept))
}

