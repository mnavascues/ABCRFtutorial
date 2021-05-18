# function to calculate likelihood for the coin toss experiment
flip_coin_likelihood <- function(p,h,t){
  likelihood <- choose(h+t,h) * p^h * (1-p)^t
  return (list(p=p,l=likelihood))
}