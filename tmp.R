library("phyclust", quietly=T)
source("readms.output.R", local =T)

n<-10
l<-10000
N<-10000
u<-1e-8
num_of_sim<-3
theta <- 2*N*u*l
msout <- ms(nsam=n,nreps=num_of_sim,opt=paste("-t",theta))
sim_data <- read.ms.output(txt=msout)


S<-function(a) {
  return(a$segsites)
}


NH<-function(a) {
  NH<-vector(length=a$nreps)
  for(i in seq_len(a$nreps)) {
    NH[i] <- nrow(unique(a$gametes[[i]]))
  }
  return(NH)
}

PI<-function(a) {
  PI<-vector(mode='numeric', length=a$nreps)
  for(i in seq_len(a$nreps)) {
    if(a$segsites[i]>0){
      for(j in seq_len(a$segsites[i])) {
         if(a$segsites[i]>1) {
          loc<-sum(a$gametes[[i]][,j])
        } else {
          loc<-sum(a$gametes[[i]])
        }
        PI[i] <- PI[i] + 2*loc*(a$nsam-loc)
      }
    }
  }
  return(PI/(a$nsam*(a$nsam-1)))
}

thetaW<-function(a) {
  a1<-sum(1/(1:(a$nsam-1)))
  return(a$segsites/a1)
}


TajimaD<-function(a) {
  a1<-sum(1/(1:(a$nsam-1)))
  a2<-sum(1/((1:(a$nsam-1))^2))
  b1<-(a$nsam+1)/3/(a$nsam-1)
  b2<-2*(a$nsam^2+a$nsam+3)/(9*a$nsam*(a$nsam-1))
  c1<-b1-1/a1
  c2<-b2-(a$nsam+2)/(a1*a$nsam)+a2/(a1^2)
  temp <- c1/a1*a$segsites + c2/(a1^2+a2)*a$segsites*(a$segsites-1)
  TajimaD <- (PI(a)-thetaW(a))/sqrt(temp) 
  return(TajimaD)
}








thetaW(sim_data)
PI(sim_data)
S(sim_data)
NH(sim_data)
TajimaD(sim_data)
