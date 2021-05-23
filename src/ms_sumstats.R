# Functions partially based on code by T. Städler, 
# B. Haubold, C. Merino, W..Stephan, and P. Pfaffelhuber
# http://guanine.evolbio.mpg.de/sampling/

# Number of segregating sites
S <- function(a) {
  return(a$segsites)
}

# Number of haplotypes
NH <- function(a) {
  NH<-vector(length=a$nreps)
  for(i in seq_len(a$nreps)) {
    NH[i] <- nrow(unique(a$gametes[[i]]))
  }
  return(NH)
}

# Number of pairwise fdifferences
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

# Theta estimator from number of segregating sites (Watterson 1975, doi:10.1016/0040-5809(75)90020-9)
thetaW<-function(a) {
  a1<-sum(1/(1:(a$nsam-1)))
  return(a$segsites/a1)
}

# Tajima's D
TajimaD<-function(a) {
  a1<-sum(1/(1:(a$nsam-1)))
  a2<-sum(1/((1:(a$nsam-1))^2))
  b1<-(a$nsam+1)/3/(a$nsam-1)
  b2<-2*(a$nsam^2+a$nsam+3)/(9*a$nsam*(a$nsam-1))
  c1<-b1-1/a1
  c2<-b2-(a$nsam+2)/(a1*a$nsam)+a2/(a1^2)
  temp <- c1/a1*a$segsites + c2/(a1^2+a2)*a$segsites*(a$segsites-1)
  TajimaD <- (PI(a)-thetaW(a))/sqrt(temp)
  TajimaD[a$segsites==0] <- 0
  return(TajimaD)
}

# Fu and Li D (Fu and Li 1993; doi:10.1093/genetics/133.3.693)
FuLiD<-function(a) {
  tS1 <- thetaS1(a)
  a1<-sum(1/(1:(a$nsam-1)))
  a2<-sum(1/((1:(a$nsam-1))^2))
  c<-(2*a$nsam*a1-4*(a$nsam-1))/((a$nsam-1)*(a$nsam-2))
  vD<-1+a1^2/(a2+a1^2)*(c-(a$nsam+1)/(a$nsam-1))
  uD<-a1-1-vD
  FLD <- vector(mode='numeric',length=a$nreps)
  for(i in 1:a$nreps) {
    if(a$segsites[i]) {
      FLD[i] <- (a$segsites[i]-a1*tS1[i])/sqrt(uD*a$segsites[i]+vD*a$segsites[i]^2)
    }else{
      FLD[i] <- 0
    }
  }
  FLD
}


# Theta estimator from derived singletons (Fu and Li 1993; doi:10.1093/genetics/133.3.693)
thetaS1<-function(a) {
  thetaS1<-vector(mode='numeric', length=a$nreps)
  for(i in 1:a$nreps) {
    if(a$segsites[i]) {
      for(j in 1:a$segsites[i]) {
        if(a$segsites[i]>1) {
          if(sum(a$gametes[[i]][,j])==1) {
            thetaS1[i]<-thetaS1[i]+1
          }
        } else {
          if(sum(a$gametes[[i]])==1) {
            thetaS1[i]<-thetaS1[i]+1
          }
        }
      }
    }
  }
  thetaS1
}
