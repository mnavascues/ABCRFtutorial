# Sampling Version 0.5: Software for Simulating Haplotypes in Non-Equilibrium Subdivided Populations
# by Thomas Städler, Bernhard Haubold, Carlos Merino, Wolfgang Stephan, and Peter Pfaffelhuber
# http://guanine.evolbio.mpg.de/sampling/


# A FUNKTION TO IMPORT OUTPUT OF MS BASED ON MSSAMPLES.R BY R. Hudson
## !!MULTIPLE DRAWS!!
ms.inp.multi <- function(nsam, ndraws, ms.output.file="ms.out") {
  txt <- scan(file=ms.output.file, what=character(0), sep="\n", quiet=TRUE)

  h <- numeric()
  result <- list()

  ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
  marker <- grep("segsites", txt)
  stopifnot(length(marker) == ndraws)

  ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
  segsites <- onlyints(txt,marker)    

  lnr1 <- 1
  lnr2 <- 0
  positions <- onlydoubles(txt,(marker+1))
  mutations <- list()
  m <- list()
  for(draw in seq(along=marker)) {
    #if(!(draw %% 10)) cat(draw, " \n ")
    if(segsites[draw] > 0) {
      haplotypes <- txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
      haplotypes <- strsplit(haplotypes, split="")
      h <- sapply(haplotypes, function(el) c(as.integer(el)))

      ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX 
      ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
      if(segsites[draw] == 1) {
        h <- as.matrix(h)
      } else {
        h<-t(h)
      }
      lnr2 <- lnr2 + segsites[draw]
      mutations[[draw]] <- positions[lnr1:lnr2]
      if(segsites[draw]) m <- mutations[[draw]]
      lnr1 <- (lnr2+1)
    } else {
      h <- matrix(nrow=nsam, ncol=0)
      m <- list()
    }
    result[[draw]] <- list(nss=segsites[draw], mutations=m,sample=h)
    stopifnot(all(dim(h) == c(nsam, segsites[draw])) ) 
  }
  cat("\n")
  result 
}

# A FUNCTION TO GET ONLY THE NR OF SEGREGATING SITES IN A TEXT VECTOR
onlyints <- function(txt,marker) {
  a <-strsplit(txt[marker], split=" ")
  b <- unlist(a)
  c <- gsub("segsites:", "", b)
  d <- as.integer(c)
  e <- na.omit(d)
  as.vector(e)
}

# A FUNCTION TO GET ONLY THE POSITIONS OF SNPS IN A TEXT VECTOR
onlydoubles <- function(txt,marker) {
  a <-strsplit(txt[marker], split=" ")
  b <- unlist(a)
  c <- gsub("positions:", "", b)
  d <- as.double(c)
  e <- na.omit(d)
  as.vector(e)
}

# produce migration matrix either according to island model or stepping stone

makeMigMatrix<-function(demes, migrationRate, migrationPattern) {
  migMatrix<-matrix(0, nrow=demes, ncol=demes)
  if(migrationPattern==1) {
    migMatrix<-migrationRate/(demes-1) + 0*(1:(demes^2))
  } else if (migrationPattern==2) {
    migMatrix<-matrix(0, nrow=demes, ncol=demes)
    l<-sqrt(demes)
    if(l == as.integer(l)) {
      for(i in 0:(demes-1)) {
        for(j in 0:(demes-1)) {
          xi<-as.integer(i/l)
          yi<-i%%l
          xj<-as.integer(j/l)
          yj<-j%%l
          if(abs(xi-xj)+abs(yi-yj)==1 || 
            (xi==xj) && abs(yj-yi)==(l-1) || 
            (yi==yj) && abs(xj-xi)==(l-1) ) {
            migMatrix[i+1,j+1]<-0.25*migrationRate
          }
        }
      }
    }
  }
  as.vector(migMatrix)
}

makeStats<-function(a, nsam) {
  ndraws<-length(a)
  S<-S(a)
  thetaW<-thetaW(a)
  thetaPi<-thetaPi(a)
  thetaS1<-thetaS1(a)
  thetaSg1<-thetaSg1(a)
  thetaFayWu<-thetaFayWu(a)
  tajimaD<-tajimaD(a, thetaW, thetaPi)
  fuliD<-fuliD(a, thetaS1)
  fayWuH<-fayWuH(a)
  MFDM<-MFDM(a)
  fst<-fst(a, nsam)
  list(nsam=nsam, S=S, thetaW=thetaW, thetaPi=thetaPi, thetaS1=thetaS1, thetaSg1=thetaSg1, thetaFayWu=thetaFayWu, tajimaD=tajimaD, fuliD=fuliD, fayWuH=fayWuH, MFDM=MFDM, fst=fst)
}

# M. Navascués: number of polymorphic sites
S<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  S<-vector(length=ndraws)
  for(i in 1:ndraws) {
    S[i]<-a[[i]]$nss
  }
  S
}

# M. Navascués: number of haplotypes
NH<-function(a) {
  ndraws<-length(a)
  NH<-vector(length=ndraws)
  for(i in 1:ndraws) {
    NH[i] <- nrow(unique(a[[i]]$sample))
  }
  NH
}



# M. Navascués: theta estimator from S (Watterson 1975)
thetaW<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  a1<-sum(1/(1:(nsam-1)))
  thetaW<-vector(length=ndraws)
  for(i in 1:ndraws) {
    thetaW[i]<-a[[i]]$nss/a1
  }
  thetaW
}

# M. Navascués: Theta estimator from pi (Tajima 1983)
thetaPi<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  thetaPi<-vector(mode='numeric', length=ndraws)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      for(j in 1:(a[[i]]$nss)) {
        if(a[[i]]$nss>1) {
          loc<-sum(a[[i]]$sample[,j])
        } else {
          loc<-sum(a[[i]]$sample)
        }
        thetaPi[i] <- thetaPi[i] + 2*loc*(nsam-loc)
      }
    }
  }
  thetaPi/(nsam*(nsam-1))
}

# M. Navascués: Theta estimator from derived singletons (Fu and Li 1993)
thetaS1<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  thetaS1<-vector(mode='numeric', length=ndraws)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      for(j in 1:a[[i]]$nss) {
        if(a[[i]]$nss>1) {
          if(sum(a[[i]]$sample[,j])==1) {
            thetaS1[i]<-thetaS1[i]+1
          }
        } else {
          if(sum(a[[i]]$sample)==1) {
            thetaS1[i]<-thetaS1[i]+1
          }
        }
      }
    }
  }
  thetaS1
}



# M. Navascués: function to calculate the maximum frequency of derived mutations (Li 2011)
MFDM<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  MFDM<-vector(mode='numeric', length=ndraws)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      MFDM[i]<-0.0
      for(j in 1:a[[i]]$nss) {
        if(a[[i]]$nss>1) {
          if(sum(a[[i]]$sample[,j])/nsam>MFDM[i]) {
            MFDM[i]<-sum(a[[i]]$sample[,j])/nsam
          }
        } else {
          MFDM[i]<-sum(a[[i]]$sample)/nsam
        }
      }
    }
  }
  MFDM
}

# M. Navascués: function to calculate SFS
SFS<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  SFS<-matrix(0,ndraws,nsam-1)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      for(j in 1:a[[i]]$nss) {
        freq <- sum(a[[i]]$sample[,j])  
        SFS[i,freq] <- SFS[i,freq]+1
      }
    }
  }
  SFS
}

# M. Navascués: function to mean He
meanHe<-function(a,length_DNA) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  meanHe<-vector(mode='numeric', length=ndraws)
  for(i in 1:ndraws) {
    meanHe[i] <- 0
    if(a[[i]]$nss) {
      
      for(j in 1:a[[i]]$nss) {
        meanHe[i] <- meanHe[i] + (a[[i]]$nss/(a[[i]]$nss-1)) * (1 - (sum(a[[i]]$sample[,j])/nsam)^2)  
        #meanHe[i] <- meanHe[i] + (1 - (sum(a[[i]]$sample[,j])/nsam)^2)  
      }
      meanHe[i] <- meanHe[i]/length_DNA
    }
  }
  meanHe
}









# M. Navascués : Theta estimator from number of non-singleton segregating sites
thetaSg1<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  a1<-sum(1/(2:(nsam-1)))
  thetaSg1<-vector(mode='numeric', length=ndraws)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      for(j in 1:a[[i]]$nss) {
        if(a[[i]]$nss>1) {
          if(sum(a[[i]]$sample[,j])>1) {
            thetaSg1[i]<-thetaSg1[i]+1
          }
        } else {
          if(sum(a[[i]]$sample)>1) {
            thetaSg1[i]<-thetaSg1[i]+1
          }
        }
      }
    }
  }
  thetaSg1/a1
}

# M. Navascués: Fay and Wu (2000) theta estimator
thetaFayWu<-function(a) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  thetaFayWu<-vector(mode='numeric', length=ndraws)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      for(j in 1:a[[i]]$nss) {
        if(a[[i]]$nss>1) {
          freq<-sum(a[[i]]$sample[,j])
          thetaFayWu[i]<-thetaFayWu[i]+freq*freq
        } else {
          freq<-sum(a[[i]]$sample)
          thetaFayWu[i]<-thetaFayWu[i]+freq*freq
        }
      }
    thetaFayWu[i]<-thetaFayWu[i]*2/(nsam*(nsam-1))
    }
  }
  thetaFayWu
}

# Tajima's D
tajimaD<-function(a, thetaW, thetaPi) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  a1<-sum(1/(1:(nsam-1)))
  a2<-sum(1/((1:(nsam-1))^2))
  b1<-(nsam+1)/3/(nsam-1)
  b2<-2*(nsam^2+nsam+3)/(9*nsam*(nsam-1))
  c1<-b1-1/a1
  c2<-b2-(nsam+2)/(a1*nsam)+a2/(a1^2)
  tajimaD<-vector(mode='numeric',length=ndraws)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      locS<-a[[i]]$nss
      loc<-c1/a1*locS + c2/(a1^2+a2)*locS*(locS-1)
      tajimaD[i]<-(thetaPi[i]-thetaW[i])/sqrt(loc)
    }
  }
  tajimaD
}

# Fu and Li (1993) D
fuliD<-function(a, thetaS1) {
  ndraws<-length(a)
  nsam<-dim(a[[1]]$sample)[1]
  a1<-sum(1/(1:(nsam-1)))
  a2<-sum(1/((1:(nsam-1))^2))
  c<-(2*nsam*a1-4*(nsam-1))/((nsam-1)*(nsam-2))
  vD<-1+a1^2/(a2+a1^2)*(c-(nsam+1)/(nsam-1))
  uD<-a1-1-vD
  fuliD<-vector(mode='numeric',length=ndraws)
  for(i in 1:ndraws) {
    if(a[[i]]$nss) {
      locS<-a[[i]]$nss
      loc<-uD*locS+vD*locS^2
      fuliD[i]<-(locS-a1*thetaS1[i])/sqrt(loc)
    }
  }
  fuliD
}

fayWuH<-function(a) {
  thetaFayWu(a)-thetaPi(a)
}

fst<-function(a, nsam) {
  ndraws<-length(a)
  nsam<-nsam[nsam>0]

  if(length(nsam[nsam==1]) || length(nsam)==1) {
    fst<-0
    cat("Computation of FST only possible if at least two sequences are present in each deme")
  } else {
    # find out break points between demes
    nsamBrk<-nsam*0
    nsamBrk[1]<-1
    for(i in 2:length(nsam))
      nsamBrk[i]<-nsamBrk[i-1]+nsam[i-1]
    nsamBrk<-c(nsamBrk, sum(nsam)+1)

    demes<-length(nsam)
    Hw<-vector(mode='numeric', length=ndraws)
    Hb<-vector(mode='numeric', length=ndraws)
    for(i in 1:ndraws) {
      if(a[[i]]$nss) {
        for(j in 1:a[[i]]$nss) {
          for(k in 1:demes) {
            lock<-sum(a[[i]]$sample[nsamBrk[k]:(nsamBrk[k+1]-1),j])
            Hw[i]<-Hw[i]+2*lock*(nsam[k]-lock)
            for(l in 1:demes) {
              if(k!=l)  
              {
                locl<-sum(a[[i]]$sample[nsamBrk[l]:(nsamBrk[l+1]-1),j])
                Hb[i]<-Hb[i] + lock*(nsam[l]-locl) + (nsam[k]-lock) * locl
              }
            }
          }
        }
      }
    }
    within<-0
    between<-0
    for(k in 1:demes) {
      within<-within + nsam[k]*(nsam[k]-1)
      for(l in 1:demes) {
        if(k!=l) {
          between<-between + nsam[k]*nsam[l]
        }
      }
    }
    Hw<-Hw/within
    Hb<-Hb/between
    fst<-1-Hw/Hb
  }
  fst
}


# make ms command for a model where all demes desccend from a single deme a time tau in the past
systemCall1<-function(ndraws, nsam, theta, rho, seqlen, demes, tau, expand, migrationRate, migrationPattern, popSizes=1, ms.output.file="ms.out") {
  
  ms.args.file="./ms.in"
  
  if(length(popSizes!=demes)) {
    popSizes = c(popSizes, 1 + 0*((length(popSizes)+1):demes))
  }

  if(migrationPattern==2) migMatrix<-makeMigMatrix(demes, migrationRate, migrationPattern)
  if(length(nsam)<demes)
    nsam<-c(nsam, 0*((length(nsam)+1):demes))

#  call<-sprintf("../MsDir/ms %d %d -t %g -r %d %d -I %d ", sum(nsam), ndraws, theta, rho, seqlen, demes)
  call <- sprintf("../MsDir/ms %d %d ", sum(nsam), ndraws)
  args <- sprintf("-t %g -r %d %d -I %d ", theta, rho, seqlen, demes)

  for( j in 1:(demes))
    args<-paste(args, sprintf(" %d ", nsam[j]))

  if(demes>1) {
    if(migrationPattern==1) {
      args<-paste(args, sprintf(" %g ", migrationRate))
    } else if(migrationPattern==2) {
      args<-paste(args, sprintf(" -ma "))
      for( i in 1:(demes^2))
        args<-paste(args, sprintf("%g", migMatrix[i]))
    }

    for( j in 1:(demes))
      args<-paste(args, sprintf(" -n %d %g", j, popSizes[j]))


    for( j in 2:demes )
      args<-paste(args, sprintf(" -ej %g %d %d ", tau, j, 1))

    args<-paste(args, sprintf("-eM %g 0", tau))
    args<-paste(args, sprintf("-en %g 1 %g", tau, expand))
  }

  args<-paste(args, sprintf(" -C "))
  if(demes > 100){
    write(args,file=ms.args.file)
    call <- paste(call, " -f ", ms.args.file, " > ", ms.output.file)
  }else{
    call<-paste(call, args, " > ", ms.output.file)
  }

#  cat(call)
  call
}

# make ms command for a model where all demes desccend from a single deme a time tau in the past; fix S
systemCall2<-function(ndraws, nsam, snps, rho, seqlen, demes, tau, expand, migrationRate, migrationPattern, popSizes=1, ms.output.file="ms.out") {
  
  ms.args.file="./ms.in"
  
  if(length(popSizes!=demes)) {
    popSizes = c(popSizes, 1 + 0*((length(popSizes)+1):demes))
  }

  if(migrationPattern==2) migMatrix<-makeMigMatrix(demes, migrationRate, migrationPattern)
  if(length(nsam)<demes)
    nsam<-c(nsam, 0*((length(nsam)+1):demes))

#  call<-sprintf("../MsDir/ms %d %d -s %g -r %d %d -I %d ", sum(nsam), ndraws, snps, rho, seqlen, demes)
  call <- sprintf("../MsDir/ms %d %d ", sum(nsam), ndraws)
  args <- sprintf("-s %g -r %d %d -I %d ", snps, rho, seqlen, demes)
  for( j in 1:(demes))
    args<-paste(args, sprintf(" %d ", nsam[j]))

  if(demes>1) {
    if(migrationPattern==1) {
      args<-paste(args, sprintf(" %g ", migrationRate))
    } else if(migrationPattern==2) {
      args<-paste(args, sprintf(" -ma "))
      for( i in 1:(demes^2))
        args<-paste(args, sprintf("%g", migMatrix[i]))
    }

    for( j in 1:(demes))
      args<-paste(args, sprintf(" -n %d %g", j, popSizes[j]))


    for( j in 2:demes )
      args<-paste(args, sprintf(" -ej %g %d %d ", tau, j, 1))

    args<-paste(args, sprintf("-eM %g 0", tau))
    argsn<-paste(args, sprintf("-en %g 1 %g", tau, expand))
  }

  args<-paste(args, sprintf(" -C "))
  if(demes > 100){
    write(args,file=ms.args.file)
    call <- paste(call, " -f ", ms.args.file, " > ", ms.output.file)
  }else{
    call<-paste(call, args, " > ", ms.output.file)
  }
#  cat(call)
  call
}

systemCallq<-function(ndraws, nsam, theta, rho, seqlen, demes, tau, expand, migrationRate, migrationPattern, snps=0, ms.output.file="ms.out") {
  call<-sprintf("../MsDir/ms %d %d -t %g -r %d %d", nsam, ndraws, theta, rho, seqlen)
  if(snps) call<-paste(call, sprintf("-s %g ", snps))
  call<-paste(call, " > ", ms.output.file)
#  cat(call)
  call
}

qTajimaD<-function(ndraws, nsam, theta, rho, seqlen, p, snps=0, twosided=TRUE) {
  call<-systemCallq(ndraws, nsam, theta, rho, seqlen, snps=snps)
  system(call)
  a<-ms.inp.multi(nsam, ndraws)
  thetaW<-thetaW(a)
  thetaPi<-thetaPi(a)
  tajimaD<-tajimaD(a, thetaW, thetaPi)
  tajimaD<-sort(tajimaD)

  if(twosided==FALSE) {
    tajimaD[as.integer(p*ndraws)]
  } else {
    c(tajimaD[as.integer(0.5*p*ndraws)], tajimaD[as.integer(ndraws-0.5*p*ndraws)])
  }
}

qFuliD<-function(ndraws, nsam, theta, rho, seqlen, p, snps=0, twosided=TRUE) {
  call<-systemCallq(ndraws, nsam, theta, rho, seqlen, snps=snps)
  system(call)
  a<-ms.inp.multi(nsam, ndraws)
  thetaS1<-thetaS1(a)
  fuliD<-fuliD(a, thetaS1)
  fuliD<-sort(fuliD)

  if(twosided==FALSE) {
    fuliD[as.integer(p*ndraws)]
  } else {
    c(fuliD[as.integer(0.5*p*ndraws)], fuliD[as.integer(ndraws-0.5*p*ndraws)])
  }
}

qFayWuH<-function(ndraws, nsam, theta, rho, seqlen, p, twosided=TRUE) {
  call<-systemCallq(ndraws, nsam, theta, rho, seqlen)
  system(call)
  a<-ms.inp.multi(nsam, ndraws)
  fayWuH<-fayWuH(a)
  fayWuH<-sort(fayWuH)
  if(twosided==FALSE) {
    fayWuH[as.integer(p*ndraws)]
  } else {
    c(fayWuH[as.integer(0.5*p*ndraws)], fayWuH[as.integer(ndraws-0.5*p*ndraws)])
  }
}


