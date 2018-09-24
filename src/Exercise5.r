##########################################################################
## Exercise 5:

seq_data         <- ms.inp.multi(sample_size, 1, ms.output.file="data/dataset2.txt")
target2_S        <- S(seq_data) 
target2_pi       <- thetaPi(seq_data) # 
target2_NH       <- NH(seq_data)
target2_SFS      <- SFS(seq_data)     
target2_TajimasD <- tajimaD(seq_data, thetaW(seq_data), target2_pi)
target2_FayWuH   <- fayWuH(seq_data)
target2_FuLiD    <- fuliD(seq_data, thetaS1(seq_data) )

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
