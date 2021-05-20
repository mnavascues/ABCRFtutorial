# tmrca
set.seed(1234)
msout <- ms(nsam = 5, 10, opts = "-t 4 -T")
msout_trees <- unlist(strsplit(msout,"\n"))
msout_trees <- msout_trees[which(msout_trees=="//")+1]

msout_trees
for (i in 1:10){
  print(get.rooted.tree.height(read.tree(text = msout_trees[i])))
}

0.145858258009+0.742465734482

0.631625115871+0.256698846817

0.108708091080+0.147990763187+0.631625115871

((s1: 0.145858258009,s3: 0.145858258009): 0.742465734482,
 (s2: 0.256698846817,(s4: 0.108708091080,s5: 0.108708091080): 0.147990763187): 0.631625115871)

# primary and derived parameters





read.ms.trees <- function( txt=NA, file.ms.output=NA ) {
  
  
}



sample_size = 50
num_of_sim = 1000000
theta_prime <- 10^runif(num_of_sim,min=-1,max=2)
msout <- ms(nsam = sample_size,
            opt = "-t tbs -T",
            tbs.matrix = cbind(theta_prime))
tmrca <- theta_prime*sapply(read.tree(text = msout[grep("//",msout)+1]),get.rooted.tree.height)
msout <- read.ms.output(txt=msout)
S_prime  <- S(msout)
PI_prime <- PI(msout)
NH_prime <- NH(msout)
TD_prime <- TajimaD(msout)
FLD_prime <- FuLiD(msout)
ref_table <- data.frame(theta=theta_prime,
                        tmrca=tmrca,
                        S=S_prime,
                        PI=PI_prime,
                        NH=NH_prime,
                        TD=TD_prime,
                        FLD=FLD_prime)
saveRDS(ref_table,file="ref_table_test.RDS")

plot(log10(tmrca),log10(PI_prime))
plot(log10(theta_prime),log10(PI_prime))


plot(log10(tmrca),log10(S_prime))
plot(log10(theta_prime),log10(S_prime))




cbind(tmrca,theta_prime)

