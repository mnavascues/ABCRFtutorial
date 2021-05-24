# functions for ABCRF
library(abcrf)
# function for weighted histogram
suppressMessages(library(weights, quietly=T))
# functions for CARTs
library(tree)

target      <- readRDS(file="data/target.RDS")
ref_table_1 <- readRDS(file="data/ref_table_1.RDS")
ref_table_2 <- readRDS(file="data/ref_table_2.RDS")



## Classification and regression trees

#Random forests are a machine learning method for regression and classification. Its name comes from the fact of being constituted of "trees". These trees are *classification and regression trees* or CARTs. So first we will make some exercises to better understand what CARTs are.

### Exercise 1

#Explore CARTs using `tree()` function from `tree` package. First make a regression tree for parameter $\theta_1$ from Tajima's $D$ and $\pi$. Then make a classification tree for model from Tajima's $D$ and Fu and Li's $D$.


# Exercise 1A
################################

regression_tree <- tree(theta1 ~ TD + PI, data=ref_table_2)

plot(regression_tree)
text(regression_tree,cex=0.75)

par(bg = 'aliceblue')
plot(ref_table_2$TD,
     ref_table_2$PI,
     xlab="Tajima's D",
     ylab=expression(pi),
     col=grey(1-ref_table_2$theta1/max(ref_table_2$theta1)),
     pch=20)
partition.tree(regression_tree,
               ordvars=c("TD","PI"),
               add=T,cex=1)

# Exercise 1B
################################

num_of_sims <- 10000
model <- c(rep("constant population size",num_of_sims),rep("population size change",num_of_sims))
sumstats <- rbind(ref_table_1[seq_len(num_of_sims),c("S","PI","NH","TD","FLD")],
                  ref_table_2[seq_len(num_of_sims),c("S","PI","NH","TD","FLD")])
ref_table <-cbind(model,sumstats)

classification_tree <- tree(model ~ TD + FLD, data=ref_table)

plot(classification_tree)
text(classification_tree,cex=0.75)

plot(ref_table$TD[which(model=="population size change")],
     ref_table$FLD[which(model=="population size change")],
     xlab="Tajima's D",
     ylab="Fu and Li's D",
     pch=20)
points(ref_table$TD[which(model=="constant population size")],
       ref_table$FLD[which(model=="constant population size")],
       col="grey",pch=20)
partition.tree(classification_tree,
               ordvars=c("TD","FLD"),
               add=T,cex=1.5,col="red")



## Forest

### Exercise 2

ref_table_1 <- ref_table_1[ref_table_1$PI!=0,] # just for plotting convenience (log scale)

regression_tree <- tree(log10(theta) ~ PI, data=ref_table_1)

plot(regression_tree)
text(regression_tree,cex=0.75)

plot(ref_table_1$PI,
     log10(ref_table_1$theta),
     xlab=expression(pi),
     ylab=expression(log[10]*theta),
     pch=20, log="x")

partition.tree(regression_tree,
               ordvars=c("PI","theta"),
               add=T,cex=1.5,col="red",lwd=2)

for (i in 1:100){
  random_sample <- sample(nrow(ref_table_1),size=500,replace=T)
  ref_table_random_sample <- ref_table_1[random_sample,]
  regression_tree_random_sample <- tree(log10(theta) ~ PI, data=ref_table_random_sample)
  partition.tree(regression_tree_random_sample,
                 ordvars=c("PI","theta"),
                 add=T,cex=1.5,col=7,lwd=1)
}



