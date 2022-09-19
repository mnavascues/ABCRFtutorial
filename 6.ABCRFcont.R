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

# note that in random forest random sub-sampling of the total reference table
# is for simulations and for summary statistics (i.e. rows AND columns)
# this is called Boostrap-AGGgregatING (=BAGGING)
# (bootstrap for the random resampling; aggregating for the averaging among models/trees)

num_of_sims <- 10000

model <- c(rep("constant",num_of_sims),rep("size change",num_of_sims))
sumstats <- rbind(ref_table_1[seq_len(num_of_sims),c("S","PI","NH","TD","FLD")],
                  ref_table_2[seq_len(num_of_sims),c("S","PI","NH","TD","FLD")])
ref_table <-cbind(model,sumstats)

model_RF <- abcrf(formula = model~.,
                  data    = ref_table,
                  lda     = F,
                  ntree   = 1000,
                  paral   = T)
# Variable Importance plot
plot(model_RF, training=ref_table)

# We get an equivalent of cross-validation from this:
# Prior error rate and confusion matrix
model_RF$prior.err
model_RF$model.rf$confusion.matrix

# can the error rate be improved by increasing the number of trees?
err.abcrf(model_RF, training = ref_table, paral = T)

# How are these errors calculated? OUT-OF-BAG
# For each simulation, there is a subset of trees from the forest that have been grown
# without the information from that simulation. These trees are used to estimate the model
# for that simulation. This is the out-of-bag (OOB) estimate.

# out-of-bag estimates:
model_RF$model.rf$predictions[1:20]
# true model:
model[1:20]


# model estimates for octomanati and rinocaracol
model_selection_result_RF <- predict(object= model_RF,
                                     obs = target,
                                     training = ref_table,
                                     ntree = 1000,
                                     paral = T,
                                     paral.predict = T)
(model_selection_result_RF)
 
# What are the votes?
# Why does the posterior probability differ from the proportion of votes?
# Not all trees are equally good, nor all branches of the trees are equally good.
# Proportion of votes is a good for revealing the best model, but it is a bad measure
# of the uncertainty of the decision.
# 
# We can calculate the "local" error (classification error for each simulation) and
# we can grow a second random forest which will learn the relationship
# between the “local” error rate and the summary statistics. The probability
# of a correct classification (posterior probability) is one minus the probability
# of an incorrect classification (error rate). This second random forest allows to
# estimate the posterior probability of the chosen model.
  

