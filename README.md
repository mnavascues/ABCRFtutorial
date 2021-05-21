---
title: Crash Course on Approximate Bayesian Computation
author: Miguel de Navascués
date: 21 May 2021
---

### Introduction

Tutorial on approximate Bayesian computation. The objective of the tutorial is to provide an insight and an intuition on how ABC works and its potential pitfalls, to give the student the tools to put in practice an ABC analysis and to be able to critically interpret results from its own or published ABC analysis. Working within R, the methodology of this tutorial is to learn from opening up the "black box". The tutorial covers from the classical ABC (i.e. Beaumont local regression algorithm) to the modern use of machine learning approaches (ABC Random Forest).

### Requirements

The code has been tested with R (3.6.3) and RStudio (1.4.1106) on Ubuntu (20.04.2 LTS). It requieres the file `readms.output.R` from the coalescent simulator program [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html) and the following R packages (tested version on brackets):

* abc (2.1)
* abcrf (1.8.1)
* learnr (0.10.1)
* phyclust (0.1-30)
* tree (1.0-40)
* weights (1.0)