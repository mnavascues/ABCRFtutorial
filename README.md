---
title: Crash Course on Approximate Bayesian Computation in Population Genetics
author: Miguel de Navascués
date: 21 May 2021
---

### Introduction

Tutorial on approximate Bayesian computation. The objective of the tutorial is to provide an insight and an intuition on how ABC works and its potential pitfalls, to give the student the tools to put in practice an ABC analysis and to be able to critically interpret results from its own or published ABC analysis. Working within R, the methodology of this tutorial is to learn from opening up the "black box". The tutorial covers from the classical ABC (i.e. Beaumont local regression algorithm) to the modern use of machine learning approaches (ABC Random Forest).

### Requirements

The code has been tested with R (3.6.3) and RStudio (1.4.1106) on Ubuntu (20.04.2 LTS). It requieres the file `readms.output.R` distributed with the coalescent simulator program [ms](http://home.uchicago.edu/~rhudson1/source/mksamples.html) and the following R packages (tested version on brackets):

* abc (2.1)
* abcrf (1.8.1)
* learnr (0.10.1)
* phyclust (0.1-30)
* tree (1.0-40)
* weights (1.0)

### Instructions

Each learning unit is contained in a RMarkdown file. Opening it with RStudio en clicking on `Run Document` button will generate an interactive file with the lessons and exercises. Currently the are the following units available:

File | Unit
---|----------
0.IntroBayesian.Rmd | Quick Introduction to Bayesian Statistics
00.Simulation.Rmd | Quick Introduction to Simulation
1.ABC.Rmd | The Three Approximations in ABC
2.GoodPractices.Rmd | Good Practices in ABC
3.ClassicalABC.Rmd | Classical ABC
4.Beyond.Rmd | Beyond model choice and parameter estimation
5.ABCRF.Rmd | ABC with Random Forests

(units on ABCRF are in progress)

### todo

* add references on abcrf
* Version in Python of some materials using msprime as simulator, tskit for calculating summary statistics and ranger for abcrf, maybe in a jupyter notebook

