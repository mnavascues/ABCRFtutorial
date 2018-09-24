##########################################################################
## Exercise 2:

# Note that it prints on screen the first 10 replicates only
likelihood_approx <- flip.coin.likelihood.approx(n     = total_tosses,
                                                 k     = heads_count,
                                                 p     = p_heads,
                                                 rep   = 1000,
                                                 trace = TRUE)
# Have a look at the code of the function:
flip.coin.likelihood.approx
