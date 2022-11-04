# Mode of parameters from 2020
 par <- read.csv("shiny90-parameter-estimates-2021.csv")
 par_med <- apply(par, 1, median)
 n_par <- length(par_med)

 # Every year new estimates are produced, we need to add one knot
 # -- UPDATE HERE --
 n_k <- length(2000:2023)
 # -- UPDATE ABOVE --

 # Starting parameters
 theta0 <- c(c(par_med[1:(n_k - 1)], par_med[(n_k - 1)]), # Females 15-24 baseline testing rate 2000 to 2020
             c(par_med[n_k:(n_par - 11)], par_med[(n_par - 11)]),  # Factor testing among diagnosed, untreated from 2010 to 2020
             par_med[(n_par - 10):(n_par - 9)], # RR for males in 2005, 2012
             par_med[(n_par - 8):(n_par - 7)],   # Factor increase among previously tested
             par_med[(n_par - 6)],               # RR re-test among PLHIV
             par_med[(n_par - 5)],              # Factor testing among diagnosed, on ART
             par_med[(n_par - 4):(n_par - 1)], # Rate ratio for 25-34
             par_med[n_par])                   # Date of inflection point for logistic growth curve OI

 # old theta - keeping for comparisons
 theta1 <- c(log(seq(0.01, 0.2, length.out = n_k)), # Females 15-24 baseline testing rate 2000 to 2020
             qlogis(rep(1.5 / 8, n_k - 10)),  # Factor testing among diagnosed, untreated from 2010 to 2020
             qlogis(rep(0.6 / 1.1, 2)), # RR for males in 2005, 2012
             qlogis(rep((1.93 - 0.95) / 7.05, 2)),   # Factor increase among previously tested
             qlogis(0.5),               # RR re-test among PLHIV
             qlogis(0.25),              # Factor testing among diagnosed, on ART
             rep(qlogis(0.9 / 5.9), 4), # Rate ratio for 25-34
             qlogis(0.5))   # Date of inflection point for logistic growth curve OI

 usethis::use_data(theta0, overwrite = TRUE)
