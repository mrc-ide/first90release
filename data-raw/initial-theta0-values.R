# Mode of parameters from 2020
par <- read.csv("shiny90-parameter-estimates-2021.csv")
par_med <- apply(par, 1, median)
n_par_2021 <- length(par_med)

n_k_2021 <- length(2000:2021)

par_med_f15to24rate <- par_med[1:n_k_2021]  # Females 15-24 baseline testing rate 2000 to 2021
par_med_diagnrr <- par_med[n_k_2021 + 1:(n_k_2021 - 10)] # Factor testing among diagnosed, untreated from 2010 to 2021
par_med_other <- par_med[n_k_2021*2 - 10 + 1:11]

stopifnot(c(par_med_f15to24rate, par_med_diagnrr, par_med_other) == par_med)
  
# Every year new estimates are produced, we need to add one knot
# -- UPDATE HERE --
n_k <- length(2000:2023)
# -- UPDATE ABOVE --

# Starting parameters
theta0 <- c(par_med_f15to24rate, rep(par_med_f15to24rate[n_k_2021], n_k - n_k_2021),  # Females 15-24 baseline testing rate 2000 to 2020
            par_med_diagnrr, rep(par_med_diagnrr[n_k_2021-10], n_k - n_k_2021), # Factor testing among diagnosed, untreated from 2010 to 2020
            par_med[(n_par_2021 - 10):(n_par_2021 - 9)], # RR for males in 2005, 2012
            par_med[(n_par_2021 - 8):(n_par_2021 - 7)],   # Factor increase among previously tested
            par_med[(n_par_2021 - 6)],               # RR re-test among PLHIV
            par_med[(n_par_2021 - 5)],              # Factor testing among diagnosed, on ART
            par_med[(n_par_2021 - 4):(n_par_2021 - 1)], # Rate ratio for 25-34
            par_med[n_par_2021])                   # Date of inflection point for logistic growth curve OI

# old theta - keeping for comparisons
theta1 <- c(log(seq(0.01, 0.2, length.out = n_k)), # Females 15-24 baseline testing rate 2000 to 2020
            qlogis(rep(1.5 / 8, n_k - 10)),  # Factor testing among diagnosed, untreated from 2010 to 2020
            qlogis(rep(0.6 / 1.1, 2)), # RR for males in 2005, 2012
            qlogis(rep((1.93 - 0.95) / 7.05, 2)),   # Factor increase among previously tested
            qlogis(0.5),               # RR re-test among PLHIV
            qlogis(0.25),              # Factor testing among diagnosed, on ART
            rep(qlogis(0.9 / 5.9), 4), # Rate ratio for 25-34
            qlogis(0.5))   # Date of inflection point for logistic growth curve OI

stopifnot(length(theta0) == n_k + n_k-10 + 11)
stopifnot(length(theta1) == n_k + n_k-10 + 11)

usethis::use_data(theta0, overwrite = TRUE)
write.table(theta0, "theta0-2023.csv", row.names = FALSE, col.names = FALSE)
