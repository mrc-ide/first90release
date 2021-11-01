
# Every year new estimates are produced, we need to add one knot
# -- UPDATE HERE --
n_k <- length(2000:2022)
# -- UPDATE ABOVE --

# Starting parameters
theta0 <- c(log(seq(0.001, 0.25, length.out = n_k)), # Females 15-24 baseline testing rate 2000 to 2020
            qlogis(rep(1.5 / 8, n_k - 10)),  # Factor testing among diagnosed, untreated from 2010 to 2020
            qlogis(rep(0.6 / 1.1, 2)), # RR for males in 2005, 2012
            qlogis(rep((1.93 - 0.95) / 7.05, 2)),   # Factor increase among previously tested
            qlogis(0.5),               # RR re-test among PLHIV
            qlogis(0.25),              # Factor testing among diagnosed, on ART
            rep(qlogis(0.9 / 5.9), 4), # Rate ratio for 25-34
            qlogis(0.5))               # Date of inflection point for logistic growth curve OI

usethis::use_data(theta0, overwrite = TRUE)
