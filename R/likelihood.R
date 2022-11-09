#' @export
ll_hts <- function(theta, fp, likdat) {
  
  fp <- create_hts_param(theta, fp)
  mod <- simmod(fp)
  
  val1 <- ll_evertest(mod, fp, likdat$evertest)
  if (!is.null(likdat$hts)) {
    val2 <- ll_prgdat(mod, fp, likdat$hts)
  } else { val2 <- 0 }
  if (!is.null(likdat$hts_pos)) {
    val3 <- ll_prgdat(mod, fp, likdat$hts_pos)
    } else { val3 <- 0 }
 
  val_prior <- lprior_hts(theta, mod, fp)

  val <- val1 + val2 + val3 + val_prior
  
  return(val)
}

# Likelihood function for ever tested
#' @export
ll_evertest <- function(mod, fp, dat) {
  mu <- evertest(mod, fp, dat)

  ## If projection_period = "calendar", interpolate proportion
  ## ever tested predictions to mid-year.
  if (fp$projection_period == "calendar") {
    dat_last <- dat
    dat_last$yidx <- dat_last$yidx - 1L
    mu_last <- evertest(mod, fp, dat_last)
    mu <- 0.5 * (mu + mu_last)
  }

  if (any(is.na(mu)) || any(mu < 0) || any(mu > 1)) {
    llk <- log(0) 
  } else {
    llk <- sum(dbinom(x = dat$nsuccess, size = dat$neff, prob = mu, log = TRUE))
  }
  
  return(llk)
}

# Likelihood function for number of test
#'@export
ll_prgdat <- function(mod, fp, dat) {
  if (as.character(dat$hivstatus[1]) == 'all') {
    mu <- total_tests(mod, df = add_ss_indices(dat, fp$ss))
    llk <- sum(dnorm(x = dat$tot, mean = mu, sd = dat$l_est_se, log = TRUE))
  } 
  if (as.character(dat$hivstatus[1]) == 'positive') {
    mu <- total_tests(mod, df = add_ss_indices(dat, fp$ss))
    llk <- sum(dnorm(x = dat$tot, mean = mu, sd = dat$l_est_se, log = TRUE))
  }
  if (!(as.character(dat$hivstatus[1]) %in% c('all', 'positive'))) {
    print('Error - HIV status is incorrect'); break
  }
  
  return(llk)
}

## -- UPDATE HERE --
## * max_year = <current_year> incremented each year

art_constraint_penalty <- function(mod, fp, max_year = 2021) {
  ind_year <- c(2000:max_year) - fp$ss$proj_start + 1L
  tot_late <- apply(attr(mod, "late_diagnoses")[,,, ind_year], 4, sum)
  tot_untreated_pop <- apply(attr(mod, "hivpop")[,,, ind_year], 4, sum)
  ratio_late_per_1000_untreated <- tot_late / tot_untreated_pop * 1000
  penalty <- sum(dnorm(x = 0, mean = ratio_late_per_1000_untreated, 
                       sd = 20, log = TRUE))
  return(penalty)
}
# Include this in ll_hts if you want to incorporate the likelihood constraint on ART.
  # val_art_penalty <- art_constraint_penalty(mod, fp, max_year = 2021)
  # val <- val1 + val2 + val3 + val_prior + val_art_penalty
  
# Function to prepare the data for input in the likelihood function.
#' @export
prepare_hts_likdat <- function(dat_evertest, dat_prg, fp) {
  
  # Testing behavior data
  dat_evertest$est <- ifelse(is.na(dat_evertest$est), NA, dat_evertest$est)
  dat_evertest <- dat_evertest[complete.cases(dat_evertest$est),]
  dat_evertest$l_est <- qlogis(dat_evertest$est)
  dat_evertest$l_est_se <- dat_evertest$se / (dat_evertest$est * (1 - dat_evertest$est))
  
  # if using binomial, we calculate Neff from SE
  dat_evertest$neff <- (dat_evertest$est * (1 - dat_evertest$est)) / dat_evertest$se^2 
  dat_evertest$neff <- ifelse(dat_evertest$est < 1e-5, dat_evertest$counts, dat_evertest$neff)
  dat_evertest$neff <- ifelse(dat_evertest$est > 0.999, dat_evertest$counts, dat_evertest$neff)
  dat_evertest$neff <- ifelse(dat_evertest$neff == Inf, dat_evertest$counts, dat_evertest$neff)
  dat_evertest$nsuccess <- round(dat_evertest$est * dat_evertest$neff, 0)
  dat_evertest$neff <- ifelse(round(dat_evertest$neff, 0) < 1, 1, round(dat_evertest$neff, 0))
  dat_evertest <- add_ss_indices(dat_evertest, fp$ss)
  
  # Verifying if input data is OK
  if(any(dat_prg$agegr != '15-99')) { stop(print('Error, HTS program data should be for the 15+ age group only /n
                                                 other age-grouping not supported at the moment'))}
  
  # We remove years with NA from programmatic data
  dat_prg_pos <- dat_prg
  dat_prg_pos <- dat_prg_pos[complete.cases(dat_prg_pos$totpos),]
  dat_prg_pos_sex <- subset(dat_prg_pos, sex != 'both')
  yr_sex <- unique(dat_prg_pos_sex$year)
  dat_prg_pos1 <- subset(dat_prg_pos, sex == 'both' & !(year %in% yr_sex))
  dat_prg_pos2 <- subset(dat_prg_pos, sex != 'both' & (year %in% yr_sex))
  dat_prg_pos <- rbind(dat_prg_pos1, dat_prg_pos2)
  
  dat_prg_pos$tot <- dat_prg_pos$totpos
  
  dat_prg_sex <- subset(dat_prg, sex != 'both')
  yr_sex <- unique(dat_prg_sex$year)
  dat_prg1 <- subset(dat_prg, sex == 'both' & !(year %in% yr_sex))
  dat_prg2 <- subset(dat_prg, sex != 'both' & (year %in% yr_sex))
  dat_prg <- rbind(dat_prg1, dat_prg2)
  dat_prg <- dat_prg[complete.cases(dat_prg$tot),]
  
  # Programmatic data - nb of test (total: VCT + PMTCT + other; ideally among 15+) 
  if (dim(dat_prg)[1] > 0) {
    dat_prg$l_est_se <- ifelse(dat_prg$sex == 'both', dat_prg$tot * 0.05, dat_prg$tot * 0.05 * sqrt(1 + 1 + 2*1)) # the sqrt() is to maintain same SE as if fitting on both sex - assuming perfect correlation
    dat_prg$hivstatus <- 'all'
    dat_prg <- add_ss_indices(dat_prg, fp$ss)
  } else { dat_prg <- NULL }

  # Programmatic data - nb pos tests
  if (dim(dat_prg_pos)[1] > 0) {
    dat_prg_pos$l_est_se <- ifelse(dat_prg_pos$sex == 'both', dat_prg_pos$tot * 0.10, dat_prg_pos$tot * 0.10 *sqrt(1 + 1 + 2*1))
    dat_prg_pos$hivstatus <- 'positive'
    dat_prg_pos <- add_ss_indices(dat_prg_pos, fp$ss)
  } else { dat_prg_pos <- NULL }
  
  return(list(evertest = dat_evertest, hts = dat_prg, hts_pos = dat_prg_pos))
}


lprior_hts <- function(theta, mod, fp) {
## Penalty to smooth testing rates among females aged 15-24 (reference group)
## We calculate penalty for RR of males on the log(rate) scale (and use same SD as for females)

  ## -- UPDATE HERE --
  ## * Extend knots by 1 year to current year
  knots <- 2000:2022 - fp$ss$proj_start + 1L
  ## -- UPDATE ABOVE --
  
  n_k1 <- length(knots)
  n_k2 <- n_k1 * 2 - 10 
    penalty_f <- log(fp$hts_rate[1,2,1, knots[-1]]) - log(fp$hts_rate[1,2,1, knots[-n_k1]])
    penalty_rr_dxunt <- theta[c((n_k1 + 2):n_k2)] - theta[c((n_k1 + 1):(n_k2 - 1))]
    penalty_rr_m <- log(plogis(theta[n_k2 + 2]) * 1.1) - log(plogis(theta[n_k2 + 1]) * 1.1)
    penalty_rr_test <- theta[n_k2 + 4] - theta[n_k2 + 3]

  lprior <- 
    ## Prior for first baseline rate for females # exp(log(0.001) + 1.96*0.25)
    dnorm(x = theta[1], mean = log(0.005), sd = 0.25, log = TRUE) +
    ## Relative testing among PLHIV diagnosed, untreated. 1.50 (95%CI: 0.14-6.00) # plogis(qlogis(1.5/8) - 1.96*1.31)*8
    dnorm(x = theta[n_k2], mean = qlogis(1.5/8), sd = 1.31, log = TRUE) +
    ## Prior for male RR 0.6 (95%CI: 0.07-1.05)# plogis(qlogis(0.6/1.1) + 1.96*1.46) * 1.1
    sum(dnorm(x = theta[n_k2 + 1], mean = qlogis(0.6/1.1), sd = 1.46, log = TRUE)) +
    ## Relative increase among previously tested. 1.93 (95%CI: 1.08-5.00) # 0.95 + plogis(qlogis((1.93 - 0.95)/7.05) - qnorm(0.975)*1.084)*7.05
    dnorm(x = theta[n_k2 + 3], mean = qlogis((1.93 - 0.95)/7.05), sd = 1.084, log = TRUE) +
    ## Relative factor for re-testing among PLHIV. 1.00 (95%CI: 0.10-1.90) # 0.05 + plogis(qlogis(0.95 / 1.90) - 1.96*1.85) * (1.95 - 0.05)
    dnorm(x = theta[n_k2 + 5], mean = qlogis(0.95 / 1.90), sd = 1.85, log = TRUE) +
    ## Relative testing among PLHIV already on ART (95%CI: 0.01-0.90) # plogis(qlogis(0.25) - 1.96*1.68)
    dnorm(x = theta[n_k2 + 6], mean = qlogis(0.25), sd = 1.68, log = TRUE) + 
    ## Prior for age (95% CI is 0.14-5.0)# 0.1 + invlogit(logit(0.9/5.9) - 1.96*1.5) * (6 - 0.1)
    sum(dnorm(x = theta[c((n_k2 + 7):(n_k2 + 10))], mean = qlogis(0.9/5.9), sd = 1.685, log = TRUE)) +
    ## RR OI diangosed for HIV relative to ART coverage 1.0 (0.3-1.7) # 0.25 + plogis(qlogis(0.5) - 1.96*1.75) * (1.75 - 0.25)
    dnorm(x = theta[n_k2 + 11], mean = qlogis(0.5), sd = 1.75, log = TRUE)
  
  prior_sd_f <- 0.205 # exp(log(0.5) - 1.96*0.205); previous of 0.35 when double-counting
  prior_sd_rr_m <- 0.26 # exp(log(0.6) + 1.96*0.26)
  prior_sd_rr_dxunt <- 0.25
  prior_sd_rr_test <- 0.25
  
  hyperprior <- 
    ## Penalty for 1st degree difference (female baseline rate)
    sum(dnorm(x = penalty_f, mean = 0, sd = prior_sd_f, log = TRUE)) + 
    ## Penalty for 1st degree difference (male RR) 
    sum(dnorm(x = penalty_rr_m, mean = 0, sd = prior_sd_rr_m, log = TRUE)) +
    ## Penalty for 1st degree difference (RR_dxunt)  
    sum(dnorm(x = penalty_rr_dxunt, mean = 0, sd = prior_sd_rr_dxunt, log = TRUE)) +
    ## Penalty for 1st degree difference (RR_test)  
    sum(dnorm(x = penalty_rr_test, mean = 0, sd = prior_sd_rr_test, log = TRUE)) 

  return(lprior + hyperprior)
}

