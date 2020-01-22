


#  Function to calculate outputs by age / sex / testing history stratifications
#' @export
prb_dx_one_yr <- function(fp, year = c(2000:2019), age = "15-24", sex = "male", test_ever = "never", dt = 0.1, version = "R") {
  
  if (version == "C") {
    val <- prb_dx_one_yr_cpp(fp, year = year, age = age, sex = sex, test_ever = test_ever, dt = dt)
  } 
  
  if (version == "R") {
    n_year <- length(year)
    val <- data.frame(year = year,
                      age = rep(paste(age, collapse = ", "), n_year),
                      sex = rep(paste(sex, collapse = ", "), n_year),
                      prb1yr = NA)
    
    if (age == "15-24") { ind_age <- 1 }
    if (age == "25-34") { ind_age <- 4 }
    if (age == "35-49") { ind_age <- 6 }
    if (age == "50-99") { ind_age <- 9 }
    
    if (sex == "male") { ind_sex <- 1
    } else {
      ind_sex = 2 }
    
    if (test_ever == "never") { ind_th <- 1 }
    if (test_ever == "ever")  { ind_th <- 2 }
    
    # CD4 progression rates
    cd4_prg <- fp$cd4_prog[, ind_age, ind_sex]
    
    # HIV-related death rates
    mort_rate <- fp$cd4_mort[, ind_age, ind_sex]
    
    # Death rates by CD4 category
    mort_rate15 <- fp$cd4_mort[, 1, ind_sex]
    mort_rate25 <- fp$cd4_mort[, 4, ind_sex]
    mort_rate35 <- fp$cd4_mort[, 6, ind_sex]
    mort_rate50 <- fp$cd4_mort[, 9, ind_sex]
    
    # For probability of being Dx before certain time (expressed in dt)
    tind6mo <- round(0.5 / dt, 0)
    tind1yr <- round(1 / dt, 0)
    tind2yr <- round(2 / dt, 0)
    tind5yr <- round(5 / dt, 0)
    
    for (n in 1:n_year) {
      # which group are we calculating this for
      ind_yr <- year[n] - fp$ss$proj_start + 1
      
      # what is the average age of incident HIV infection
      avg_age15 <- 14 + weighted.mean(1:10, w = fp$infections[1:10, ind_sex, ind_yr])
      avg_age25 <- 14 + weighted.mean(11:20, w = fp$infections[11:20, ind_sex, ind_yr])
      avg_age35 <- 14 + weighted.mean(21:35, w = fp$infections[21:35, ind_sex, ind_yr])
      
      # Testing rates by CD4 category
      test_rate15 <- fp$diagn_rate[, 1, ind_sex, ind_th, ind_yr]
      test_rate25 <- fp$diagn_rate[, 4, ind_sex, ind_th, ind_yr]
      test_rate35 <- fp$diagn_rate[, 6, ind_sex, ind_th, ind_yr]
      test_rate50 <- fp$diagn_rate[, 9, ind_sex, ind_th, ind_yr]
      
      # Creating vector for time and testing rates or death rates
      time_int <- seq(0, 35, by = dt)[-1]
      vec_l <- length(time_int)
      
      if (age == '15-24') {
        ind15 <- round((25 - avg_age15) / dt, 0)
        ind25 <- ind15 + 10 / dt
        ind35 <- ind25 + 15 / dt
        test_rate_t <- rbind( do.call("rbind", replicate(ind15, test_rate15, simplify = FALSE)),
                              do.call("rbind", replicate((ind25 - ind15), test_rate25, simplify = FALSE)),
                              do.call("rbind", replicate((ind35 - ind25), test_rate35, simplify = FALSE)),
                              do.call("rbind", replicate((vec_l - ind35), test_rate50, simplify = FALSE))) 
        mort_rate_t <- rbind( do.call("rbind", replicate(ind15, mort_rate15, simplify = FALSE)),
                              do.call("rbind", replicate((ind25 - ind15), mort_rate25, simplify = FALSE)),
                              do.call("rbind", replicate((ind35 - ind25), mort_rate35, simplify = FALSE)),
                              do.call("rbind", replicate((vec_l - ind35), mort_rate50, simplify = FALSE))) 
      }
      if (age == '25-34') {
        ind25 <- round((35 - avg_age25) / dt, 0)
        ind35 <- ind25 + 15 / dt
        test_rate_t <- rbind( do.call("rbind", replicate(ind25, test_rate25, simplify = FALSE)),
                              do.call("rbind", replicate((ind35 - ind25), test_rate35, simplify = FALSE)),
                              do.call("rbind", replicate((vec_l - ind35), test_rate50, simplify = FALSE))) 
        mort_rate_t <- rbind( do.call("rbind", replicate(ind25, mort_rate25, simplify = FALSE)),
                              do.call("rbind", replicate((ind35 - ind25), mort_rate35, simplify = FALSE)),
                              do.call("rbind", replicate((vec_l - ind35), mort_rate50, simplify = FALSE))) 
      }
      if (age == '35-49') {
        ind35 <- round((50 - avg_age35) / dt, 0)
        test_rate_t <- rbind( do.call("rbind", replicate(ind35, test_rate35, simplify = FALSE)),
                              do.call("rbind", replicate((vec_l - ind35), test_rate50, simplify = FALSE))) 
        mort_rate_t <- rbind( do.call("rbind", replicate(ind35, mort_rate35, simplify = FALSE)),
                              do.call("rbind", replicate((vec_l - ind35), mort_rate50, simplify = FALSE))) 
      }
      if (age == '50-99') {
        test_rate_t <- do.call("rbind", replicate(vec_l, test_rate50, simplify = FALSE))
        mort_rate_t <- do.call("rbind", replicate(vec_l, mort_rate50, simplify = FALSE))
      }
      
      # Initialize life table
      time_int_prb_one_year <- seq(0, 2, by = dt)[-1]
      vec_l_prb_one_year <- length(time_int_prb_one_year)
      N <- 1000
      X <- array(data = 0, dim = c(vec_l_prb_one_year, 7))
      nb_dx <- array(data = 0, dim = c(vec_l_prb_one_year, 7))
      X[1, 1] <- N * fp$cd4_initdist[1, ind_age, ind_sex]
      X[1, 2] <- N * fp$cd4_initdist[2, ind_age, ind_sex]
      Xt <- X[1, ]
      # For mean or median time between infection and diagnosis or HIV-death
      for (i in 2:vec_l_prb_one_year) {
        X[i, 1] <- Xt[1] + (-(cd4_prg[1] + mort_rate_t[i, 1] + test_rate_t[i, 1]) * Xt[1]) * dt
        X[i, 2] <- Xt[2] + (cd4_prg[1] * Xt[1] - (cd4_prg[2] + mort_rate_t[i, 2] + test_rate_t[i, 2]) * Xt[2]) * dt
        X[i, 3] <- Xt[3] + (cd4_prg[2] * Xt[2] - (cd4_prg[3] + mort_rate_t[i, 3] + test_rate_t[i, 3]) * Xt[3]) * dt
        X[i, 4] <- Xt[4] + (cd4_prg[3] * Xt[3] - (cd4_prg[4] + mort_rate_t[i, 4] + test_rate_t[i, 4]) * Xt[4]) * dt
        X[i, 5] <- Xt[5] + (cd4_prg[4] * Xt[4] - (cd4_prg[5] + mort_rate_t[i, 5] + test_rate_t[i, 5]) * Xt[5]) * dt
        X[i, 6] <- Xt[6] + (cd4_prg[5] * Xt[5] - (cd4_prg[6] + mort_rate_t[i, 6] + test_rate_t[i, 6]) * Xt[6]) * dt
        X[i, 7] <- Xt[7] + (cd4_prg[6] * Xt[6] - (mort_rate_t[i, 7] + test_rate_t[i, 7]) * Xt[7]) * dt
        Xt <- X[i,]
        
        nb_dx[i, ] <- (test_rate_t[i, ] * X[i - 1, ]) * dt
      }
      
      # Probability of being diagnosed (among those not dying)
      val$prb1yr[val$year == year[n]] <- sum(nb_dx[1:tind1yr, ]) / N
    }
  }
  return(val)
}


# Pooling the 16 stratifications together.
#' @export
pool_prb_dx_one_yr <- function(mod, fp, year = c(2000:2019), 
                            age = c("15-24", "25-34", "35-49", "50-99"),
                            sex = c("male", "female")) {
  
  n_year <- length(year)
  val <- data.frame(year = year, age = rep(paste(age, collapse = "+"), n_year),
                    sex = rep(paste(sex, collapse = "+"), n_year),
                    prb1yr = NA)
  
  for (n in 1:n_year) {
    ind_yr <- year[n] - fp$ss$proj_start + 1
    
    time_mat <- expand.grid(year = year[n], 
                            age = c("15-24", "25-34", "35-49", "50-99"),
                            sex = c("male", "female"), 
                            test_ever = c("never", "ever"))
    time_mat$prb1 <- NA
    time_mat$w <- NA
    index_age <- data.frame(x = c(1, 11, 21, 36),
                            y = c(10, 20, 35, 66))
    
    info_need <- expand.grid(year = year[n],
                             agegr = c("15-24", "25-34", "35-49", "50-99"),
                             sex = c("female", "male"),
                             hivstatus = "negative")
    denom <- denom_lifetable(mod, fp, add_ss_indices(info_need, fp$ss))
    denom$w_n <- denom$never / (denom$never + denom$ever)
    denom$w_e <- denom$ever / (denom$never + denom$ever)
    
    for (i in 1:16) {
      if (time_mat$age[i] == "15-24") { index_agei <- index_age[1, ] }
      if (time_mat$age[i] == "25-34") { index_agei <- index_age[2, ] }
      if (time_mat$age[i] == "35-49") { index_agei <- index_age[3, ] }
      if (time_mat$age[i] == "50-99") { index_agei <- index_age[4, ] }
      if (time_mat$sex[i] == "male") { index_sexi <- 1 }
      if (time_mat$sex[i] == "female") { index_sexi <- 2 }
      
      result_i <- prb_dx_one_yr(fp, year = year[n], age = as.character(time_mat$age[i]),
                                sex = as.character(time_mat$sex[i]), test_ever = as.character(time_mat$test_ever[i]))
      time_mat$prb1[i] <- result_i$prb1yr

      if (time_mat$test_ever[i] == 'never') {
        time_mat$w[i] <- sum(fp$infections[index_agei$x:index_agei$y, index_sexi, ind_yr]) *
          denom$w_n[denom$agegr == time_mat$age[i] & denom$sex == time_mat$sex[i]]
      } else {
        time_mat$w[i] <- sum(fp$infections[index_agei$x:index_agei$y, index_sexi, ind_yr]) *
          denom$w_e[denom$agegr == time_mat$age[i] & denom$sex == time_mat$sex[i]]
      } }
    time_mat_sel <- time_mat[time_mat$age %in% age & time_mat$sex %in% sex, ]
    val$prb1yr[val$year == year[n]] <- weighted.mean(time_mat_sel$prb1, w = time_mat_sel$w)    
  }
  return(val)
}


#' @export
simul_pool_prb_dx_one_yr <- function(samp, mod, fp, year = c(2010:2019),
                               age = c("15-24", "25-34", "35-49", "50-99"),
                               sex = c("male", "female")) {
  
  n_year <- length(year)
  
  n_samp <- nrow(samp)
  val_lst <- list()
  for (i in 1:n_samp) {
    fp_i <- create_hts_param(samp[i, ], fp)
    mod_i <- simmod(fp_i)
    val_lst[[i]] <- pool_prb_dx_one_yr(mod_i, fp_i, 
                        year = year, age = age, sex = sex)
  }
  
  val_ <- data.table::rbindlist(val_lst)
  
  val <- data.frame(year = year, age = rep(paste(age, collapse = "+"), n_year),
                    sex = rep(paste(sex, collapse = "+"), n_year),
                    prb1yr = NA, prb1yr_lci = NA, prb1yr_uci = NA)
  
  for (n in 1:n_year) {
    uncertainty <- quantile(val_$prb1yr[val_$year == year[n]], probs = c(0.5, 0.025, 0.975))
    val$prb1yr[val$year == year[n]] <- uncertainty[1]
    val$prb1yr_lci[val$year == year[n]] <- uncertainty[2]
    val$prb1yr_uci[val$year == year[n]] <- uncertainty[3]
  }
  return(val)
}


# Get the total number of never and ever tested among HIV-
#' @export
denom_lifetable <- function(mod, fp, df_ind) {
  
  ever <- never <- numeric(length(df_ind$haidx))
  
  for(i in seq_along(df_ind$haidx)) {
    
    haidx <- df_ind$haidx[i] + 1:df_ind$hagspan[i] - 1
    sidx <- if(df_ind$sidx[i] == 0) 1:2 else df_ind$sidx[i]
    
    paidx <- fp$ss$agfirst.idx[df_ind$haidx[i]] + 1:sum(fp$ss$h.ag.span[haidx]) - 1
    
    if(df_ind$hvidx[i] %in% c(0, 1)){ # testing among HIV-
      
      pop_hivn_ha <- apply(mod[paidx, sidx, fp$ss$hivn.idx, df_ind$yidx[i],
                               drop = FALSE],
                           2:3, fastmatch::ctapply, fp$ss$ag.idx[paidx], sum)
      tested_ha <- attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivn.idx,
                                           df_ind$yidx[i], drop = FALSE]
      
      never[i] <- sum(c(pop_hivn_ha) - c(tested_ha))  # among untested HIV- population
      ever[i] <- sum(tested_ha)    # among previously tested HIV- population
    }
    val <- data.frame(df_ind[, !(colnames(df_ind) %in% 
                                   c("agestart", "aidx", "agspan", 
                                     "haidx", "hagspan", "sidx", 
                                     "hvidx", "yidx"))], never, ever)
  }
  return(val)
}
