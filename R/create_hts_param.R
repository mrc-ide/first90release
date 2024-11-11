#' Create HIV testing rates from parameter inputs
#' 
#' @details
#'
#' This function creates two arrays `hts_rate` and `diagn_rate`.
#' The `hts_rate` array summarizes HIV testing rates for the HIV negative population stratified by:
#' * HIV aggregated age groups (hAG: 1-9)
#' * sex (1 = male; 2 = female)
#' * testing history (1 = never tested; 2 = previously tested negative)
#' * year
#'
#' Array `diagn_rate` summarizes HIV testing rates for the HIV positive population stratified by:
#' * CD4 stage (hDS: 1-7)
#' * HIV age group (hAG: 1-9)
#' * sex (1 = male; 2 = female)
#' * testing and treatment history (1 = never tested; 2 = previously tested negative; 3 = diagnosed, not treated; 4 = on ART)
#' * year
#' @export
create_hts_param <- function(theta, fp) {
  # Incidence of opportunistic infection (as Holmes 2006 and Anglaret 2005, cited in Johson et al. 2015)
  # 1= >500; 2=350-499; 3=250-349; 4=200-249; 5=100-199; 6=50-99; 7= <50.
  # theta_fx <- c(0.05, 0.12, 0.27, 0.27, rep(0.9, 3))
    theta_fx <- c(   0,    0, 0.27, 0.27, rep(0.9, 3)) 
    
  ## We name the parameters to be fitted
  ## Every year new estimates are produced, we need to add one knot
  ## 
  ## # HOTFIX 14 Dec 2020:
  ##   * Adding knots each year makes software incompatible with .shiny90.zip
  ##     outputs saved with previous version.
  ##   * There is not reliable information about which version saved easily
  ##     accessible to what gets passed through here.
  ##   * Hotfix: use the length of the parameter vector to determine which
  ##     year and number of knots.
  ##     * THIS WILL FAIL if other parameters are changed.
  ##       * Added a pretty rigid check on length of parameter vector in the
  ##         block below.
  ##       * Added package test to ensure we know if something changes.
  ##
  ## -- UPDATE HERE --
  max_knot_year <- if (length(theta) == 41) {
                     2019
                   } else if (length(theta) == 43) {
                     2020
                   } else if (length(theta) == 45) {
                     2021
                   } else if (length(theta) == 47) {
                     2022
                   } else if (length(theta) == 49) {
                     2023
                   } else {
                     stop("Unexpected length of parameter vector.")
                   }
  
  knots <- c(1995, 2000:max_knot_year) - fp$ss$proj_start + 1L
  knots_rr_dxunt <- c(2010:max_knot_year) - fp$ss$proj_start + 1L
  
  ## -- UPDATE ABOVE --
  
  n_k1 <- length(knots) - 1 # we remove 1995
  n_k2 <- n_k1 + length(knots_rr_dxunt)
  rate_f <- exp(theta[1:n_k1])
  rr_dxunt <- stats::plogis(theta[(n_k1 + 1):n_k2]) * 8 # Range is 0-8
  rr_m <- stats::plogis(theta[(n_k2 + 1):(n_k2 + 2)]) * 1.1
  rr_test <- 0.95 + stats::plogis(theta[(n_k2 + 3):(n_k2 + 4)]) * 7.05 # Range is 1-8
  rr_plhiv <- 0.05 + stats::plogis(theta[n_k2 + 5]) * (1.95 - 0.05)
  rr_dxart <- stats::plogis(theta[n_k2 + 6])
  rr_25m <- 0.1 + stats::plogis(theta[n_k2 + 7]) * (6 - 0.1)
  rr_35m <- 0.1 + stats::plogis(theta[n_k2 + 8]) * (6 - 0.1)
  rr_25f <- 0.1 + stats::plogis(theta[n_k2 + 9]) * (6 - 0.1)
  rr_35f <- 0.1 + stats::plogis(theta[n_k2 + 10]) * (6 - 0.1)
  pr_oidx <- 0.25 + (stats::plogis(theta[n_k2 + 11]) * (1.75 - 0.25))
  
  # Age function for males and females
  agefn_m <- c(rep(1, 3), rep(rr_25m, 2), rep(rr_35m, 3), rr_35m * 0.8112) # last age group is 50+
  agefn_f <- c(rep(1, 3), rep(rr_25f, 2), rep(rr_35f, 3), rr_35f * 0.8183) # last age group is 50+
  
  # Time trends in testing rates
  fp$t_hts_start <- as.integer(1995 - fp$ss$proj_start + 1L)
  base_rate_f <- stats::approx(knots, c(0, rate_f), seq_len(fp$ss$PROJ_YEARS), rule = 2)$y
  
  knots_rr_m <- c(2005, 2012) - fp$ss$proj_start + 1L
  rr_m <- stats::approx(knots_rr_m, rr_m, seq_len(fp$ss$PROJ_YEARS), rule = 2)$y
  base_rate_m <- base_rate_f * rr_m
  
  # Retest with time 
  knots_rr_test <- c(2005, 2010, 2015) - fp$ss$proj_start + 1L
  rate_rr_test <- exp(stats::approx(knots_rr_test, log(c(1, rr_test)), seq_len(fp$ss$PROJ_YEARS), rule = 2)$y)
  
  ## Testing rate for HIV negative population:
  ## 1: never tested
  ## 2: previously tested negative
  hts_rate <- array(1.0, c(fp$ss$hAG, fp$ss$NG, fp$ss$pDS, fp$ss$PROJ_YEARS))
  
  hts_rate_m <- sweep(hts_rate[,1,,, drop = FALSE], 4, base_rate_m, "*")
  hts_rate_f <- sweep(hts_rate[,2,,, drop = FALSE], 4, base_rate_f, "*")
  hts_rate[,1,,] <- hts_rate_m
  hts_rate[,2,,] <- hts_rate_f
  # We add the age trends
  hts_rate_am <- sweep(hts_rate_m[,,,, drop = FALSE], 1, agefn_m, "*")
  hts_rate_af <- sweep(hts_rate_f[,,,, drop = FALSE], 1, agefn_f, "*")
  hts_rate[,1,,] <- hts_rate_am
  hts_rate[,2,,] <- hts_rate_af
  
  # Relative increase among previously tested negative
  hts_rate[,,2,] <- sweep(hts_rate[,,2,, drop = FALSE], 4, rate_rr_test, "*")

  ## Testing rate for HIV positive population
  ## 1: never tested
  ## 2: previously tested negative
  ## 3: diagnosed
  ## 4: on ART
  diagn_rate <- array(1.0, c(fp$ss$hDS, fp$ss$hAG, fp$ss$NG, 4, fp$ss$PROJ_YEARS))
  
  diagn_rate_m <- sweep(diagn_rate[,,1,,, drop = FALSE], 5, base_rate_m, "*")
  diagn_rate_f <- sweep(diagn_rate[,,2,,, drop = FALSE], 5, base_rate_f, "*")
  diagn_rate[,,1,,] <- diagn_rate_m 
  diagn_rate[,,2,,] <- diagn_rate_f
  # We add the age trends
  diagn_rate_am <- sweep(diagn_rate_m[,,,,, drop = FALSE], 2, agefn_m, "*")
  diagn_rate_af <- sweep(diagn_rate_f[,,,,, drop = FALSE], 2, agefn_f, "*")
  diagn_rate[,,1,,] <- diagn_rate_am
  diagn_rate[,,2,,] <- diagn_rate_af
  
  # Testing rate due to incidence of OI (in never tested and previously tested neg)
  mod_hts <- simmod(fp)
  art_m <- colSums(attr(mod_hts, "artpop")[, , 1:8, 1,, drop = FALSE], , 4)
  hiv_m <- colSums(attr(mod_hts, "hivpop")[, 1:8, 1,, drop = FALSE], , 3)
  oi_m <- pmin(art_m / (art_m + hiv_m) * pr_oidx, 0.95)
  oi_m[is.na(oi_m)] <- 0
  art_f <- colSums(attr(mod_hts, "artpop")[, , 1:8, 2,, drop = FALSE], , 4)
  hiv_f <- colSums(attr(mod_hts, "hivpop")[, 1:8, 2,, drop = FALSE], , 3)
  oi_f <- pmin(art_f / (art_f + hiv_f) * pr_oidx, 0.95)
  oi_f[is.na(oi_f)] <- 0

  diagn_rate_oi <- array(1.0, c(fp$ss$hDS, fp$ss$hAG, fp$ss$NG, 4, fp$ss$PROJ_YEARS))
  diagn_rate_oi[,,1,,] <- sweep(diagn_rate_oi[,,1,,, drop = FALSE], 5, oi_m, "*")
  diagn_rate_oi[,,2,,] <- sweep(diagn_rate_oi[,,2,,, drop = FALSE], 5, oi_f, "*")
  diagn_rate_oi <- sweep(diagn_rate_oi[,,,,], 1, theta_fx, "*")
  
  # Testing among HIV+ never tested (unaware)
  diagn_rate[,,,1,] <- diagn_rate[,,,1,] * rr_plhiv + diagn_rate_oi[,,,1,]
  
  # Relative testing among PLHIV previously tested negative (unaware) 
  diagn_rate[,,,2,] <- sweep(diagn_rate[,,,2,, drop = FALSE], 5, rate_rr_test, "*")
  diagn_rate[,,,2,] <- diagn_rate[,,,2,] * rr_plhiv + diagn_rate_oi[,,,2,]
  
  ## Relative testing among aware
  # Re-testing rate among PLHIV aware (not on ART)
  rate_dxunt <- stats::approx(knots_rr_dxunt, rr_dxunt, seq_len(fp$ss$PROJ_YEARS), rule = 2)$y
  
  diagn_rate_dxunt <- sweep(diagn_rate[,,,3,, drop = FALSE], 5, rate_dxunt, "*")
  diagn_rate[,,,3,] <- diagn_rate_dxunt 
  diagn_rate[,,,3,] <- diagn_rate[,,,3,] + diagn_rate_oi[,,,3,]
  
  ## Relative testing among already on ART (aware)
  diagn_rate_dxart <- sweep(diagn_rate[,,,4,, drop = FALSE], 5, rate_dxunt, "*")
  diagn_rate[,,,4,] <- diagn_rate_dxart * rr_dxart
  
  fp$hts_rate <- hts_rate
  fp$diagn_rate <- diagn_rate
  
  fp
}
