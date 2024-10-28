# --- Individuals functions to plot PJNZ inputs ----
#' @export
get_pjnz_summary_data <- function(fp) {
  mod <- simmod(fp)
  start <- fp$ss$proj_start
  end <- start + fp$ss$PROJ_YEARS - 1L

  list(
    year = start:end,
    pop = apply(mod[,,,], 4, FUN=sum),
    prevalence = attr(mod, "prev15to49"),
    incidence = attr(mod, "incid15to49"),
    plhiv = apply(attr(mod, "hivpop")[,,,], 4, FUN=sum) + apply(attr(mod, "artpop")[,,,,], 5, FUN=sum),
    art_coverage = colSums(fp$art15plus_num)
  )
}

#' @export
#'
## -- UPDATE HERE --
## * update yr_pred to current year
plot_pjnz_prv <- function(pjnz_summary, yr_pred = 2023) {
  pjnz_summary <- stats::na.omit(data.frame(year = pjnz_summary[["year"]], prv = pjnz_summary[["prevalence"]]*100))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  plot(pjnz_summary$prv ~ pjnz_summary$year,
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$prv) * 1.2),
       main = 'HIV Prevalence (15-49 years)', ylab = 'HIV Prevalence (%)', xlab = 'Year', col='grey40')
  graphics::polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$prv)), rev(pjnz_summary$prv)),
          col = grDevices::rgb(155, 0, 50, 150, max = 255), border = NA)
}

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_pjnz_inc <- function(pjnz_summary, yr_pred = 2023) {
  pjnz_summary <- stats::na.omit(data.frame(year = pjnz_summary[["year"]], inc = pjnz_summary[["incidence"]]*1000))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  pjnz_summary <- subset(pjnz_summary, year >= 2000)
  plot(pjnz_summary$inc ~ pjnz_summary$year,
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$inc) * 1.1),
       main='HIV Incidence (15-49 years)', ylab = 'HIV Incidence (per 1000)', xlab='Year', col='grey40')
  graphics::polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$inc)), rev(pjnz_summary$inc)),
          col = grDevices::rgb(255, 155, 0, 200, max = 255), border = NA)
  }

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_pjnz_pop <- function(pjnz_summary, yr_pred = 2023) {
  pjnz_summary <- stats::na.omit(data.frame(year = pjnz_summary[["year"]], pop = pjnz_summary[["pop"]]/1000))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  plot(pjnz_summary$pop ~ pjnz_summary$year,
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$pop) * 1.2),
       main = 'Total Population (15-49 years)', ylab = 'Population (in 1,000s)', xlab='Year', col='grey40')
  graphics::polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$pop)), rev(pjnz_summary$pop)),
          col = grDevices::rgb(90, 170, 240, 150, max = 255), border = NA)
}

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_pjnz_plhiv <- function(pjnz_summary, yr_pred = 2023) {
  pjnz_summary <- stats::na.omit(data.frame(year = pjnz_summary[["year"]], plhiv = pjnz_summary[["plhiv"]]/1000))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  plot(pjnz_summary$plhiv ~ pjnz_summary$year,
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$plhiv) * 1.2),
       main = 'Total Number of PLHIV (15-49 years)', ylab = 'Population of PLHIV (in 1,000s)', xlab='Year', col='grey40')
  graphics::polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$plhiv)), rev(pjnz_summary$plhiv)),
          col = grDevices::rgb(250, 50, 10, 150, max = 255), border = NA)
}


# ---- Single function ----

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_pjnz <- function(fp, yr_pred = 2023) {
  summary <- get_pjnz_summary_data(fp)
  graphics::par(mfrow=c(2,2))
  plot_pjnz_prv(summary, yr_pred)
  plot_pjnz_inc(summary, yr_pred)
  plot_pjnz_pop(summary, yr_pred)
  plot_pjnz_plhiv(summary, yr_pred)
}

combine_rows <- function(prgdat) {

    prg_f <- prgdat[which(prgdat$sex == 'female' & !is.na(prgdat$tot)),]
    prg_m <- prgdat[which(prgdat$sex == 'male' & !is.na(prgdat$tot)),]
    prg_t <- prgdat[which(prgdat$sex == 'both' & !is.na(prgdat$tot)),]

    uni_year <- unique(prg_f$year)
    tot_prg <- NULL
    if (dim(prg_f)[1] > 0) {
        for (i in 1:length(uni_year)) {
            tot_prg.i <- prg_f[prg_f$year == uni_year[i],]
            tot_prg.i$sex <- 'both'
            tot_prg.i$tot <- ifelse(length(prg_f$tot[prg_f$year == uni_year[i]]) == 0, NA,
                                    prg_f$tot[prg_f$year == uni_year[i]]) +
                                ifelse(length(prg_m$tot[prg_m$year == uni_year[i]]) == 0, NA,
                                    prg_m$tot[prg_m$year == uni_year[i]])

            tot_prg.i$totpos <- ifelse(length(prg_f$totpos[prg_f$year == uni_year[i]]) == 0, NA,
                                        prg_f$totpos[prg_f$year == uni_year[i]]) +
                                ifelse(length(prg_m$totpos[prg_m$year == uni_year[i]]) == 0, NA,
                                        prg_m$totpos[prg_m$year == uni_year[i]])

            tot_prg.i[, c('vct','vctpos','anc','ancpos')] <- NA
            tot_prg <- rbind(tot_prg, tot_prg.i)
        }
    }
    to_add <- subset(prg_t, !(year %in% tot_prg$year))
    rbind(tot_prg, to_add)
}

# --- Individuals functions to plot input data ----
#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_input_tot <- function(prgdat, fp, yr_pred = 2023) {
    start <- fp$ss$proj_start
    mod <- simmod(fp)
    pop <- apply(mod[1:35,,,], 4, FUN=sum)

    prg_t <- subset(prgdat, sex == 'both' & !is.na(tot))
    prg_f <- subset(prgdat, sex == 'female' & !is.na(tot))
    prg_m <- subset(prgdat, sex == 'male' & !is.na(tot))

    if (dim(prg_t)[1] > 0 | dim(prg_f)[1] > 0 | dim(prg_m)[1] > 0) {
        tot_prg <- combine_rows(prgdat)
        plot(I(prg_t$tot/1000) ~ prg_t$year,
           ylim = c(0, max(tot_prg$tot/1000, na.rm = TRUE)*1.1),
           xlim = c(2005, yr_pred),
           lwd = 2, col = 'steelblue3', xlab = 'Year', ylab = 'Total number of tests (in 1,000s)', pch = 16,
           main = expression(bold(paste("Total ", N^o, " of Tests Performed (VCT & ANC)"))))
        graphics::lines(I(prg_t$tot/1000) ~ prg_t$year, lwd = 2, col = 'steelblue4')
        p_test_pop <- round(tot_prg$tot/pop[(tot_prg$year - start)]*100, 0)
        p_test_pop_r <- c(min(p_test_pop, na.rm = TRUE), max(p_test_pop, na.rm = TRUE))
        graphics::text(paste('No Test/Pop = ', p_test_pop_r[1], '-', p_test_pop_r[2], '%'),
           x = 2005, y = max(tot_prg$tot/1000, na.rm = TRUE)/20, pos = 4)
        graphics::lines(I(prg_f$tot/1000) ~ prg_f$year, lty = 2, lwd = 2, col = 'steelblue4')
        graphics::points(I(prg_f$tot/1000) ~ prg_f$year, pch = 'f')
        graphics::lines(I(prg_m$tot/1000) ~ prg_m$year, lty = 2, lwd = 2, col = 'steelblue4')
        graphics::points(I(prg_m$tot/1000) ~ prg_m$year, pch = 'm')
  }
}

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_input_totpos <- function(prgdat, fp, yr_pred = 2023) {
    start <- fp$ss$proj_start
    mod <- simmod(fp)
    plhiv <- apply(attr(mod, "hivpop")[,1:8,,], 4, FUN = sum) +
    apply(attr(mod, "artpop")[,,1:8,,], 5, FUN = sum)

    if(sum(prgdat$totpos, na.rm = TRUE) > 0){
        prg_t <- subset(prgdat, sex == 'both' & !is.na(totpos))
        prg_f <- subset(prgdat, sex == 'female' & !is.na(totpos))
        prg_m <- subset(prgdat, sex == 'male' & !is.na(totpos))

        if (dim(prg_t)[1] > 0 | dim(prg_f)[1] > 0 | dim(prg_m)[1] > 0) {

            tot_prg <- combine_rows(prgdat)

            plot(I(prg_t$totpos/1000) ~ prg_t$year,
                ylim = c(0, max(tot_prg$totpos/1000, na.rm = TRUE)*1.1),
                xlim = c(2005, yr_pred), pch=17,
                lwd = 2, col = 'violetred4', xlab = 'Year',
                ylab = 'Total of HIV+ tests (in 1,000s)',
                main = expression(bold(paste("Total ", N^o, " of HIV+ Tests (VCT & ANC)"))))

            graphics::lines(I(prg_t$totpos/1000) ~ prg_t$year, lwd = 2, col = 'violetred4')
            yield <- round(tot_prg$totpos/tot_prg$tot * 100, 1)
            y_range <- c(min(yield, na.rm = TRUE), max(yield, na.rm = TRUE))
            p_pos_pop <- round(tot_prg$totpos/plhiv[(tot_prg$year - start)]*100, 0)
            p_pos_pop_r <-  c(min(p_pos_pop, na.rm = TRUE), max(p_pos_pop, na.rm = TRUE))
            graphics::text(paste('Positivity = ', y_range[1], '-', y_range[2], '%', sep = ''),
                 x = 2005, y = max(tot_prg$totpos/1000, na.rm = TRUE)/6, pos = 4)
            graphics::text(paste('No Positive/PLHIV = ', p_pos_pop_r[1], '-',
                       p_pos_pop_r[2], '%', sep=''), x = 2005,
                 y = max(tot_prg$totpos/1000, na.rm = TRUE)/20, pos = 4)

            graphics::lines(I(prg_f$totpos/1000) ~ prg_f$year, lty = 2, lwd = 2, col = 'violetred4')
            graphics::points(I(prg_f$totpos/1000) ~ prg_f$year, pch = 'f')
            graphics::lines(I(prg_m$totpos/1000) ~ prg_m$year, lty = 2, lwd = 2, col = 'violetred4')
            graphics::points(I(prg_m$totpos/1000) ~ prg_m$year, pch = 'm')
        }
    }
}

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_input_anctot <- function(prgdat, fp, yr_pred = 2023) {

  if (sum(prgdat$anc, na.rm = TRUE) > 0) {
    prgdat <- subset(prgdat, sex != 'male')
      plot(I(prgdat$anc/1000) ~ prgdat$year,
           ylim = c(0, max(prgdat$anc/1000, na.rm = TRUE)*1.1), xlim = c(2005, yr_pred),
           lwd = 2, col = 'darkorange2', xlab = 'Year',
           ylab = 'Total number of ANC tests (in 1,000s)', pch=16,
           main = expression(bold(paste("Total ", N^o, " of HIV Tests Performed at ANC"))))
      graphics::lines(I(prgdat$anc/1000) ~ prgdat$year, lwd = 2, col = 'darkorange3')
    }
}

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_input_ancpos <- function(prgdat, fp, yr_pred = 2023) {

  if (sum(prgdat$ancpos, na.rm = TRUE) > 0) {
    prgdat <- subset(prgdat, sex != 'male')
    plot(I(prgdat$ancpos/1000) ~ prgdat$year,
         ylim = c(0, max(prgdat$ancpos/1000, na.rm = TRUE)*1.1), xlim = c(2005, yr_pred ),
         lwd = 2, col = 'deeppink2', xlab='Year',
         ylab = 'Total number of positive ANC tests (in 1,000s)', pch = 17,
         main = expression(bold(paste("Total ", N^o, " of HIV+ Tests at ANC"))))
    graphics::lines(I(prgdat$ancpos/1000) ~ prgdat$year, lwd = 2, col = 'deeppink2')
  }
}

# ---- Single Function Inputs Data ----

#' @export
## -- UPDATE HERE --
## * update yr_pred to current year
plot_inputdata <- function(prgm_dat, fp, yr_pred = 2023) {
  graphics::par(mfrow = c(2,2))
  plot_input_tot(prgm_dat, fp, yr_pred)
  plot_input_totpos(prgm_dat, fp, yr_pred)
  plot_input_anctot(prgm_dat, fp, yr_pred)
  plot_input_ancpos(prgm_dat, fp, yr_pred)
}


plot_prior <- function(plotprior = TRUE, sd = c(1.625, 1.85, 1.31, 1.68)){
  if (plotprior == TRUE){
  x <- logit(seq(0.00001, 1, 0.001))
  graphics::par(mfrow = c(1,4))
  rr_test <- stats::dnorm(x, mean = logit(0.7/7.2), sd = sd[1])
  rr_plhiv <- stats::dnorm(x, mean = logit(0.5), sd = sd[2])
  rr_dxunt <- stats::dnorm(x, mean = logit(1.5/8), sd = sd[3])
  rr_dxart <- stats::dnorm(x, mean = logit(0.25), sd = sd[4])
  plot(rr_test ~ I(0.8 + stats::plogis(x) * 7.2), type = 'l', xlab = 'RR re-testing', ylab = 'Density', xlim = c(0,8))
    graphics::abline(v = 0.8 + 0.7/7.2 * 7.2, lty = 2, col = 'grey50')
    lci <- round(0.8 + stats::plogis(logit(0.7/7.2) - stats::qnorm(0.975) * sd[1]) * 7.2, 1)
    uci <- round(0.8 + stats::plogis(logit(0.7/7.2) + stats::qnorm(0.975) * sd[1]) * 7.2, 1)
    graphics::text(x = 8, y = max(rr_test), paste('Prior: ', 1.5, '(95%: ', lci, '-', uci, ')'), pos = 2)
  plot(rr_plhiv ~ I(0.05 + stats::plogis(x) * 1.9), type = 'l', xlab = 'RR PLHIV', ylab = 'Density', xlim = c(0,2))
    graphics::abline(v = 0.05 + 0.5 * 1.9, lty = 2, col = 'grey50')
    lci <- round(0.05 + stats::plogis(logit(0.5) - stats::qnorm(0.975) * sd[2]) * 1.9, 1)
    uci <- round(0.05 + stats::plogis(logit(0.5) + stats::qnorm(0.975) * sd[2]) * 1.9, 1)
    graphics::text(x = 2, y = max(rr_plhiv), paste('Prior: ', 1.0, '(95%: ', lci, '-', uci, ')'), pos = 2)
  plot(rr_dxunt ~ I(stats::plogis(x) * 8), type = 'l', xlab = 'RR re-testing for untreated diagnosed', ylab = 'Density', xlim = c(0,8))
    graphics::abline(v = 1.5/8 * 8, lty = 2, col = 'grey50')
    lci <- round(stats::plogis(logit(1.5/8) - stats::qnorm(0.975) * sd[3]) * 8, 1)
    uci <- round(stats::plogis(logit(1.5/8) + stats::qnorm(0.975) * sd[3]) * 8, 1)
    graphics::text(x = 8, y = max(rr_test), paste('Prior: ', 1, '(95%: ', lci, '-', uci, ')'), pos = 2)
  plot(rr_dxart ~ I(stats::plogis(x)), type = 'l', xlab = 'RR re-testing on ART', ylab = 'Density', xlim = c(0,1))
    graphics::abline(v = 0.25, lty = 2, col = 'grey50')
    lci <- round(stats::plogis(logit(0.25) - stats::qnorm(0.975) * sd[4]), 1)
    uci <- round(stats::plogis(logit(0.25) + stats::qnorm(0.975) * sd[4]), 1)
    graphics::text(x = 8, y = max(rr_test), paste('Prior: ', 0.25, '(95%: ', lci, '-', uci, ')'), pos = 2)
  }
}

plot_oi <- function(plotprior = TRUE, par = c(logit(0.5), 0.56)) {
  if (plotprior == TRUE){
    r <- 0.5
    tdx <- seq(1995, 2017, by = 0.1)
    date_oidx <- 2005 + stats::plogis(stats::rnorm(2500, mean = par[1], sd = par[2])) * (2015 - 2005)
    plot(I(stats::rnorm(2)) ~ 1, pch = '', xlab = 'Year', ylab = '% OI Diagnosed',
         ylim = c(0, 1), xlim = c(1995, 2017))
    for (i in 1:length(date_oidx)) {
      pr_oi <- 0.95 / (1 + exp(-r * (tdx - date_oidx[i])))
      graphics::lines(pr_oi ~ tdx, col = grDevices::rgb(150, 100, 200, 25, max = 255), lwd = 1)
    }
    pr_oimin <- 0.95 / (1 + exp(-r * (tdx - 2005)))
    pr_oimax <- 0.95 / (1 + exp(-r * (tdx - 2015)))
    pr_median <-  0.95 / (1 + exp(-r * (tdx - stats::median(date_oidx))))
    graphics::lines(pr_oimin ~ tdx, col = grDevices::rgb(100, 100, 100, 150, max = 255), lwd = 1, lty = 2)
    graphics::lines(pr_oimax ~ tdx, col = grDevices::rgb(100, 100, 100, 150, max = 255), lwd = 1, lty = 2)
    graphics::lines(pr_median ~ tdx, col = 'steelblue4', lwd = 2, lty = 1)
    }
}

#' @export
optimized_par <- function(opt, param = NULL) {
  n_k2 <- length(opt$par) - 11
  n_k1 <- n_k2 - 10
  if (is.null(opt$hessian)) {
    print('You forgot to specify Hessian = TRUE; 95%CrI not calculated.')
    se <- rep(NA, length(opt$par))
    } else {
    # We first test is system is singular (TEMPORARY FIX)
      if (any(diag(opt$hessian) == 0)) {
        index <- which(diag(opt$hessian) == 0)
        opt$hessian[index, index] <- opt$hessian[index, index] - 1e-3
      }
    # From the hessian, simulate the model
    vcova <- solve(-opt$hessian)
    # Test if positive semi-definite
    eS <- eigen(vcova, symmetric = TRUE)
    ev <- eS$values
    if (!all(ev >= -1e-06 * abs(ev[1L]))) {
        vcova <- nearPD(vcova, corr = FALSE)$mat }
    se <- sqrt(diag(as.matrix(vcova))) }

  param_names <- c("RRm_05", "RRm_12", "RR_Test10","RR_Test15",
                   "RR_PLHIV", "RR_DxUnt10",  "RR_DxUnt17",
                   "RR_DxART", "RR_25p_m", "RR_35p_m", "RR_25p_f", "RR_35p_f",
                   "RR OI Dx")

  description <- c('RR testing: men in 2005',
                   'RR testing: men in 2012',
                   'RR re-testing 2010',
                   'RR re-testing 2015',
                   'RR testing: PLHIV unaware',
                   'RR re-testing: PLHIV aware (not ART) 2010',
## -- UPDATE HERE --
## * update label to current year
                   'RR re-testing: PLHIV aware (not ART) 2023',
                   'RR re-testing: PLHIV on ART (*RR not ART)',
                   'RR among 25-34 men',
                   'RR among 35+ men',
                   'RR among 25-34 women',
                   'RR among 35+ women',
                   'RR OI Dx (ART Cov)')

  pt <- c(round(stats::plogis(opt$par[n_k2 + 1]) * 1.1, 2),
          round(stats::plogis(opt$par[n_k2 + 2]) * 1.1, 2),
          round(0.95 + stats::plogis(opt$par[n_k2 + 3]) * 7.05, 2),
          round(0.95 + stats::plogis(opt$par[n_k2 + 4]) * 7.05, 2),
          round(0.05 + stats::plogis(opt$par[n_k2 + 5]) * (1.95 - 0.05), 2),
          round(stats::plogis(opt$par[n_k1 + 1]) * 8, 2),
          round(stats::plogis(opt$par[n_k2 - 1]) * 8, 2),
          round(stats::plogis(opt$par[n_k2 + 6]), 2),
          round(0.1 + stats::plogis(opt$par[n_k2 + 7]) * (6 - 0.1), 2),
          round(0.1 + stats::plogis(opt$par[n_k2 + 8]) * (6 - 0.1), 2),
          round(0.1 + stats::plogis(opt$par[n_k2 + 9]) * (6 - 0.1), 2),
          round(0.1 + stats::plogis(opt$par[n_k2 + 10]) * (6 - 0.1), 2),
          round(0.25 + (stats::plogis(opt$par[n_k2 + 11])) * (1.75 - 0.25), 1))

  if (is.null(param)) {
  lci <- c(round(stats::plogis(opt$par[n_k2 + 1] - stats::qnorm(0.975) * se[n_k2 + 1]) * 1.1, 2),
           round(stats::plogis(opt$par[n_k2 + 2] - stats::qnorm(0.975) * se[n_k2 + 2]) * 1.1, 2),
           round(0.95 + stats::plogis(opt$par[n_k2 + 3] - stats::qnorm(0.975) * se[n_k2 + 3]) * 7.05, 2),
           round(0.95 + stats::plogis(opt$par[n_k2 + 4] - stats::qnorm(0.975) * se[n_k2 + 4]) * 7.05, 2),
           round(0.05 + stats::plogis(opt$par[n_k2 + 5] - stats::qnorm(0.975) * se[n_k2 + 5]) * (1.95 - 0.05), 2),
           round(stats::plogis(opt$par[n_k1 + 1] - stats::qnorm(0.975) * se[n_k1 + 1]) * 8, 2),
           round(stats::plogis(opt$par[n_k2 - 1] - stats::qnorm(0.975) * se[n_k2 - 1]) * 8, 2),
           round(stats::plogis(opt$par[n_k2 + 6] - stats::qnorm(0.975) * se[n_k2 + 6]), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 7] - stats::qnorm(0.975) * se[n_k2 + 7]) * (6 - 0.1), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 8] - stats::qnorm(0.975) * se[n_k2 + 8]) * (6 - 0.1), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 9] - stats::qnorm(0.975) * se[n_k2 + 9]) * (6 - 0.1), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 10] - stats::qnorm(0.975) * se[n_k2 + 10]) * (6 - 0.1), 2),
           round(0.25 + (stats::plogis(opt$par[n_k2 + 11] - stats::qnorm(0.975) * se[n_k2 + 11])) * (1.75 - 0.25), 1))

  uci <- c(round(stats::plogis(opt$par[n_k2 + 1] + stats::qnorm(0.975) * se[n_k2 + 1]) * 1.1, 2),
           round(stats::plogis(opt$par[n_k2 + 2] + stats::qnorm(0.975) * se[n_k2 + 2]) * 1.1, 2),
           round(0.95 + stats::plogis(opt$par[n_k2 + 3] + stats::qnorm(0.975) * se[n_k2 + 3]) * 7.05, 2),
           round(0.95 + stats::plogis(opt$par[n_k2 + 4] + stats::qnorm(0.975) * se[n_k2 + 4]) * 7.05, 2),
           round(0.05 + stats::plogis(opt$par[n_k2 + 5] + stats::qnorm(0.975) * se[n_k2 + 5]) * (1.95 - 0.05), 2),
           round(stats::plogis(opt$par[n_k1 + 1] + stats::qnorm(0.975) * se[n_k1 + 1]) * 8, 2),
           round(stats::plogis(opt$par[n_k2 - 1] + stats::qnorm(0.975) * se[n_k2 - 1]) * 8, 2),
           round(stats::plogis(opt$par[n_k2 + 6] + stats::qnorm(0.975) * se[n_k2 + 6]), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 7] + stats::qnorm(0.975) * se[n_k2 + 7]) * (6 - 0.1), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 8] + stats::qnorm(0.975) * se[n_k2 + 8]) * (6 - 0.1), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 9] + stats::qnorm(0.975) * se[n_k2 + 9]) * (6 - 0.1), 2),
           round(0.1 + stats::plogis(opt$par[n_k2 + 10] + stats::qnorm(0.975) * se[n_k2 + 10]) * (6 - 0.1), 2),
           round(0.25 + (stats::plogis(opt$par[n_k2 + 11] + stats::qnorm(0.975) * se[n_k2 + 11])) * (1.75 - 0.25), 1))
  }
  if (!is.null(param)){
    est <- rbind(stats::plogis(param[, n_k2 + 1]) * 1.1,
            stats::plogis(param[, n_k2 + 2]) * 1.1,
            0.95 + stats::plogis(param[, n_k2 + 3]) * 7.05,
            0.95 + stats::plogis(param[, n_k2 + 4]) * 7.05,
            0.05 + stats::plogis(param[, n_k2 + 5]) * (1.95 - 0.05),
            stats::plogis(param[, n_k1 + 1]) * 8,
            stats::plogis(param[, n_k2 - 1]) * 8,
            stats::plogis(param[, n_k2 + 6]),
            0.1 + stats::plogis(param[, n_k2 + 7]) * (6 - 0.1),
            0.1 + stats::plogis(param[, n_k2 + 8]) * (6 - 0.1),
            0.1 + stats::plogis(param[, n_k2 + 9]) * (6 - 0.1),
            0.1 + stats::plogis(param[, n_k2 + 10]) * (6 - 0.1),
            0.25 + (stats::plogis(param[, n_k2 + 11])) * (1.75 - 0.25))
    lci <- round(apply(est, 1, stats::quantile, probs = 0.025), 2)
    uci <- round(apply(est, 1, stats::quantile, probs = 0.975), 2)

  }
  RR_opt <- data.frame(Parameter_Name = description, Estimate = pt, LCI = lci, UCI = uci)

  return(RR_opt)
}
