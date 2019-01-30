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
plot_pjnz_prv <- function(pjnz_summary, yr_pred = 2018) {
  pjnz_summary <- na.omit(data.frame(year = pjnz_summary[["year"]], prv = pjnz_summary[["prevalence"]]*100))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  plot(pjnz_summary$prv ~ pjnz_summary$year, 
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$prv) * 1.2),
       main = 'HIV Prevalence (15-49 years)', ylab = 'HIV Prevalence (%)', xlab = 'Year', col='grey40')
  polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$prv)), rev(pjnz_summary$prv)),
          col = rgb(155, 0, 50, 150, max = 255), border = NA)
}

#' @export
plot_pjnz_inc <- function(pjnz_summary, yr_pred = 2018) {
  pjnz_summary <- na.omit(data.frame(year = pjnz_summary[["year"]], inc = pjnz_summary[["incidence"]]*1000))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  pjnz_summary <- subset(pjnz_summary, year >= 2000)
  plot(pjnz_summary$inc ~ pjnz_summary$year, 
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$inc) * 1.1),
       main='HIV Incidence (15-49 years)', ylab = 'HIV Incidence (per 1000)', xlab='Year', col='grey40')
  polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$inc)), rev(pjnz_summary$inc)),
          col = rgb(255, 155, 0, 200, max = 255), border = NA)
  }

#' @export
plot_pjnz_pop <- function(pjnz_summary, yr_pred = 2018) {
  pjnz_summary <- na.omit(data.frame(year = pjnz_summary[["year"]], pop = pjnz_summary[["pop"]]/1000))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  plot(pjnz_summary$pop ~ pjnz_summary$year, 
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$pop) * 1.2),
       main = 'Total Population (15-49 years)', ylab = 'Population (in 1,000s)', xlab='Year', col='grey40')
  polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$pop)), rev(pjnz_summary$pop)),
          col = rgb(90, 170, 240, 150, max = 255), border = NA)
}

#' @export
plot_pjnz_plhiv <- function(pjnz_summary, yr_pred = 2018) {
  pjnz_summary <- na.omit(data.frame(year = pjnz_summary[["year"]], plhiv = pjnz_summary[["plhiv"]]/1000))
  pjnz_summary$year <- pjnz_summary$year + 0.5
  plot(pjnz_summary$plhiv ~ pjnz_summary$year, 
       type = 'l', lwd = 2, xlim = c(2000, yr_pred),
       ylim = c(0, max(pjnz_summary$plhiv) * 1.2),
       main = 'Total Number of PLHIV (15-49 years)', ylab = 'Population of PLHIV (in 1,000s)', xlab='Year', col='grey40')
  polygon(x = c(pjnz_summary$year, rev(pjnz_summary$year)),
          y = c(rep(0, length(pjnz_summary$plhiv)), rev(pjnz_summary$plhiv)),
          col = rgb(250, 50, 10, 150, max = 255), border = NA)
}


# ---- Single function ----
plot_pjnz <- function(fp, yr_pred = 2018) {
  summary <- get_pjnz_summary_data(fp)
  par(mfrow=c(2,2))
  plot_pjnz_prv(summary, yr_pred)
  plot_pjnz_inc(summary, yr_pred)
  plot_pjnz_pop(summary, yr_pred)
  plot_pjnz_plhiv(summary, yr_pred)
}


combine_rows <- function(prgdat, col_name){

    prg_b <- prgdat[which(prgdat$sex == "both"), c(col_name,"year")]
    prg_f <- prgdat[which(prgdat$sex == "female"), ]
    prg_m <- prgdat[which(prgdat$sex == "male"), ]

    prg_disagg <- rbind(prg_f, prg_m)

    if (sum(prg_disagg[[col_name]], na.rm = TRUE) > 0) {

        prg_re_agg <- aggregate(get(col_name)~year, data = prg_disagg, function(...){
            if (length(...) != 2){
                NA
            }
            else {
                sum(...)
            }
        })

        prg_re_agg <- setNames(prg_re_agg, c("year", col_name))
        prg_re_agg <- na.omit(prg_re_agg)

        if (nrow(prg_re_agg) > 0) {
            prg_re_agg$sex <- 2
        }
        if (nrow(prg_b) > 0) {
            prg_b$sex <- 1
        }

        prg_t <- rbind(prg_b, prg_re_agg)

        if (nrow(prg_t)>0){
            prg_t <- prg_t[order(prg_t[,'year'],prg_t[,'sex']),]
            prg_t <- prg_t[!duplicated(prg_t$year),]
        }

        prg_t
    }
    else {
        prg_b
    }

}

# --- Individuals functions to plot input data ----
#' @export
plot_input_tot <- function(prgdat, fp, yr_pred = 2018) {
    start <- fp$ss$proj_start
    mod <- simmod(fp)
    pop <- apply(mod[1:35,,,], 4, FUN=  sum)

    prg_t <- combine_rows(prgdat, "tot")

    if (sum(prg_t$tot, na.rm = TRUE) > 0){
        prg_f <- prgdat[which(prgdat$sex == "female"), ]
        prg_m <- prgdat[which(prgdat$sex == "male"), ]

        plot(I(prg_t$tot/1000) ~ prg_t$year,
           ylim = c(0, max(prg_t$tot/1000, na.rm = TRUE)*1.1),
           xlim  = c(2005, yr_pred),
           lwd = 2, col = 'steelblue3', xlab = 'Year', ylab = 'Total number of tests (in 1,000s)', pch = 16,
           main = expression(bold(paste("Total ", N^o, " of Tests Performed (VCT & ANC)"))))

        lines(I(prg_t$tot/1000) ~ prg_t$year, lwd = 2, col = 'steelblue4')
        p_test_pop <- round(prg_t$tot/pop[(prg_t$year - start)]*100, 0)
        p_test_pop_r <- c(min(p_test_pop, na.rm=TRUE), max(p_test_pop, na.rm = TRUE))
        text(paste('No Test/Pop = ', p_test_pop_r[1], '-', p_test_pop_r[2], '%'),
                x = 2005, y = max(prg_t$tot/1000, na.rm = TRUE)/20, pos = 4)

        lines(I(prg_f$tot/1000) ~ prg_f$year, lty = 2, lwd = 2, col = 'steelblue4')
        points(I(prg_f$tot/1000) ~ prg_f$year, pch = 'f')
        lines(I(prg_m$tot/1000) ~ prg_m$year, lty = 2, lwd = 2, col = 'steelblue4')
        points(I(prg_m$tot/1000) ~ prg_m$year, pch = 'm')
    }

}

#' @export
plot_input_totpos <- function(prgdat, fp, yr_pred = 2018) {
    start <- fp$ss$proj_start
    mod <- simmod(fp)
    plhiv <- apply(attr(mod, "hivpop")[,1:8,,], 4, FUN = sum) +
    apply(attr(mod, "artpop")[,,1:8,,], 5, FUN = sum)

    prg_t <- combine_rows(prgdat, "totpos")

    if(sum(prg_t$totpos, na.rm = TRUE) > 0){

        prg_f <- subset(prgdat, sex == 'female')
        prg_m <- subset(prgdat, sex == 'male')

        plot(I(prg_t$totpos/1000) ~ prg_t$year,
             ylim = c(0, max(prg_t$totpos/1000, na.rm=TRUE)*1.1),
             xlim = c(2005, yr_pred), pch=17,
             lwd = 2, col = 'violetred4', xlab = 'Year',
             ylab = 'Total of HIV+ tests (in 1,000s)',
             main = expression(bold(paste("Total ", N^o, " of HIV+ Tests (VCT & ANC)"))))

        lines(I(prg_t$totpos/1000) ~ prg_t$year, lwd = 2, col = 'violetred4')
        yield <- round(prg_t$totpos/prg_t$tot * 100, 1)
        y_range <- c(min(yield, na.rm = TRUE), max(yield, na.rm = TRUE))
        p_pos_pop <- round(prg_t$totpos/plhiv[(prg_t$year - start)]*100, 0)
        p_pos_pop_r <-  c(min(p_pos_pop, na.rm = TRUE), max(p_pos_pop, na.rm = TRUE))
        text(paste('Positivity = ', y_range[1], '-', y_range[2], '%', sep = ''),
             x = 2005, y = max(prg_t$totpos/1000, na.rm = TRUE)/6, pos = 4)
        text(paste('No Positive/PLHIV = ', p_pos_pop_r[1], '-',
                   p_pos_pop_r[2], '%', sep=''), x = 2005,
             y = max(prg_t$totpos/1000, na.rm = TRUE)/20, pos = 4)

        lines(I(prg_f$totpos/1000) ~ prg_f$year, lty = 2, lwd = 2, col = 'violetred4')
        points(I(prg_f$totpos/1000) ~ prg_f$year, pch = 'f')
        lines(I(prg_m$totpos/1000) ~ prg_m$year, lty = 2, lwd = 2, col = 'violetred4')
        points(I(prg_m$totpos/1000) ~ prg_m$year, pch = 'm')
    }
}

#' @export
plot_input_anctot <- function(prgdat, fp, yr_pred = 2018) {

    prgdat <- prgdat[which(prgdat$sex != 'male'),]
    if (sum(prgdat$anc, na.rm = TRUE) > 0) {
        plot(I(prgdat$anc/1000) ~ prgdat$year,
           ylim = c(0, max(prgdat$anc/1000, na.rm = TRUE)*1.1), xlim = c(2005, yr_pred),
           lwd = 2, col = 'darkorange2', xlab = 'Year',
           ylab = 'Total number of ANC tests (in 1,000s)', pch=16,
           main = expression(bold(paste("Total ", N^o, " of HIV Tests Performed at ANC"))))

        lines(I(prgdat$anc/1000) ~ prgdat$year, lwd = 2, col = 'darkorange3')
    }
}

#' @export
plot_input_ancpos <- function(prgdat, fp, yr_pred = 2018) {

    prgdat <- prgdat[which(prgdat$sex != 'male'),]
    if (sum(prgdat$ancpos, na.rm = TRUE) > 0) {
        plot(I(prgdat$ancpos/1000) ~ prgdat$year,
             ylim = c(0, max(prgdat$ancpos/1000, na.rm = TRUE)*1.1), xlim = c(2005, yr_pred ),
             lwd = 2, col = 'deeppink2', xlab='Year',
             ylab = 'Total number of positive ANC tests (in 1,000s)', pch = 17,
             main = expression(bold(paste("Total ", N^o, " of HIV+ Tests at ANC"))))
        lines(I(prgdat$ancpos/1000) ~ prgdat$year, lwd = 2, col = 'deeppink2')
  }
}

# ---- Single Function Inputs Data ----
plot_inputdata <- function(prgm_dat, fp, yr_pred = 2018) {
  par(mfrow = c(2,2))
  plot_input_tot(prgm_dat, fp, yr_pred)
  plot_input_totpos(prgm_dat, fp, yr_pred)
  plot_input_anctot(prgm_dat, fp, yr_pred)
  plot_input_ancpos(prgm_dat, fp, yr_pred) 
}


plot_prior <- function(plotprior = TRUE, sd = c(1.625, 1.85, 1.31, 1.68)){
  if (plotprior == TRUE){
  x <- logit(seq(0.00001, 1, 0.001))
  par(mfrow = c(1,4))
  rr_test <- dnorm(x, mean = logit(0.7/7.2), sd = sd[1])
  rr_plhiv <- dnorm(x, mean = logit(0.5), sd = sd[2])
  rr_dxunt <- dnorm(x, mean = logit(1.5/8), sd = sd[3])
  rr_dxart <- dnorm(x, mean = logit(0.25), sd = sd[4])
  plot(rr_test ~ I(0.8 + plogis(x) * 7.2), type = 'l', xlab = 'RR re-testing', ylab = 'Density', xlim = c(0,8))
    abline(v = 0.8 + 0.7/7.2 * 7.2, lty = 2, col = 'grey50')
    lci <- round(0.8 + plogis(logit(0.7/7.2) - qnorm(0.975) * sd[1]) * 7.2, 1)
    uci <- round(0.8 + plogis(logit(0.7/7.2) + qnorm(0.975) * sd[1]) * 7.2, 1)
    text(x = 8, y = max(rr_test), paste('Prior: ', 1.5, '(95%: ', lci, '-', uci, ')'), pos = 2)
  plot(rr_plhiv ~ I(0.05 + plogis(x) * 1.9), type = 'l', xlab = 'RR PLHIV', ylab = 'Density', xlim = c(0,2))
    abline(v = 0.05 + 0.5 * 1.9, lty = 2, col = 'grey50')
    lci <- round(0.05 + plogis(logit(0.5) - qnorm(0.975) * sd[2]) * 1.9, 1)
    uci <- round(0.05 + plogis(logit(0.5) + qnorm(0.975) * sd[2]) * 1.9, 1)
    text(x = 2, y = max(rr_plhiv), paste('Prior: ', 1.0, '(95%: ', lci, '-', uci, ')'), pos = 2)    
  plot(rr_dxunt ~ I(plogis(x) * 8), type = 'l', xlab = 'RR re-testing for untreated diagnosed', ylab = 'Density', xlim = c(0,8))
    abline(v = 1.5/8 * 8, lty = 2, col = 'grey50')
    lci <- round(plogis(logit(1.5/8) - qnorm(0.975) * sd[3]) * 8, 1)
    uci <- round(plogis(logit(1.5/8) + qnorm(0.975) * sd[3]) * 8, 1)
    text(x = 8, y = max(rr_test), paste('Prior: ', 1, '(95%: ', lci, '-', uci, ')'), pos = 2)
  plot(rr_dxart ~ I(plogis(x)), type = 'l', xlab = 'RR re-testing on ART', ylab = 'Density', xlim = c(0,1))
    abline(v = 0.25, lty = 2, col = 'grey50')
    lci <- round(plogis(logit(0.25) - qnorm(0.975) * sd[4]), 1)
    uci <- round(plogis(logit(0.25) + qnorm(0.975) * sd[4]), 1)
    text(x = 8, y = max(rr_test), paste('Prior: ', 0.25, '(95%: ', lci, '-', uci, ')'), pos = 2)
  }
}

plot_oi <- function(plotprior = TRUE, par = c(logit(0.5), 0.56)) {
  if (plotprior == TRUE){
    r <- 0.5
    tdx <- seq(1995, 2017, by = 0.1)
    date_oidx <- 2005 + plogis(rnorm(2500, mean = par[1], sd = par[2])) * (2015 - 2005)
    plot(I(rnorm(2)) ~ 1, pch = '', xlab = 'Year', ylab = '% OI Diagnosed', 
         ylim = c(0, 1), xlim = c(1995, 2017))
    for (i in 1:length(date_oidx)) {
      pr_oi <- 0.95 / (1 + exp(-r * (tdx - date_oidx[i])))
      lines(pr_oi ~ tdx, col = rgb(150, 100, 200, 25, max = 255), lwd = 1)
    }
    pr_oimin <- 0.95 / (1 + exp(-r * (tdx - 2005)))
    pr_oimax <- 0.95 / (1 + exp(-r * (tdx - 2015)))
    pr_median <-  0.95 / (1 + exp(-r * (tdx - median(date_oidx))))
    lines(pr_oimin ~ tdx, col = rgb(100, 100, 100, 150, max = 255), lwd = 1, lty = 2)
    lines(pr_oimax ~ tdx, col = rgb(100, 100, 100, 150, max = 255), lwd = 1, lty = 2)
    lines(pr_median ~ tdx, col = 'steelblue4', lwd = 2, lty = 1)
    }
}

#' @export
optimized_par <- function(opt, param = NULL){
  require(Matrix)
  if (is.null(opt$hessian)) { 
    print('You forgot to specify Hessian = TRUE; 95%CI not calculated.') 
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
                   'RR re-testing: PLHIV aware (not ART) 2017', 
                   'RR re-testing: PLHIV on ART (*RR not ART)', 
                   'RR among 25-34 men',
                   'RR among 35+ men',
                   'RR among 25-34 women',
                   'RR among 35+ women', 
                   'RR OI Dx (ART Cov)')
  
  pt <- c(round(plogis(opt$par[21]) * 1.1, 2),
          round(plogis(opt$par[22]) * 1.1, 2),
          round(0.95 + plogis(opt$par[23]) * 7.05, 2),
          round(0.95 + plogis(opt$par[24]) * 7.05, 2), 
          round(0.05 + plogis(opt$par[25]) * (1.95 - 0.05), 2), 
          round(plogis(opt$par[26]) * 8, 2),
          round(plogis(opt$par[33]) * 8, 2),
          round(plogis(opt$par[36]), 2),
          round(0.1 + plogis(opt$par[37]) * (6 - 0.1), 2),
          round(0.1 + plogis(opt$par[38]) * (6 - 0.1), 2),
          round(0.1 + plogis(opt$par[39]) * (6 - 0.1), 2),
          round(0.1 + plogis(opt$par[40]) * (6 - 0.1), 2),
          round(0.25 + (plogis(opt$par[41])) * (1.75 - 0.25), 1))
  
  if (is.null(param)) {
  lci <- c(round(plogis(opt$par[21] - qnorm(0.975)*se[21]) * 1.1, 2),
           round(plogis(opt$par[22] - qnorm(0.975)*se[22]) * 1.1, 2),
           round(0.95 + plogis(opt$par[23] - qnorm(0.975)*se[23]) * 7.05, 2),
           round(0.95 + plogis(opt$par[24] - qnorm(0.975)*se[24]) * 7.05, 2),
           round(0.05 + plogis(opt$par[25] - qnorm(0.975)*se[25]) * (1.95 - 0.05), 2), 
           round(plogis(opt$par[26] - qnorm(0.975)*se[26]) * 8, 2),
           round(plogis(opt$par[33] - qnorm(0.975)*se[33]) * 8, 2),
           round(plogis(opt$par[36] - qnorm(0.975)*se[36]), 2),
           round(0.1 + plogis(opt$par[37] - qnorm(0.975)*se[37]) * (6 - 0.1), 2),
           round(0.1 + plogis(opt$par[38] - qnorm(0.975)*se[38]) * (6 - 0.1), 2),
           round(0.1 + plogis(opt$par[39] - qnorm(0.975)*se[39]) * (6 - 0.1), 2),
           round(0.1 + plogis(opt$par[40] - qnorm(0.975)*se[40]) * (6 - 0.1), 2),
           round(0.25 + (plogis(opt$par[41] - qnorm(0.975)*se[41])) * (1.75 - 0.25), 1))
  
  uci <- c(round(plogis(opt$par[21] + qnorm(0.975)*se[21]) * 1.1, 2),
           round(plogis(opt$par[22] + qnorm(0.975)*se[22]) * 1.1, 2),
           round(0.95 + plogis(opt$par[23] + qnorm(0.975)*se[23]) * 7.05, 2),
           round(0.95 + plogis(opt$par[24] + qnorm(0.975)*se[24]) * 7.05, 2),
           round(0.05 + plogis(opt$par[25] + qnorm(0.975)*se[25]) * (1.95 - 0.05), 2), 
           round(plogis(opt$par[26] + qnorm(0.975)*se[26]) * 8, 2),
           round(plogis(opt$par[33] + qnorm(0.975)*se[33]) * 8, 2),
           round(plogis(opt$par[36] + qnorm(0.975)*se[36]), 2),
           round(0.1 + plogis(opt$par[37] + qnorm(0.975)*se[37]) * (6 - 0.1), 2),
           round(0.1 + plogis(opt$par[38] + qnorm(0.975)*se[38]) * (6 - 0.1), 2),
           round(0.1 + plogis(opt$par[39] + qnorm(0.975)*se[39]) * (6 - 0.1), 2),
           round(0.1 + plogis(opt$par[40] + qnorm(0.975)*se[40]) * (6 - 0.1), 2),
           round(0.25 + (plogis(opt$par[41] + qnorm(0.975)*se[41])) * (1.75 - 0.25), 1))
  }
  if (!is.null(param)){
    est <- rbind(plogis(param[,21]) * 1.1,
            plogis(param[,22]) * 1.1,
            0.95 + plogis(param[,23]) * 7.05, 
            0.95 + plogis(param[,24]) * 7.05, 
            0.05 + plogis(param[,25]) * (1.95 - 0.05), 
            plogis(param[,26]) * 8,
            plogis(param[,33]) * 8,
            plogis(param[,36]),
            0.1 + plogis(param[,37]) * (6 - 0.1),
            0.1 + plogis(param[,38]) * (6 - 0.1),
            0.1 + plogis(param[,39]) * (6 - 0.1),
            0.1 + plogis(param[,40]) * (6 - 0.1),
            0.25 + (plogis(param[,41])) * (1.75 - 0.25))
    lci <- round(apply(est, 1, quantile, probs = 0.025), 2)
    uci <- round(apply(est, 1, quantile, probs = 0.975), 2)    
    
  }
  RR_opt <- data.frame(Parameter_Name = description, Estimate = pt, LCI = lci, UCI = uci)
  
  return(RR_opt)
}





