#' Calculate number adn proportion of retests conducted or number tested by year
#' funtion based on number test in the outputs.R script
#' 
#' @param mod model output of class 'eppasm'
#' @param fp parameter inputs (class 'specfp')
#' @param df a data.frame with indices for prediction. See [evertest()] for more information.
#'
#' @return a data.frame consisting of the number of tests, and population size corresponding to rows of df
#'
#' @details
#' 
#' Number of tests (or number tested in past 12 months) are approximated by mid-year
#' counts and annual testing rates within each stratum.
#'
#' @export

number_retests <- function(mod, fp, df){
  tests <- unaware <- aware <- art <- pop <- retests <- numeric(length(df$haidx))
  
  for(i in seq_along(df$haidx)) {
    
    haidx <- df$haidx[i] + 1:df$hagspan[i] - 1
    sidx <- if(df$sidx[i] == 0) 1:2 else df$sidx[i]
    
    paidx <- fp$ss$agfirst.idx[df$haidx[i]] + 1:sum(fp$ss$h.ag.span[haidx]) - 1
    
    if(df$hvidx[i] %in% c(0, 1)){ # testing among HIV-
      
      pop_hivn_ha <- apply(mod[paidx, sidx, fp$ss$hivn.idx, df$yidx[i],
                               drop = FALSE],
                           2:3, fastmatch::ctapply, fp$ss$ag.idx[paidx], sum)
      tested_ha <- attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivn.idx,
                                           df$yidx[i], drop = FALSE]
      
      tests[i] <- tests[i] +
        sum((c(pop_hivn_ha) - c(tested_ha)) * fp$hts_rate[haidx, sidx, 1, df$yidx[i]]) + # among untested HIV- population
        sum(tested_ha * fp$hts_rate[haidx, sidx, 2, df$yidx[i], drop = FALSE])   # among previously tested HIV- population
      retests[i] <- retests[i] +
        sum(tested_ha * fp$hts_rate[haidx, sidx, 2, df$yidx[i], drop = FALSE])
      
      pop[i] <- pop[i] + sum(pop_hivn_ha)
    }
    
    if(df$hvidx[i] %in% c(0, 2)){ # testing among HIV+
      hivpop_ha_hm <- attr(mod, "hivpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1,
                                           sidx, df$yidx[i], drop = FALSE]
      diagn_ha_hm <- attr(mod, "diagnpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1,
                                            sidx, df$yidx[i], drop = FALSE]
      artpop_ha_hm <- colSums(attr(mod, "artpop")[ , , df$haidx[i] + 1:df$hagspan[i] - 1,
                                                   sidx, df$yidx[i], drop = FALSE])
      
      # Proportion ever tested in undiagnosed
      undiagnosed_ha_hm <- hivpop_ha_hm - diagn_ha_hm
      testneg_ha <- attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivp.idx,
                                            df$yidx[i], drop = FALSE]
      prop_testneg <- 1 / colSums(hivpop_ha_hm - diagn_ha_hm) * c(testneg_ha)
      
      naive_ha_hm <- sweep(undiagnosed_ha_hm, 2:4, 1 - prop_testneg, "*")
      testneg_ha_hm <- sweep(undiagnosed_ha_hm, 2:4, prop_testneg, "*")
      
      tests[i] <- tests[i] +
        # Tested for PLHIV never tested
        sum(naive_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 1, df$yidx[i]])) +
        # Tested for PLHIV (unaware) ever tested
        sum(testneg_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 2, df$yidx[i]])) +
        # Test among PLHIV (aware, not on ART)
        sum(diagn_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 3, df$yidx[i]])) +
        # Test among PLHIV (on ART)
        sum(artpop_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 4, df$yidx[i]])) +
        # Late diagnoses
        sum(attr(mod, "late_diagnoses")[, haidx, sidx, df$yidx[i]])
      
      retests[i] <- retests[i] +
        sum(diagn_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 3, df$yidx[i]])) +
        sum(artpop_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 4, df$yidx[i]]))
      
      pop[i] <- pop[i] +
        sum(hivpop_ha_hm) +
        sum(artpop_ha_hm)
      
      art[i] <- art[i] +
        sum(artpop_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 4, df$yidx[i]]))
      
      aware[i] <- aware[i] +
        sum(diagn_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 3, df$yidx[i]])) 
      
      unaware[i] <- unaware[i] +
        sum(naive_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 1, df$yidx[i]])) +
        sum(testneg_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 2, df$yidx[i]])) +
        sum(attr(mod, "late_diagnoses")[, haidx, sidx, df$yidx[i]])
    }
  }
  val <- data.frame(pop = pop, retests = retests, tests = tests,
                    art = art, aware = aware, unaware = unaware)
}

#---- Test retest ----
# HIV-
#' @export
plot_retest_test_neg <- function(mod, fp, likdat, cnt, relative = F, yr_pred = 2018) {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_retest <- expand.grid(year = 2000:end_date, 
                                outcome = "numbertests", 
                                agegr = c("15-99"), 
                                sex = c("both"), 
                                hivstatus = 'negative')
  out_retest <- subset(out_retest, year <= yr_pred)
  
  out_retest$pop <- number_retests(mod, fp,
                                   add_ss_indices(out_retest, fp$ss))$pop
  out_retest$retests <- number_retests(mod, fp,
                                       add_ss_indices(out_retest, fp$ss))$retests
  out_retest$tests <- number_retests(mod, fp,
                                     add_ss_indices(out_retest, fp$ss))$tests
  out_retest$first_test <- out_retest$tests - out_retest$retests
  
  col1 <- rgb(230, 200, 0, 200, max=255)
  col2 <- rgb(50, 100, 155, 200, max=255)

  if(relative == F) {
    plot(out_retest$tests/1000 ~ out_retest$year, type = 'l', col = col1,
         main = expression(bold(paste("Total ", N^o, " of HIV- Tests (in 1,000)"))),
         xlim = c(2000, yr_pred),  ylim = c(0, max(out_retest$tests/1000)*1.2),
         ylab=expression(paste(N^o, " of HIV- Tests (in 1,000)")), xlab = 'Year', lty = 1, lwd = 1)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(rep(0, length = length(out_retest$tests)), rev(out_retest$first_test/1000)),
            col = col1, border = NA)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(out_retest$tests/1000, rev(out_retest$first_test/1000)),
            col = col2, border = NA)
    legend('topleft', inset = 0.02, legend = c('Repeat Testers','First-Time Testers'), 
           col = c(col2, col1), bty = 'n', pch = c(15,15), pt.cex = 2, cex = 0.9)
  } else if(relative == T) {
    plot((out_retest$tests/out_retest$tests)*100 ~ out_retest$year, type='l', ylim=c(0,100), col=col1, 
         main = "Proportion of HIV- Tests \n Among Never and Ever Tested",
         xlim = c(2000, yr_pred), ylab='Proportion of Total HIV- Tests (%)', xlab='Year', lty=1, lwd=1)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c((out_retest$retest/out_retest$tests)*100, rep(100, length = length(out_retest$tests))),
            col = col1, border = NA)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(rep(0, length = length(out_retest$tests)), rev((out_retest$retest/out_retest$tests)*100)),
            col = col2, border = NA)
    legend('topleft', inset = 0.02, legend = c('% Repeat Testers','% First-Time Testers'), 
           col = c(col2, col1), bg = 'grey90', box.col = 'grey90', pch = c(15,15), pt.cex = 2, cex = 0.8)
    } else {
    print('True for relative scale, False for absolute scale')
  }
}

# HIV+
#' @export
plot_retest_test_pos <- function(mod, fp, likdat, cnt, relative = F, yr_pred = 2018, plot_legend = TRUE) {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_retest <- expand.grid(year = 2000:end_date, 
                            outcome = "numbertests", 
                            agegr = c("15-99"), 
                            sex = c("both"), 
                            hivstatus = 'positive')
  out_retest <- subset(out_retest, year <= yr_pred)
  out_retest$pop <- number_retests(mod, fp, add_ss_indices(out_retest, fp$ss))$pop
  out_retest$retests <- number_retests(mod, fp, add_ss_indices(out_retest, fp$ss))$retests
  out_retest$tests <- number_retests(mod, fp, add_ss_indices(out_retest, fp$ss))$tests
  out_retest$art <- number_retests(mod, fp, add_ss_indices(out_retest, fp$ss))$art
  out_retest$aware <- number_retests(mod, fp, add_ss_indices(out_retest, fp$ss))$aware
  out_retest$newdiag <- out_retest$tests - out_retest$retests
  out_retest$newuna <- out_retest$tests - out_retest$art
  
  ylim <- c(0, max(out_retest$tests/1000) * 1.3)
  
  col1 <- rgb(230, 180, 205, 200, max=255)
  col2 <- rgb(150, 0, 50, 200, max=255)
  col3 <- rgb(100, 0, 100, 200, max=255)
    
  if(relative == F){
    plot(out_retest$tests/1000 ~ out_retest$year, type = 'n', ylim = ylim,
         main = "Distribution of HIV+ Tests by \n Awareness and ART Status",
         xlim = c(2000, yr_pred),
         ylab = expression(paste(N^o, " of HIV+ Tests (in 1,000)")), xlab = 'Year')
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(out_retest$newuna/1000, rev(out_retest$tests/1000)),
            col = col1, border = NA)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(out_retest$newdiag/1000, rev(out_retest$newuna/1000)),
            col = col2, border = NA)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(rep(0, length = length(out_retest$tests)), rev(out_retest$newdiag/1000)),
            col = col3, border = NA)
    if (plot_legend) { 
    legend('topleft', inset = 0.02, legend = c('Retests - PLHIV on ART','Retests - Aware Not on ART','New Diagnoses'), 
           col=c(col1, col2, col3), bty = 'n', pch = c(15,15, 15), pt.cex = 2, cex = 0.9) }
  } else if(relative == T) {
    plot((out_retest$tests/out_retest$tests)*100 ~ out_retest$year, type='n', ylim = c(0,100), 
         main = "Distribution of HIV+ Tests by \n Awareness and ART Status",
         xlim = c(2000, yr_pred), ylab = "Proportion of Total HIV+ Tests (%)", xlab = 'Year', lty=1, lwd=1)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c((out_retest$newuna/out_retest$tests)*100, rev((out_retest$tests/out_retest$tests)*100)),
            col = col1, border = NA)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c((out_retest$newdiag/out_retest$tests)*100, rev((out_retest$newuna/out_retest$tests)*100)),
            col = col2, border = NA)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(rep(0, length = length(out_retest$tests)), rev((out_retest$newdiag/out_retest$tests)*100)),
            col = col3, border = NA)
    legend('bottomleft', inset = 0.02, legend = c('% Retests - PLHIV on ART','% Retests - Aware Not on ART','% New Diagnoses'), 
           col = c(col1, col2, col3), bg = 'grey90', box.col = 'grey90', pch = c(15,15, 15), pt.cex = 2, cex = 0.8)
  } else {
    print('TRUE for relative scale; FALSE for absolute scale')
  }
}

#' @export
plot_prv_pos_yld <- function(mod, fp, likdat, cnt, yr_pred = 2018, plot_legend = TRUE) {
  
  start <- fp$ss$proj_start
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  
  yr <- 2000:yr_pred
  prv <- colSums(mod[ , , 2L, yr - start + 1L],,2) / colSums(mod[ , , , yr - start + 1L],,3)
  
  out_test <- expand.grid(year = 2000:end_date, 
                             outcome = "numbertests", 
                             agegr = c("15-99"), 
                             sex = c("both"), 
                             hivstatus = 'all')
  
  out_postest <- expand.grid(year = 2000:end_date, 
                            outcome = "numbertests", 
                            agegr = c("15-99"), 
                            sex = c("both"), 
                            hivstatus = 'positive')
  
  out_test <- subset(out_test, year <= yr_pred)
  out_postest <- subset(out_postest, year <= yr_pred)
  
  # total nb of tests (denominator of yield)
  out_test$tot <- number_tests(mod, fp, add_ss_indices(out_test, fp$ss))$tests
  # total number positive tests (numerator of yield)
  out_postest$tot <- number_tests(mod, fp, add_ss_indices(out_postest, fp$ss))$tests
  
  out_postest$tests <- number_retests(mod, fp,add_ss_indices(out_postest, fp$ss))$tests
  out_postest$retests <- number_retests(mod, fp, add_ss_indices(out_postest, fp$ss))$retests
  out_postest$art <- number_retests(mod, fp, add_ss_indices(out_postest, fp$ss))$art
  out_postest$aware <- number_retests(mod, fp, add_ss_indices(out_postest, fp$ss))$aware
  out_postest$newdiag <- out_postest$tests - out_postest$retests
  out_postest$newuna <- out_postest$tests - out_postest$art
  
  yld <- out_postest$tot / out_test$tot
  ndx <- out_postest$newdiag / out_test$tot
    
  ylim <- c(0, ifelse(max(c(prv, yld, ndx) * 1.3) > 1, 1, max(c(prv, yld, ndx) * 1.3))) * 100
  
  col1 <- rgb(255, 155, 50, 250, max = 255)
  col2 <- rgb(255, 0, 130, 250, max = 255)
  col3 <- rgb(100, 0, 100, 250, max = 255)
  
    plot(I(prv * 100) ~ yr, type = 'l', ylim = ylim, 
         main = "HIV Prevalence (15+), Positivity, and \n Yield of New HIV Diagnoses",
         xlim = c(2000, yr_pred), col = col1, lty = 2, lwd = 1.5, 
         ylab = "Proportion (%)", xlab = 'Year')
    lines(I(yld * 100) ~ out_test$year, col = col2, lwd = 1.5)
    lines(I(ndx * 100) ~ out_postest$year, col = col3, lwd = 1.5)
    if (plot_legend) { 
    legend('topright', inset = 0.01, legend = c('HIV Prevalence (15+)','Positivity','Yield of New Diagnoses'), 
           col=c(col1, col2, col3), bty = 'n', lty = c(2, 1, 1), lwd = 1.5, cex = 0.7) }
 
}
  
  
  

# LATERLATER (not functional function)
plot_retest_test_gr <- function(mod, fp, likdat, cnt, relative = F, yr_pred = 2018) {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_retest <- expand.grid(year=2000:end_date, 
                            outcome = "numbertests", 
                            agegr=c("15-24", "25-34", "35-49"), 
                            sex=c("male", "female"), 
                            hivstatus='negative')
  out_retest$pop <- number_retests(mod, fp,
                                   add_ss_indices(out_retest, fp$ss))$pop
  out_retest$retests <- number_retests(mod, fp,
                                       add_ss_indices(out_retest, fp$ss))$retests
  out_retest$tests <- number_retests(mod, fp,
                                     add_ss_indices(out_retest, fp$ss))$tests
  
  out_1524f <- subset(out_retest, agegr == '15-24' & outcome == 'numbertests' & hivstatus=='negative' & sex=='female')
  out_2534f <- subset(out_retest, agegr == '25-34' & outcome == 'numbertests' & hivstatus=='negative' & sex=='female')
  out_3549f <- subset(out_retest, agegr == '35-49' & outcome == 'numbertests' & hivstatus=='negative' & sex=='female')
  out_1524m <- subset(out_retest, agegr == '15-24' & outcome == 'numbertests' & hivstatus=='negative' & sex=='male')
  out_2534m <- subset(out_retest, agegr == '25-34' & outcome == 'numbertests' & hivstatus=='negative' & sex=='male')
  out_3549m <- subset(out_retest, agegr == '35-49' & outcome == 'numbertests' & hivstatus=='negative' & sex=='male')
  
  if(relative == F) {
    par(mfrow=c(1,2))
    # Females
    plot(out_1524f$tests ~ out_1524f$year, type='l', col='#408000', main="Number of Tests and retests females",
         xlim=c(2000, yr_pred), ylab='Number Tested Among HIV-', xlab='Year', lty=1, lwd=1)
#    polygon(x = c(out_1524f$year, rev(out_1524f$year)),
#            y = c(out_1524f$tests, rev(out_1524f$retests)),
#            col = rgb(211, 211, 211, 100, max = 255),
#            border = NA)
    lines(out_1524f$retests ~ out_1524f$year, col='#408000', lty=2, lwd=1)
    lines(out_2534f$tests ~ out_2534f$year, col='#000080', lty=1, lwd=1)
#    polygon(x = c(out_2534f$year, rev(out_2534f$year)),
#            y = c(out_2534f$tests, rev(out_2534f$retests)),
#            col = rgb(211, 211, 211, 100, max = 255),
#            border = NA)
    lines(out_2534f$retests ~ out_2534f$year, col='#000080', lty=2, lwd=1)
    lines(out_3549f$tests ~ out_3549f$year, col='#804000', lty=1, lwd=1)
#    polygon(x = c(out_3549f$year, rev(out_3549f$year)),
#            y = c(out_3549f$tests, rev(out_3549f$retests)),
#            col = rgb(211, 211, 211, 100, max = 255),
#            border = NA)
    lines(out_3549f$retests ~ out_3549f$year, col='#804000', lty=2, lwd=1)
#    legend('topleft', inset = 0.02, legend=c('Number of tests','Number of retests'), col='navy', lty=1:2, lwd=1, pch=c(15,16))
    # Males
    plot(out_1524m$tests ~ out_1524m$year, type='l', col='#B030A0', main="Number of Tests and retests males",
         xlim=c(2000, yr_pred), ylab='Number Tested Among HIV-', xlab='Year', lty=1, lwd=1)
    polygon(x = c(out_1524m$year, rev(out_1524m$year)),
            y = c(out_1524m$tests, rev(out_1524m$retests)),
            col = rgb(211, 211, 211, 100, max = 255),
            border = NA)
    lines(out_1524m$retests ~ out_1524m$year, col='#B030A0', lty=2, lwd=1)
    lines(out_2534m$tests ~ out_1524m$year, col='#B03060', lty=1, lwd=1)
    polygon(x = c(out_2534m$year, rev(out_2534m$year)),
            y = c(out_2534m$tests, rev(out_2534m$retests)),
            col = rgb(211, 211, 211, 100, max = 255),
            border = NA)
    lines(out_2534m$retests ~ out_2534f$year, col='#B03060', lty=2, lwd=1)
    lines(out_3549m$tests ~ out_3549m$year, col='#B04030', lty=1, lwd=1)
    polygon(x = c(out_3549m$year, rev(out_3549m$year)),
            y = c(out_3549m$tests, rev(out_3549m$retests)),
            col = rgb(211, 211, 211, 100, max = 255),
            border = NA)
    lines(out_3549m$retests ~ out_3549m$year, col='#B04030', lty=2, lwd=1)
#    legend('topleft', inset = 0.02, legend=c('Number of tests','Number of retests'), col='navy', lty=1:2, lwd=1, pch=c(15,16))
  } else if(relative == T) {
    plot((out_retest$tests/out_retest$tests)*100 ~ out_retest$year, type='l', ylim=c(0,100), col='navy', main="Percentage of retests",
         xlim=c(2000, yr_pred), ylab='Proportion Tested Among HIV-', xlab='Year', lty=1, lwd=1)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c((out_retest$tests/out_retest$tests)*100, rev((out_retest$retests/out_retest$tests)*100)),
            col = rgb(211, 211, 211, 100, max = 255),
            border = NA)
    polygon(x = c(out_retest$year, rev(out_retest$year)),
            y = c(rep(0, length = length(out_retest$tests)), rev((out_retest$retests/out_retest$tests)*100)),
            col = rgb(211, 211, 211, 220, max = 255),
            border = NA)
    lines((out_retest$retests/out_retest$tests)*100 ~ out_retest$year, col='navy', lty=2, lwd=1)
    legend('bottomright', inset = 0.02, legend=c('Total of tests','Proportion of retests'), col='navy', lty=1:2, lwd=1, pch=c(15,16))
  } else {
    print('True for relative scale, False for absolute scale')
  }
}
