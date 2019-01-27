#' @export
get_out_evertest <- function(mod, fp,
                             agegr = c("15-24", '25-34','35-49', '15-49'),
                             sex = c("both", "female", "male"),
                             hivstatus = c("all", "negative", "positive") ) {
  
    end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
    out_evertest <- expand.grid(year = 2000:end_date,
                                outcome = "evertest",
                                agegr = agegr,
                                sex = sex,
                                hivstatus = hivstatus)

    out_evertest$value <- evertest(mod, fp, add_ss_indices(out_evertest, fp$ss))
    
    out_evertest
}

#' @export
get_out_nbtest <- function(mod, fp) {
  
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_nbtest <- expand.grid(year = 2000:end_date,
                            outcome = "numbertests",
                            agegrp = "15-99",
                            sex = "both",
                            hivstatus = 'all')
  out_nbtest$hivstatus <- as.character(out_nbtest$hivstatus)
  out_nbtest$value <- total_tests(mod, add_ss_indices(out_nbtest, fp$ss))
  out_nbtest
}

#' @export
get_out_nbtest_sex <- function(mod, fp) {
  
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_nbtest <- expand.grid(year = 2000:end_date,
                            outcome = "numbertests",
                            agegrp = "15-99",
                            sex = c("female", "male"),
                            hivstatus = 'all')
  out_nbtest$hivstatus <- as.character(out_nbtest$hivstatus)
  out_nbtest$value <- total_tests(mod, add_ss_indices(out_nbtest, fp$ss))
  out_nbtest
}

#' @export
get_out_nbtest_pos <- function(mod, fp) {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_nbtest_pos <- expand.grid(year = 2000:end_date,
                                outcome = "numbertests",
                                agegrp = "15-99",
                                sex = "both",
                                hivstatus = 'positive')
  out_nbtest_pos$hivstatus <- as.character(out_nbtest_pos$hivstatus)
  out_nbtest_pos$value <- total_tests(mod, add_ss_indices(out_nbtest_pos, fp$ss))
  
  out_nbtest_pos
}

#' @export
get_out_nbtest_pos_sex <- function(mod, fp) {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_nbtest_pos <- expand.grid(year = 2000:end_date,
                                outcome = "numbertests",
                                agegrp = "15-99",
                                sex = c("female", "male"),
                                hivstatus = 'positive')
  out_nbtest_pos$hivstatus <- as.character(out_nbtest_pos$hivstatus)
  out_nbtest_pos$value <- total_tests(mod, add_ss_indices(out_nbtest_pos, fp$ss))
  
  out_nbtest_pos
}

#' @export
get_out_aware <- function(mod, fp, agegr = '15-49',
                          sex = c("both", "female", "male")) {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_aware <- expand.grid(year = 2000:end_date,
                     outcome = 'aware',
                     agegr = agegr,
                     sex = sex,
                     hivstatus = 'positive')
  out_aware$value <- diagnosed(mod, fp, add_ss_indices(out_aware, fp$ss))
  
  out_aware
}

#' @export
get_out_nbaware <- function(mod, fp, agegr = '15-49',
                          sex = c("both", "female", "male")) {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_nbaware <- expand.grid(year = 2000:end_date,
                           outcome = 'number aware',
                           agegr = agegr,
                           sex = sex,
                           hivstatus = 'positive')
  out_nbaware$value <- nb_aware(mod, fp, add_ss_indices(out_nbaware, fp$ss))
  
  out_nbaware
}

#' @export
get_out_art <- function(mod, fp,
                        gender = 'both') {
  end_date <- fp$ss$proj_start + fp$ss$PROJ_YEARS - 1L
  out_artcov <- data.frame(year = 1970:end_date, outcome = "artcov", agegr = "15-49",
             sex = gender, hivstatus = "positive", value = artcov15to49(mod, sex = gender))
  out_artcov <- subset(out_artcov)
  
  out_artcov
}

#' @export
get_out_pregprev <- function(mod, fp) {
  ## births by age and by HIV aggregation age group
  births_a <- apply(mod[fp$ss$p.fert.idx, fp$ss$f.idx, , ], c(1, 3), sum)
  births_ha <- apply(births_a, 2, fastmatch::ctapply, fp$ss$ag.idx[fp$ss$p.fert.idx], sum)

  ## Weights for proportion of births based on FRRs
  hivneg_w <- apply(mod[fp$ss$p.fert.idx, fp$ss$f.idx, fp$ss$hivn.idx, ], 2, fastmatch::ctapply, fp$ss$ag.idx[fp$ss$p.fert.idx], sum)
  undiagn_w <- colSums((attr(mod, "hivpop") - attr(mod, "diagnpop"))[, fp$ss$h.fert.idx, fp$ss$f.idx, ])
  diagn_w <- colSums(attr(mod, "diagnpop")[, fp$ss$h.fert.idx, fp$ss$f.idx, ])
  artpop_w <- colSums(attr(mod, "artpop")[, , fp$ss$h.fert.idx, fp$ss$f.idx, ],,2)
  
  tot_w <- (hivneg_w + undiagn_w + diagn_w + artpop_w)
  
  ## Number of births by HIV, awareness, and ART status
  hivneg_births <- colSums(births_ha * hivneg_w / tot_w)
  undiagn_births <- colSums(births_ha * undiagn_w / tot_w)
  diagn_births <- colSums(births_ha * diagn_w / tot_w)
  artpop_births <- colSums(births_ha * artpop_w / tot_w)
  
  births <- colSums(births_ha)
  hivp_births <- births - hivneg_births
  
  data.frame(year = fp$ss$proj_start - 1 + seq_len(fp$ss$PROJ_YEARS),
             prev = hivp_births / births,
             aware = (diagn_births + artpop_births) / hivp_births,
             artcov = artpop_births / hivp_births)
}


# ---- Individuals functions to plot model outputs ----
#' @export
plot_out_nbtest <- function(mod, fp, likdat, cnt, simul = NULL, yr_pred = 2018) {
    require(data.table)
  
  # if fitting with HTS program data stratified by sex, we add both sex back
  ld <- likdat$hts
  if (!is.null(ld)) {
    likdat$hts <- aggregate(ld$tot, by = list(ld$year), FUN = sum)
    colnames(likdat$hts) <- c('year','tot') 
  }

  # redact <- c('Namibia','Uganda','Zambia','Zimbabwe')
    redact <- c('XXX')
  
  cnt_to_plot <- ifelse(cnt == "Cote d'Ivoire", "Côte d'Ivoire",
                        ifelse(cnt == "Swaziland", "eSwatini", cnt)) 
  start <- fp$ss$proj_start
  mod <- simmod(fp)
  plhiv <- apply(attr(mod, "hivpop")[,1:8,,], 4, FUN=sum) + 
    apply(attr(mod, "artpop")[,,1:8,,], 5, FUN=sum)
  pop <- apply(mod[1:35,,,], 4, FUN=sum)

  out_nbtest <- subset(get_out_nbtest(mod, fp), year <= yr_pred)
  out_nbtest$year <- out_nbtest$year + 0.5
  
  if (!is.null(ld) & any(ld$sex != 'both')) {
    ld_sex <- subset(ld, sex != 'both')
    ld_sex_pred <- total_tests(mod, df = add_ss_indices(ld_sex, fp$ss))
  }

  
    # Decide if we plot CI or not
    if (!is.null(simul)){
        ci <- getCI(simul$number.test)
        ci <- subset(ci, sex == 'both' & year <= yr_pred)
        ci[is.na(ci)] <- 0
        ci$year <- ci$year + 0.5
        plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    # Total tested
    if (plot.ci == T){
        plot(I(out_nbtest$value/1000) ~ out_nbtest$year, pch = '',
          ylim = c(0, max(out_nbtest$value, na.rm = TRUE)/1000 * 1.25),
          cex = 1, ylab = expression(paste(N^o, " of Tests (in 1,000)")), 
          xlab = 'Year', xlim = c(2000, yr_pred), 
          main = expression(bold(paste("Total ", N^o, " of Tests"))))
        polygon(x = c(ci$year, rev(ci$year)),
          y = c(ci$upper/1000, rev(ci$lower/1000)),
          col = col_ci, border = NA)
        lines(I(out_nbtest$value/1000) ~ out_nbtest$year, lwd = 1, col = 'seagreen3')
        text(x = 2000, y = max(out_nbtest$value, na.rm = TRUE)/1000 * 1.15, cnt_to_plot, cex = 1.25, pos = 4)
        
        if(length(likdat$hts$tot) > 0) {
          if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
          points(I(likdat$hts$tot/1000) ~ I(likdat$hts$year + 0.5), pch = pchpt) 
          p_test_pop <- round(likdat$hts$tot/pop[likdat$hts$year - start]*100, 0)
          p_test_pop_r <- c(min(p_test_pop, na.rm = TRUE), max(p_test_pop, na.rm=TRUE))
          text(paste('No Test/Pop = ', p_test_pop_r[1], '-', p_test_pop_r[2], '%'),
               x = yr_pred, y = max(likdat$hts$tot/1000, na.rm = TRUE)/25, pos = 2)  }
    } else {
      plot(I(out_nbtest$value/1000) ~ out_nbtest$year, pch='',
           ylim = c(0, max(out_nbtest$value, na.rm=TRUE)/1000 + max(out_nbtest$value, na.rm = TRUE)/1000*0.25),
           cex = 1, ylab = expression(paste(N^o, " of Tests (in 1,000)")), 
           xlab = 'Year', xlim=c(2000, yr_pred), 
           main = expression(bold(paste("Total ", N^o, " of Tests"))))
        lines(I(out_nbtest$value/1000) ~ out_nbtest$year, lwd = 1, col = 'seagreen3')
        text(x = 2000, y = max(out_nbtest$value, na.rm = TRUE)/1000, cnt_to_plot, cex = 1.5, pos = 2)

        if(length(likdat$hts$tot) > 0) {
          if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
          points(I(likdat$hts$tot/1000) ~ I(likdat$hts$year + 0.5), pch = pchpt) 
          p_test_pop <- round(likdat$hts$tot/pop[(likdat$hts$year - start)]*100, 0)
          p_test_pop_r <- c(min(p_test_pop, na.rm = TRUE), max(p_test_pop, na.rm=TRUE))
          text(paste('No Test/Pop = ', p_test_pop_r[1], '-', p_test_pop_r[2], '%'),
               x = yr_pred, y = max(likdat$hts$tot/1000, na.rm = TRUE)/20, pos = 2) }
    }
}

plot_out_nbtest_sex <- function(mod, fp, likdat, cnt, simul = NULL, yr_pred = 2018) {
  require(data.table)
  
  # redact <- c('Namibia','Uganda','Zambia','Zimbabwe')
  redact <- c('XXX')
  
  cnt_to_plot <- ifelse(cnt == "Cote d'Ivoire", "Côte d'Ivoire",
                        ifelse(cnt == "Swaziland", "eSwatini", cnt)) 
  start <- fp$ss$proj_start
  mod <- simmod(fp)
  plhiv <- apply(attr(mod, "hivpop")[,1:8,,], 4, FUN=sum) + 
    apply(attr(mod, "artpop")[,,1:8,,], 5, FUN=sum)
  pop <- apply(mod[1:35,,,], 4, FUN=sum)
  
  out_nbtest <- subset(get_out_nbtest_sex(mod, fp), year <= yr_pred)
  out_nbtest$year <- out_nbtest$year + 0.5
  out_nbtest_f <- subset(out_nbtest, sex == 'female')
  out_nbtest_m <- subset(out_nbtest, sex == 'male')
  max_ylim <- max(c(out_nbtest_f$value, out_nbtest_m$value), na.rm = TRUE)
  
  lik_sex <- likdat$hts
  if (!is.null(lik_sex)){
    lik_f <- subset(lik_sex, sex == 'female')
    lik_m <- subset(lik_sex, sex == 'male')
  }
  # Decide if we plot CI or not
  if (!is.null(simul)){
    ci <- getCI(simul$number.test)
    ci <- subset(ci, sex != 'both' & year <= yr_pred)
    ci[is.na(ci)] <- 0
    ci$year <- ci$year + 0.5
    plot.ci <- TRUE
    ci_f <- subset(ci, sex == 'female')
    ci_m <- subset(ci, sex == 'male')
  }
  if (is.null(simul)) plot.ci <- FALSE
  
  col_ci <- rgb(150, 150, 150, 125, max = 255)
  
  # Total tested
  if (plot.ci == T){
    plot(I(out_nbtest_f$value/1000) ~ out_nbtest_f$year, pch = '',
         ylim = c(0, max_ylim /1000 * 1.25),
         cex = 1, ylab = expression(paste(N^o, " of Tests (in 1,000)")), 
         xlab = 'Year', xlim = c(2000, yr_pred), 
         main = expression(bold(paste("Total ", N^o, " of Tests"))))
    polygon(x = c(ci_f$year, rev(ci_f$year)),
            y = c(ci_f$upper/1000, rev(ci_f$lower/1000)),
            col = col_ci, border = NA)
    polygon(x = c(ci_m$year, rev(ci_m$year)),
            y = c(ci_m$upper/1000, rev(ci_m$lower/1000)),
            col = col_ci, border = NA)
    lines(I(out_nbtest_f$value/1000) ~ out_nbtest_f$year, lwd = 1, col = 'seagreen3')
    lines(I(out_nbtest_m$value/1000) ~ out_nbtest_m$year, lwd = 1, col = 'seagreen3')
    text(x = 2000, y = max(out_nbtest$value, na.rm = TRUE)/1000 * 1.15, cnt_to_plot, cex = 1.25, pos = 4)
    
    if(length(lik_f$tot) > 0 | length(lik_m$tot) > 0) {
      if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
      points(I(lik_f$tot/1000) ~ I(lik_f$year + 0.5), pch = pchpt) 
      points(I(lik_m$tot/1000) ~ I(lik_m$year + 0.5), pch = pchpt + 1) 
 }
  } else {
    plot(I(out_nbtest_f$value/1000) ~ out_nbtest_f$year, pch = '',
         ylim = c(0, max_ylim /1000 * 1.25),
         cex = 1, ylab = expression(paste(N^o, " of Tests (in 1,000)")), 
         xlab = 'Year', xlim = c(2000, yr_pred), 
         main = expression(bold(paste("Total ", N^o, " of Tests"))))
    lines(I(out_nbtest_f$value/1000) ~ out_nbtest_f$year, lwd = 1, col = 'seagreen3')
    lines(I(out_nbtest_m$value/1000) ~ out_nbtest_m$year, lwd = 1, col = 'seagreen3')
    text(x = 2000, y = max(out_nbtest$value, na.rm = TRUE)/1000, cnt_to_plot, cex = 1.5, pos = 2)
    
    if(length(likdat$hts$tot) > 0) {
      if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
      points(I(lik_f$tot/1000) ~ I(lik_f$year + 0.5), pch = pchpt) 
      points(I(lik_m$tot/1000) ~ I(lik_m$year + 0.5), pch = pchpt + 1)  }
  }
}


#' @export
plot_out_nbpostest <- function(mod, fp, likdat, cnt, simul = NULL, yr_pred = 2018) {
  require(data.table)
  
  # if fitting with HTS program data stratified by sex, we add both sex back
  ld <- likdat$hts
  if (!is.null(ld)) {
    likdat$hts <- aggregate(ld$tot, by = list(ld$year), FUN = sum)
    colnames(likdat$hts) <- c('year','tot')
    likdat$hts$totpos <- aggregate(ld$totpos, by = list(ld$year), FUN = sum)$x
  }
  ldp <- likdat$hts_pos
  if (!is.null(ldp)) {
  likdat$hts_pos <- aggregate(ldp$tot, by = list(ldp$year), FUN = sum)
  colnames(likdat$hts_pos) <- c('year','tot')
  }
  
  # redact <- c('Namibia','Uganda','Zambia','Zimbabwe')
  redact <- c('XXX')
  
  start <- fp$ss$proj_start
  mod <- simmod(fp)
  plhiv <- apply(attr(mod, "hivpop")[,1:8,,], 4, FUN=sum) + 
    apply(attr(mod, "artpop")[,,1:8,,], 5, FUN=sum)
  pop <- apply(mod[1:35,,,], 4, FUN=sum)

    out_nbtest_pos <- get_out_nbtest_pos(mod, fp)
    out_nbtest_pos <- subset(out_nbtest_pos, year <= yr_pred)
    out_nbtest_pos$year <- out_nbtest_pos$year + 0.5  
    
    # Decide if we plot CI or not
    if (!is.null(simul)){
        ci <- getCI(simul$number.test.pos)
        ci <- subset(ci, sex == 'both' & year <= yr_pred)
        ci[is.na(ci)] <- 0
        ci$year <- ci$year + 0.5
        plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    # Postive tests
    if (plot.ci == T){
        plot(I(out_nbtest_pos$value/1000) ~ out_nbtest_pos$year, pch = '',
          ylim = c(0, max(out_nbtest_pos$value, na.rm = TRUE)/1000 * 1.25),
          cex = 1, ylab = expression(paste(N^o, " of Positive Tests (in 1,000)")), 
          xlab='Year', xlim = c(2000, yr_pred), 
          main = expression(bold(paste("Total ", N^o, " of Positive Tests"))))
        polygon(x = c(ci$year, rev(ci$year)),
          y = c(ci$upper/1000, rev(ci$lower/1000)),
          col = col_ci, border = NA)
        lines(I(out_nbtest_pos$value/1000) ~ out_nbtest_pos$year, lwd = 1, col = 'orangered2')
        if (!is.null(likdat$hts_pos$tot)) {
          if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
          points(I(likdat$hts_pos$tot/1000) ~ I(likdat$hts_pos$year + 0.5), pch = pchpt) 
          yield <- round(likdat$hts$totpos / likdat$hts$tot * 100, 1)
          y_range <- c(min(yield, na.rm = TRUE), max(yield, na.rm = TRUE))
          p_pos_pop <- round(likdat$hts_pos$tot / plhiv[(likdat$hts_pos$year - start)] * 100, 0)
          p_pos_pop_r <-  c(min(p_pos_pop, na.rm = TRUE), max(p_pos_pop, na.rm = TRUE))
          text(paste('Positivity = ', y_range[1], '-', y_range[2], '%', sep=''), 
               x = yr_pred, y = max(likdat$hts_pos$tot/1000, na.rm = TRUE)/6, pos = 2)
          text(paste('No Positive/PLHIV = ', p_pos_pop_r[1], '-', 
                     p_pos_pop_r[2], '%', sep=''), x = yr_pred, 
               y = max(likdat$hts_pos$tot/1000, na.rm = TRUE)/20, pos = 2)
        }
    } else {
        plot(I(out_nbtest_pos$value/1000) ~ out_nbtest_pos$year, pch = '',
          ylim = c(0, max(out_nbtest_pos$value, na.rm=TRUE)/1000 * 1.25),
          cex = 1, ylab = expression(paste(N^o, " of Positive Tests (in 1,000)")), 
          xlab = 'Year', xlim = c(2000, yr_pred), 
          main = expression(bold(paste("Total ", N^o, " of Positive Tests"))))
        lines(I(out_nbtest_pos$value/1000) ~ out_nbtest_pos$year, lwd = 1, col = 'orangered2')
        if (!is.null(likdat$hts_pos$tot)) {
          if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
            points(I(likdat$hts_pos$tot/1000) ~ I(likdat$hts_pos$year + 0.5), pch = pchpt) 
          yield <- round(likdat$hts$totpos / likdat$hts$tot * 100, 1)
          y_range <- c(min(yield, na.rm = TRUE), max(yield, na.rm = TRUE))
          p_pos_pop <- round(likdat$hts_pos$tot / plhiv[(likdat$hts_pos$year - start)] * 100, 0)
          p_pos_pop_r <-  c(min(p_pos_pop, na.rm = TRUE), max(p_pos_pop, na.rm = TRUE))
          text(paste('Positivity = ', y_range[1], '-', y_range[2], '%', sep=''), 
               x = yr_pred, y = max(likdat$hts_pos$tot/1000, na.rm = TRUE)/6, pos = 2)
          text(paste('No Positive/PLHIV = ', p_pos_pop_r[1], '-', 
                     p_pos_pop_r[2], '%', sep=''), x = yr_pred, 
               y = max(likdat$hts_pos$tot/1000, na.rm = TRUE)/20, pos = 2)
        }
    }
}

#' @export
plot_out_nbpostest_sex <- function(mod, fp, likdat, cnt, simul = NULL, yr_pred = 2018) {
  require(data.table)
  
  # redact <- c('Namibia','Uganda','Zambia','Zimbabwe')
  redact <- c('XXX')
  
  start <- fp$ss$proj_start
  mod <- simmod(fp)
  plhiv <- apply(attr(mod, "hivpop")[,1:8,,], 4, FUN=sum) + 
    apply(attr(mod, "artpop")[,,1:8,,], 5, FUN=sum)
  pop <- apply(mod[1:35,,,], 4, FUN=sum)
  
  out_nbtest_pos <- subset(get_out_nbtest_pos_sex(mod, fp), year <= yr_pred)
  out_nbtest_pos$year <- out_nbtest_pos$year + 0.5
  out_nbtest_pos_f <- subset(out_nbtest_pos, sex == 'female')
  out_nbtest_pos_m <- subset(out_nbtest_pos, sex == 'male')
  max_ylim <- max(c(out_nbtest_pos_f$value, out_nbtest_pos_m$value), na.rm = TRUE)
  
  lik_sex <- likdat$hts_pos
  if (!is.null(lik_sex)){
    lik_f <- subset(lik_sex, sex == 'female')
    lik_m <- subset(lik_sex, sex == 'male')
  }
  # Decide if we plot CI or not
  if (!is.null(simul)){
    ci <- getCI(simul$number.test.pos)
    ci <- subset(ci, sex != 'both' & year <= yr_pred)
    ci[is.na(ci)] <- 0
    ci$year <- ci$year + 0.5
    plot.ci <- TRUE
    ci_f <- subset(ci, sex == 'female')
    ci_m <- subset(ci, sex == 'male')
  }
  if (is.null(simul)) plot.ci <- FALSE
  
  col_ci <- rgb(150, 150, 150, 125, max = 255)
  
  # Postive tests
  if (plot.ci == T){
    plot(I(out_nbtest_pos_f$value/1000) ~ out_nbtest_pos_f$year, pch = '',
         ylim = c(0, max_ylim/1000 * 1.25),
         cex = 1, ylab = expression(paste(N^o, " of Positive Tests (in 1,000)")), 
         xlab = 'Year', xlim = c(2000, yr_pred), 
         main = expression(bold(paste("Total ", N^o, " of Positive Tests"))))
    polygon(x = c(ci_f$year, rev(ci_f$year)),
            y = c(ci_f$upper/1000, rev(ci_f$lower/1000)),
            col = col_ci, border = NA)
    polygon(x = c(ci_m$year, rev(ci_m$year)),
            y = c(ci_m$upper/1000, rev(ci_m$lower/1000)),
            col = col_ci, border = NA)
    lines(I(out_nbtest_pos_f$value/1000) ~ out_nbtest_pos_f$year, lwd = 1, col = 'orangered2')
    lines(I(out_nbtest_pos_m$value/1000) ~ out_nbtest_pos_m$year, lwd = 1, col = 'orangered2')
    
    if (length(lik_f) > 0 | length(lik_m) > 0) {
      if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
      points(I(lik_f$totpos/1000) ~ I(lik_f$year + 0.5), pch = pchpt) 
      points(I(lik_m$totpos/1000) ~ I(lik_m$year + 0.5), pch = pchpt + 1) 
    }
  } else {
    plot(I(out_nbtest_pos_f$value/1000) ~ out_nbtest_pos_f$year, pch = '',
         ylim = c(0, max_ylim/1000 * 1.25),
         cex = 1, ylab = expression(paste(N^o, " of Positive Tests (in 1,000)")), 
         xlab = 'Year', xlim = c(2000, yr_pred), 
         main = expression(bold(paste("Total ", N^o, " of Positive Tests"))))
    lines(I(out_nbtest_pos_f$value/1000) ~ out_nbtest_pos_f$year, lwd = 1, col = 'orangered2')
    lines(I(out_nbtest_pos_m$value/1000) ~ out_nbtest_pos_m$year, lwd = 1, col = 'orangered2')
    if (!is.null(likdat$hts_pos$tot)) {
      if (cnt %in% redact) { pchpt <- '' } else { pchpt <- 16 }
      points(I(lik_f$totpos/1000) ~ I(lik_f$year + 0.5), pch = pchpt) 
      points(I(lik_m$totpos/1000) ~ I(lik_m$year + 0.5), pch = pchpt + 1) 
    }
  }
}

#' @export
plot_out_evertestneg <- function(mod, fp, likdat, cnt, survey_hts, out_evertest, simul = NULL, plot_legend = TRUE, yr_pred = 2018) {
    
    out_evertest <- subset(out_evertest, year <= yr_pred)
    out_evertest$year <- out_evertest$year + 0.5
    survey_hts$year <- survey_hts$year + 0.5
    
    # Decide if we plot CI or not
    if (!is.null(simul)){
      ci <- getCI(simul$ever.test)
      ci <- subset(ci, year <= yr_pred)
      ci[is.na(ci)] <- 0
      ci$year <- ci$year + 0.5
      plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    # Plots of Ever Test by Sex (among negatives)
    out_f <- subset(out_evertest, agegr == '15-49' & outcome == 'evertest' & sex =='female' & hivstatus == 'negative')
    out_m <- subset(out_evertest, agegr == '15-49' & outcome == 'evertest' & sex =='male' & hivstatus == 'negative')
    dat_f <- subset(survey_hts, country == cnt & outcome == 'evertest' & agegr == '15-49' & sex == 'female' & hivstatus == 'negative')
    dat_m <- subset(survey_hts, country == cnt & outcome == 'evertest' & agegr == '15-49' & sex == 'male' & hivstatus == 'negative')
    
    if (plot.ci == T){
        ci_f <- subset(ci, agegr == '15-49' & outcome == 'evertest' & sex =='female' & hivstatus == 'negative')
        ci_m <- subset(ci, agegr == '15-49' & outcome == 'evertest' & sex =='male' & hivstatus == 'negative')

        plot(I(out_f$value*100) ~ out_f$year, pch='', ylim=c(0,100), col='maroon', main="Negative Ever Tested",
          xlim=c(2000, yr_pred), ylab='Proportion Ever Tested Among Susceptibles', xlab='Year', lwd=1)
        polygon(x = c(ci_f$year, rev(ci_f$year)),
          y = c(I(ci_f$upper*100), rev(I(ci_f$lower*100))),
          col = col_ci, border = NA)
        polygon(x = c(ci_m$year, rev(ci_m$year)),
          y = c(I(ci_m$upper*100), rev(I(ci_m$lower*100))),
          col = col_ci, border = NA)
        lines(I(out_f$value*100) ~ out_f$year, col = 'maroon', lwd = 1)
        lines(I(out_m$value*100) ~ out_m$year, col = 'navy', lwd = 1)
        points(I(dat_f$est*100) ~ dat_f$year, pch = 15, col = 'maroon')
        points(I(dat_m$est*100) ~ dat_m$year, pch = 16, col = 'navy')
        segments(x0 = dat_f$year, y0 = I(dat_f$ci_l*100), 
                 x1 = dat_f$year, y1 = I(dat_f$ci_u*100), lwd = 1, col = 'maroon')
        segments(x0 = dat_m$year, y0 = I(dat_m$ci_l*100), 
                 x1 = dat_m$year, y1 = I(dat_m$ci_u*100), lwd = 1, col = 'navy')
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
                 legend = c('Women (15-49 years)','Men (15-49 years)'), 
                 col = c('maroon','navy'), bty = 'n', lwd = 1, pch = c(15,16)) }
    } else {
        plot(I(out_f$value*100) ~ out_f$year, type='l', ylim=c(0,100), 
             col = 'maroon', main = "Negative Ever Tested", xlab='Year', lwd=1,
             xlim = c(2000, yr_pred), ylab = 'Proportion Ever Tested Among Susceptibles')
        lines(I(out_m$value*100) ~ out_m$year, col = 'navy', lwd = 1)
        points(I(dat_f$est*100) ~ dat_f$year, pch = 15, col = 'maroon')
        points(I(dat_m$est*100) ~ dat_m$year, pch = 16, col = 'navy')
        segments(x0 = dat_f$year, y0 = I(dat_f$ci_l*100), 
                 x1 = dat_f$year, y1 = I(dat_f$ci_u*100), lwd = 1, col = 'maroon')
        segments(x0 = dat_m$year, y0 = I(dat_m$ci_l*100), 
                 x1 = dat_m$year, y1 = I(dat_m$ci_u*100), lwd = 1, col = 'navy')
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
                 legend = c('Women (15-49 years)','Men (15-49 years)'), 
                 col = c('maroon','navy'), bty = 'n', lwd = 2, pch = c(15,16)) }
    }
}

#' @export
plot_out_evertestpos <- function(mod, fp, likdat, cnt, survey_hts, out_evertest, simul = NULL, plot_legend = TRUE, yr_pred = 2018) {

  out_evertest <- subset(out_evertest, year <= yr_pred)
  out_evertest$year <- out_evertest$year + 0.5
  survey_hts$year <- survey_hts$year + 0.5
  
    # Decide if we plot CI or not
    if (!is.null(simul)){
        ci <- getCI(simul$ever.test)
        ci <- subset(ci, year <= yr_pred)
        ci[is.na(ci)] <- 0
        ci$year <- ci$year + 0.5
        plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    # Plots of Ever Test by Sex (among PLHIV)
    out_f <- subset(out_evertest, agegr == '15-49' & outcome == 'evertest' & sex =='female' & hivstatus == 'positive')
    out_m <- subset(out_evertest, agegr == '15-49' & outcome == 'evertest' & sex =='male' & hivstatus == 'positive')
    dat_f <- subset(survey_hts, country == cnt & agegr == '15-49' & sex == 'female' & hivstatus == 'positive' & outcome == 'evertest')
    dat_m <- subset(survey_hts, country == cnt & agegr == '15-49' & sex == 'male' & hivstatus == 'positive' & outcome == 'evertest')
    
    if (plot.ci == T){
        ci_f <- subset(ci, agegr == '15-49' & outcome == 'evertest' & sex =='female' & hivstatus == 'positive')
        ci_m <- subset(ci, agegr == '15-49' & outcome == 'evertest' & sex =='male' & hivstatus == 'positive')

        plot(I(out_f$value*100) ~ out_f$year, pch='', ylim=c(0,100), col='maroon', main='PLHIV Ever Tested',
          xlim=c(2000, yr_pred), ylab='Proportion PLHIV Ever Tested', xlab='Year', lwd=1)
          polygon(x = c(ci_f$year, rev(ci_f$year)),
          y = c(I(ci_f$upper*100), rev(I(ci_f$lower*100))),
          col = col_ci, border = NA)
        polygon(x = c(ci_m$year, rev(ci_m$year)),
          y = c(I(ci_m$upper*100), rev(I(ci_m$lower*100))),
          col = col_ci, border = NA)
        lines(I(out_f$value*100) ~ out_f$year, col = 'maroon', lwd = 1)
        lines(I(out_m$value*100) ~ out_m$year, col = 'navy', lwd = 1)
        points(I(dat_f$est*100) ~ I(dat_f$year-0.1), pch = 15, col = 'maroon')
        points(I(dat_m$est*100) ~ I(dat_m$year+0.1), pch = 16, col = 'navy')
        segments(x0 = dat_f$year-0.1, y0 = I(dat_f$ci_l*100), 
                 x1 = dat_f$year-0.1, y1 = I(dat_f$ci_u*100), lwd = 1, col = 'maroon')
        segments(x0 = dat_m$year+0.1, y0 = I(dat_m$ci_l*100), 
                 x1 = dat_m$year+0.1, y1 = I(dat_m$ci_u*100), lwd = 1, col = 'navy')
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
                 legend = c('Women (15-49 years)','Men (15-49 years)'), 
                 col = c('maroon','navy'), bty = 'n', lwd = 1, pch = c(15,16)) }
    } else {
        plot(I(out_f$value*100) ~ out_f$year, type='l', ylim = c(0,100), col = 'maroon', main = 'PLHIV Ever Tested',
          xlim = c(2000, yr_pred), ylab = 'Proportion PLHIV Ever Tested', xlab = 'Year', lwd = 1)
        lines(I(out_m$value*100) ~ out_m$year, col = 'navy', lwd = 1)
        points(I(dat_f$est*100) ~ I(dat_f$year-0.1), pch = 15, col = 'maroon')
        points(I(dat_m$est*100) ~ I(dat_m$year+0.1), pch = 16, col = 'navy')
        segments(x0 = dat_f$year-0.1, y0 = I(dat_f$ci_l*100), 
                 x1 = dat_f$year-0.1, y1 = I(dat_f$ci_u*100), lwd = 1, col = 'maroon')
        segments(x0 = dat_m$year+0.1, y0 = I(dat_m$ci_l*100), 
                 x1 = dat_m$year+0.1, y1 = I(dat_m$ci_u*100), lwd = 1, col = 'navy')
        if (plot_legend) {  
          legend(x = 2000, y = 100,
                 legend = c('Women (15-49 years)','Men (15-49 years)'), 
                 col = c('maroon','navy'), bty = 'n', lwd = 2, pch = c(15,16)) }
    }
}

#' @export
plot_out_evertest <- function(mod, fp, likdat, cnt, survey_hts, out_evertest, simul = NULL, plot_legend= TRUE, yr_pred = 2018) {

  out_evertest <- subset(out_evertest, year <= yr_pred)
  out_evertest$year <- out_evertest$year + 0.5
  survey_hts$year <- survey_hts$year + 0.5
  
    # Decide if we plot CI or not
    if (!is.null(simul)){
        ci <- getCI(simul$ever.test)
        ci <- subset(ci, year <= yr_pred)
        ci[is.na(ci)] <- 0
        ci$year <- ci$year + 0.5
        plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    # Plots of Ever Test by Sex (overall)
    out_f <- subset(out_evertest, agegr == '15-49' & outcome == 'evertest' & sex =='female' & hivstatus == 'all')
    out_m <- subset(out_evertest, agegr == '15-49' & outcome == 'evertest' & sex =='male' & hivstatus == 'all')
    dat_f <- subset(survey_hts, country == cnt & agegr == '15-49' & sex == 'female' & hivstatus == 'all' & outcome == 'evertest')
    dat_m <- subset(survey_hts, country == cnt & agegr == '15-49' & sex == 'male' & hivstatus == 'all' & outcome == 'evertest')
    
    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    if (plot.ci == T){
        ci_f <- subset(ci, agegr == '15-49' & outcome == 'evertest' & sex =='female' & hivstatus == 'all')
        ci_m <- subset(ci, agegr == '15-49' & outcome == 'evertest' & sex =='male' & hivstatus == 'all')

        plot(I(out_f$value*100) ~ out_f$year, pch='', ylim=c(0,100), col='maroon', main='Population Ever Tested',
          xlim=c(2000, yr_pred), ylab='Proportion Ever Tested', xlab='Year', lwd=1)
        polygon(x = c(ci_f$year, rev(ci_f$year)),
          y = c(I(ci_f$upper*100), rev(I(ci_f$lower*100))),
          col = col_ci, border = NA)
        lines(I(out_f$value*100) ~ out_m$year, col='maroon', lwd=1)
        polygon(x = c(ci_m$year, rev(ci_m$year)),
          y = c(I(ci_m$upper*100), rev(I(ci_m$lower*100))),
          col = col_ci, border = NA)
        lines(I(out_m$value*100) ~ out_m$year, col = 'navy', lwd = 1)
        points(I(dat_f$est*100) ~ I(dat_f$year + 0.1), pch = 15, col='maroon')
        points(I(dat_m$est*100) ~ I(dat_m$year - 0.1), pch = 16, col = 'navy')
        segments(x0 = dat_f$year + 0.1, y0 = I(dat_f$ci_l*100), 
                 x1 = dat_f$year + 0.1, y1 = I(dat_f$ci_u*100), lwd = 1, col = 'maroon')
        segments(x0 = dat_m$year - 0.1, y0 = I(dat_m$ci_l*100), 
                 x1 = dat_m$year - 0.1, y1 = I(dat_m$ci_u*100), lwd = 1, col = 'navy')
        if (plot_legend) { 
              legend(x = 2000, y = 100, 
                     legend = c('Women (15-49 years)','Men (15-49 years)'), 
                     col = c('maroon','navy'), bty = 'n', lwd = 1, pch = c(15,16)) }
    } else {
        plot(I(out_f$value*100) ~ out_f$year, type = 'l', ylim = c(0,100), 
          col = 'maroon', main = 'Population Ever Tested',
          xlim = c(2000, yr_pred), ylab = 'Proportion Ever Tested', 
          xlab = 'Year', lwd = 1)
        lines(I(out_m$value*100) ~ out_m$year, col = 'navy', lwd = 1)
        points(I(dat_f$est*100) ~ I(dat_f$year + 0.1), pch = 15, col = 'maroon')
        points(I(dat_m$est*100) ~ I(dat_m$year - 0.1), pch = 16, col = 'navy')
        segments(x0 = dat_f$year + 0.1, y0 = I(dat_f$ci_l*100), 
                 x1 = dat_f$year + 0.1, y1 = I(dat_f$ci_u*100), lwd = 1, col = 'maroon')
        segments(x0 = dat_m$year - 0.1, y0 = I(dat_m$ci_l*100), 
                 x1 = dat_m$year - 0.1, y1 = I(dat_m$ci_u*100), lwd = 1, col = 'navy')
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
                 legend = c('Women (15-49 years)','Men (15-49 years)'), 
                 col = c('maroon','navy'), bty = 'n', lwd = 2, pch = c(15,16)) }
    }
}

#' @export
plot_out_90s <- function(mod, fp, likdat, cnt, out_evertest, survey_hts, simul = NULL, plot_legend = TRUE, yr_pred = 2018) {

    phia_aware <- subset(survey_hts, country == cnt & agegr == '15-49' &
                  sex == 'both' & outcome == 'aware')
    phia_art <- subset(survey_hts, country == cnt & agegr == '15-49' &
                  sex == 'both' & outcome == 'art_both')
    phia_aware$year <- phia_aware$year + 0.5
    phia_art$year <- phia_art$year + 0.5
    
    if (all(is.na(phia_art$est))) {
      phia_art <- subset(survey_hts, country == cnt & agegr == '15-49' &
                  sex == 'both' & outcome == 'art_self')    
      phia_art$year <- phia_art$year + 0.5
    }
    svy_evertest <- subset(survey_hts, country == cnt & agegr == '15-49' & hivstatus == 'positive' &
                              sex == 'both' & outcome == 'evertest')
    svy_evertest$year <- svy_evertest$year + 0.5
      
    out_evertest <- subset(out_evertest, year <= yr_pred)
    out_aware <- subset(get_out_aware(mod, fp), year <= yr_pred)
    out_art <- subset(get_out_art(mod, fp), year >= 2000 & year <= yr_pred)
    out_evertest$year <- out_evertest$year + 0.5
    out_aware$year <- out_aware$year + 0.5
    out_art$year <- out_art$year + 0.5

    # Decide if we plot CI or not
    if (!is.null(simul)){
        ci <- lapply(simul[names(simul) != 'param'], getCI)
        ci <- lapply(ci, function(x) subset(x, year <= yr_pred))
        ci$ever.test[is.na(ci$ever.test)] <- 0
        ci$diagnoses[is.na(ci$diagnoses)] <- 0
        ci$diagnoses$year <- ci$diagnoses$year + 0.5
        ci$ever.test$year <- ci$ever.test$year + 0.5
        ci$number.test$year <- ci$number.test$year + 0.5
        ci$number.test.pos$year <- ci$number.test.pos$year + 0.5
        plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    # First 90
    out_f <- subset(out_aware, agegr == '15-49' & outcome == 'aware' & sex =='female')
    out_m <- subset(out_aware, agegr == '15-49' & outcome == 'aware' & sex =='male')
    out_a <- subset(out_aware, agegr == '15-49' & outcome == 'aware' & sex =='both')
    out_ever_all <- subset(out_evertest, agegr == '15-49' & outcome == 'evertest' &
        sex == 'both' & hivstatus == 'positive')
    
    if (plot.ci == TRUE){
        ci_f <- subset(ci$diagnoses, agegr == '15-49' & outcome == 'aware' & sex =='female')
        ci_m <- subset(ci$diagnoses, agegr == '15-49' & outcome == 'aware' & sex =='male')
        ci_a <- subset(ci$diagnoses, agegr == '15-49' & outcome == 'aware' & sex =='both')
        ci_ever_all <- subset(ci$ever.test, agegr == '15-49' & outcome == 'evertest' &
            sex == 'both' & hivstatus == 'positive')

        plot(I(out_a$value*100) ~ out_a$year, pch = '', ylim = c(0,100), 
             col = 'darkslategrey', main = 'PLHIV Ever Tested, 1st 90, and ART Coverage',
             xlim = c(2000, yr_pred), ylab = 'Proportion', xlab='Year', lwd = 1)
        polygon(x = c(ci_a$year, rev(ci_a$year)),
          y = c(I(ci_a$upper*100), rev(I(ci_a$lower*100))),
          col = col_ci, border = NA)
        lines(I(out_a$value*100) ~ out_a$year, col = 'darkslategrey', lwd = 1)
        lines(I(out_art$value*100) ~ out_art$year, col = 'deepskyblue2', lwd = 1, lty = 1)
        polygon(x = c(ci_ever_all$year, rev(ci_ever_all$year)),
          y = c(I(ci_ever_all$upper*100), rev(I(ci_ever_all$lower*100))),
          col = col_ci,border = NA)
        lines(I(out_ever_all$value*100) ~ out_ever_all$year, col = 'orange3', lwd = 1, lty = 1)
        legend(x = 2000, y = 100, legend = c('PLHIV Ever Tested','PLHIV Aware','ART Coverage'),
          col = c('orange3','darkslategrey','deepskyblue2'), pch = c(25, 1, 2), lwd = 1, bty = 'n')
        # We plot ever tested
        if (dim(svy_evertest)[1] > 0) {
          points(I(svy_evertest$est*100) ~ svy_evertest$year, pch = 25, bg = 'grey15') 
          segments(x0 = svy_evertest$year, y0 = I(svy_evertest$ci_l*100), 
                   x1 = svy_evertest$year, y1 = I(svy_evertest$ci_u*100), lwd = 1, col = 'bisque4') }
        # We plot PHIA if any
        if (dim(phia_aware)[1] > 0) {
          points(I(phia_aware$est*100) ~ phia_aware$year, pch = 1)
          segments(x0 = phia_aware$year, y0 = I(phia_aware$ci_l*100), 
                   x1 = phia_aware$year, y1 = I(phia_aware$ci_u*100), lwd = 1, col = 'bisque4')
          text('Empty symbols are \n direct survey measures \n (not included in fitting)', x = 2015, y = 12.5, cex = 0.8) 
          if (any(cnt == c('Swaziland','Zambia','Zimbabwe'))) {
            text('(not ARV-corrected)', x = 2015, y = 1, cex = 0.8) }
            }
          if (dim(phia_art)[1] > 0) {
            points(I(phia_art$est*100) ~ phia_art$year, pch = 2)
            segments(x0 = phia_art$year, y0 = I(phia_art$ci_l*100), 
                   x1 = phia_art$year, y1 = I(phia_art$ci_u*100), lwd = 1, col = 'bisque4') }
    } else {
        plot(I(out_a$value*100) ~ out_a$year, type = 'l', ylim = c(0,100), 
             col = 'darkslategrey', main='PLHIV Ever Tested, 1st 90, and ART Coverage',
             xlim = c(2000, yr_pred), ylab = 'Proportion', xlab = 'Year', lwd = 1)
        lines(I(out_art$value*100) ~ out_art$year, col = 'deepskyblue2', lwd = 1, lty = 1)
        lines(I(out_ever_all$value*100) ~ out_ever_all$year, col='orange3', lwd = 1, lty = 1)
        legend(x = 2000, y = 100, legend = c('PLHIV Ever Tested','PLHIV Aware','ART Coverage'),
               col = c('orange3','darkslategrey','deepskyblue2'), pch = c(25, 1, 2), lwd = 1, bty = 'n')
        # We plot ever tested
        if (dim(svy_evertest)[1] > 0) {
          points(I(svy_evertest$est*100) ~ svy_evertest$year, pch = 25, bg = 'grey15') 
          segments(x0 = svy_evertest$year, y0 = I(svy_evertest$ci_l*100), 
                   x1 = svy_evertest$year, y1 = I(svy_evertest$ci_u*100), lwd = 1, col = 'bisque4') }
        # We plot PHIA if any
        if (dim(phia_aware)[1]>0) {
          points(I(phia_aware$est*100) ~ phia_aware$year, pch = 1)
          segments(x0 = phia_aware$year, y0 = I(phia_aware$ci_l*100), 
                   x1 = phia_aware$year, y1 = I(phia_aware$ci_u*100), lwd = 1, col = 'bisque4')
          text('X-Validation', x = 2015, y = 12.5, cex = 0.8) 
          if (any(cnt == c('Swaziland','Zambia','Zimbabwe'))) {
            text('(not ARV-corrected)', x = 2015, y = 1, cex = 0.8) }
          }
        if (dim(phia_art)[1]>0) {
          points(I(phia_art$est*100) ~ phia_art$year, pch = 2)
          segments(x0 = phia_art$year, y0 = I(phia_art$ci_l*100), 
                   x1 = phia_art$year, y1 = I(phia_art$ci_u*100), lwd = 1, col = 'bisque4') }
    }
}

#' @export
plot_out_evertest_fbyage <- function(mod, fp, likdat, cnt, survey_hts, out_evertest, 
                            simul = NULL, plot_legend = TRUE, yr_pred = 2018) {

  out_evertest <- subset(out_evertest, year <= yr_pred)
  out_evertest$year <- out_evertest$year + 0.5
  survey_hts$year <- survey_hts$year + 0.5
  
    # Decide if we plot CI or not
    if (!is.null(simul)){
        ci <- getCI(simul$ever.test)
        ci <- subset(ci, year <= yr_pred)
        ci[is.na(ci)] <- 0
        ci$year <- ci$year + 0.5
        plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    # For females
    out_1524f <- subset(out_evertest, agegr == '15-24' & outcome == 'evertest' & hivstatus == 'all' & sex == 'female')
    out_2534f <- subset(out_evertest, agegr == '25-34' & outcome == 'evertest' & hivstatus == 'all' & sex == 'female')
    out_3549f <- subset(out_evertest, agegr == '35-49' & outcome == 'evertest' & hivstatus == 'all' & sex == 'female')
    datage15f <- subset(survey_hts, country == cnt & hivstatus == "all" &
        agegr == "15-24" & sex == "female" & outcome == "evertest")
    datage25f <- subset(survey_hts, country == cnt & hivstatus == "all" &
        agegr == "25-34" & sex == "female" & outcome == "evertest")
    datage35f <- subset(survey_hts, country == cnt & hivstatus == "all" &
        agegr == "35-49" & sex == "female" & outcome == "evertest")
    
    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    if (plot.ci == T){
        ci_1524f <- subset(ci, agegr == '15-24' & outcome == 'evertest' & hivstatus=='all' & sex=='female')
        ci_2534f <- subset(ci, agegr == '25-34' & outcome == 'evertest' & hivstatus=='all' & sex=='female')
        ci_3549f <- subset(ci, agegr == '35-49' & outcome == 'evertest' & hivstatus=='all' & sex=='female')
        plot(I(out_1524f$value*100) ~ out_1524f$year, pch = '', ylim = c(0,100), 
             xlim=c(2000, yr_pred),
             ylab = 'Proportion ever tested', xlab = 'Year', 
             main = 'Women Ever Tested, by Age')
        polygon(x = c(ci_1524f$year, rev(ci_1524f$year)),
          y = c(I(ci_1524f$upper*100), rev(I(ci_1524f$lower*100))),
          col = col_ci, border = NA)
        polygon(x = c(ci_2534f$year, rev(ci_2534f$year)),
          y = c(I(ci_2534f$upper*100), rev(I(ci_2534f$lower*100))),
          col = col_ci, border = NA)
        polygon(x = c(ci_3549f$year, rev(ci_3549f$year)),
          y = c(I(ci_3549f$upper*100), rev(I(ci_3549f$lower*100))),
          col = col_ci, border = NA)
        lines(I(out_1524f$value*100) ~ out_1524f$year, col = 'maroon', lwd = 1, lty = 3)
        lines(I(out_2534f$value*100) ~ out_2534f$year, col = 'maroon', lwd = 1, lty = 2)
        lines(I(out_3549f$value*100) ~ out_3549f$year, col = 'maroon', lwd = 1, lty = 1)
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
                 legend = c('15-24 years', '25-34 years','35-49 years'), 
                 col = c('maroon','maroon'),
                 lwd = 2, lty = c(3,2,1), pch = c(16,15,17), pt.cex = 0.7, bty = 'n') }
        points(I(datage15f$est*100) ~ I(datage15f$year-0.1), pch = 16, col = 'goldenrod3', cex = 0.8)
        points(I(datage25f$est*100) ~ datage25f$year, pch=15, col = 'goldenrod3', cex = 0.8)
        points(I(datage35f$est*100) ~ I(datage35f$year+0.1), pch = 17, col = 'goldenrod3', cex = 0.8)
        segments(x0 = datage15f$year-0.1, y0 = I(datage15f$ci_l*100), 
                 x1 = datage15f$year-0.1, y1 = I(datage15f$ci_u*100), col = 'goldenrod3', lwd=1)
        segments(x0 = datage25f$year, y0 = I(datage25f$ci_l*100), 
                 x1 = datage25f$year, y1 = I(datage25f$ci_u*100), col = 'goldenrod3', lwd = 1)
        segments(x0 = datage35f$year+0.1, y0 = I(datage35f$ci_l*100), 
                 x1 = datage35f$year+0.1, y1 = I(datage35f$ci_u*100), col = 'goldenrod3', lwd = 1)
    } else {
        plot(I(out_1524f$value*100) ~ out_1524f$year, type = 'l', col = 'maroon',
             ylim = c(0,100), xlim = c(2000, yr_pred), lwd = 1, lty = 3,
             ylab = 'Proportion ever tested', xlab = 'Year', 
             main = 'Women Ever Tested, by Age')
        lines(I(out_2534f$value*100) ~ out_2534f$year, col = 'maroon', lwd = 1, lty = 2)
        lines(I(out_3549f$value*100) ~ out_3549f$year, col = 'maroon', lwd = 1, lty = 1)
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
               legend = c('15-24 years', '25-34 years','35-49 years'), 
               col = c('maroon','maroon'),
               lwd = 2, lty = c(3,2,1), pch = c(16,15,17), pt.cex = 0.7, bty = 'n') }
        points(I(datage15f$est*100) ~ I(datage15f$year-0.1), pch = 16, col = 'goldenrod3', cex = 0.8)
        points(I(datage25f$est*100) ~ datage25f$year, pch = 15, col = 'goldenrod3', cex = 0.8)
        points(I(datage35f$est*100) ~ I(datage35f$year+0.1), pch = 17, col = 'goldenrod3', cex = 0.8)
        segments(x0 = datage15f$year-0.1, y0 = I(datage15f$ci_l*100), 
                 x1 = datage15f$yea-0.1, y1 = I(datage15f$ci_u*100), col = 'goldenrod3', lwd = 1)
        segments(x0 = datage25f$year, y0 = I(datage25f$ci_l*100), 
                 x1 = datage25f$year, y1 = I(datage25f$ci_u*100), col = 'goldenrod3', lwd = 1)
        segments(x0 = datage35f$year+0.1, y0 = I(datage35f$ci_l*100), 
                 x1 = datage35f$year+0.1, y1 = I(datage35f$ci_u*100), col = 'goldenrod3', lwd = 1)
    }
}

#' @export
plot_out_evertest_mbyage <- function(mod, fp, likdat, cnt, survey_hts, out_evertest, simul = NULL, plot_legend = TRUE, yr_pred = 2018) {

  out_evertest <- subset(out_evertest, year <= yr_pred)
  out_evertest$year <- out_evertest$year + 0.5
  survey_hts$year <- survey_hts$year + 0.5
  
    # Decide if we plot CI or not
    if (!is.null(simul)){
      ci <- getCI(simul$ever.test)
      ci <- subset(ci, year <= yr_pred)
      ci[is.na(ci)] <- 0
      ci$year <- ci$year + 0.5
      plot.ci <- TRUE
    }
    if (is.null(simul)) plot.ci <- FALSE

    out_1524m <- subset(out_evertest, agegr == '15-24' & outcome == 'evertest' & hivstatus=='all' & sex=='male')
    out_2534m <- subset(out_evertest, agegr == '25-34' & outcome == 'evertest' & hivstatus=='all' & sex=='male')
    out_3549m <- subset(out_evertest, agegr == '35-49' & outcome == 'evertest' & hivstatus=='all' & sex=='male')
    datage15m <- subset(survey_hts, country == cnt & hivstatus == "all" &
        agegr == "15-24" & sex == "male" & outcome == "evertest")
    datage25m <- subset(survey_hts, country == cnt & hivstatus == "all" &
        agegr == "25-34" & sex == "male" & outcome == "evertest")
    datage35m <- subset(survey_hts, country == cnt & hivstatus == "all" &
        agegr == "35-49" & sex == "male" & outcome == "evertest")

    col_ci <- rgb(150, 150, 150, 125, max = 255)
    
    if (plot.ci == T){
        ci_1524m <- subset(ci, agegr == '15-24' & outcome == 'evertest' & hivstatus=='all' & sex=='male')
        ci_2534m <- subset(ci, agegr == '25-34' & outcome == 'evertest' & hivstatus=='all' & sex=='male')
        ci_3549m <- subset(ci, agegr == '35-49' & outcome == 'evertest' & hivstatus=='all' & sex=='male')

        plot(I(out_1524m$value*100) ~ out_1524m$year, pch = '', 
          ylim = c(0,100), xlim = c(2000, yr_pred),
          ylab = 'Proportion ever tested', xlab = 'Year',
          main = 'Men Ever Tested, by Age')
        polygon(x = c(ci_1524m$year, rev(ci_1524m$year)),
          y = c(I(ci_1524m$upper*100), rev(I(ci_1524m$lower*100))),
          col = col_ci, border = NA)
        polygon(x = c(ci_2534m$year, rev(ci_2534m$year)),
          y = c(I(ci_2534m$upper*100), rev(I(ci_2534m$lower*100))),
          col = col_ci, border = NA)
        polygon(x = c(ci_3549m$year, rev(ci_3549m$year)),
          y = c(I(ci_3549m$upper*100), rev(I(ci_3549m$lower*100))),
          col = col_ci, border = NA)
        lines(I(out_1524m$value*100) ~ I(out_1524m$year-0.1), col = 'navy', lwd = 1, lty = 3)
        lines(I(out_2534m$value*100) ~ out_2534m$year, col = 'navy', lwd = 1, lty = 2)
        lines(I(out_3549m$value*100) ~ I(out_3549m$year+0.1), col = 'navy', lwd = 1, lty = 1)
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
                 legend = c('15-24', '25-34', '35-49'), col = c('navy','navy'),
                  lwd = 2, lty = c(3,2,1), pch = c(16,15,17), pt.cex = 0.7, bty = 'n') }
        points(I(datage15m$est*100) ~ I(datage15m$year-0.1), pch = 16, col = 'goldenrod3', cex = 0.8)
        points(I(datage25m$est*100) ~ datage25m$year, pch = 15, col = 'goldenrod3', cex = 0.8)
        points(I(datage35m$est*100) ~ I(datage25m$year+0.1), pch = 17, col = 'goldenrod3', cex = 0.8)
        segments(x0 = datage15m$year-0.1, y0 = I(datage15m$ci_l*100), 
                 x1 = datage15m$year-0.1, y1 = I(datage15m$ci_u*100), col = 'goldenrod3', lwd = 1)
        segments(x0 = datage25m$year, y0 = I(datage25m$ci_l*100), 
                 x1=datage25m$year, y1 = I(datage25m$ci_u*100), col = 'goldenrod3', lwd = 1)
        segments(x0 = datage35m$year+0.1, y0 = I(datage35m$ci_l*100), 
                 x1 = datage35m$year+0.1, y1 = I(datage35m$ci_u*100), col = 'goldenrod3', lwd = 1)
    } else {
        plot(I(out_1524m$value*100) ~ out_1524m$year, type = 'l', col = 'navy', 
             ylim = c(0,100), xlim = c(2000, yr_pred), lwd = 1, lty = 3,
             ylab = 'Proportion ever tested', xlab = 'Year', 
             main = 'Men Ever Tested, by Age')
        lines(I(out_2534m$value*100) ~ out_2534m$year, col = 'navy', lwd = 1, lty = 2)
        lines(I(out_3549m$value*100) ~ out_3549m$year, col = 'navy', lwd = 1, lty = 1)
        if (plot_legend) { 
          legend(x = 2000, y = 100, 
                 legend = c('15-24', '25-34', '35-49'), col = c('navy','navy'),
                 lwd = 2, lty = c(3,2,1), pch = c(16,15,17), pt.cex = 0.7, bty = 'n') }
        points(I(datage15m$est*100) ~ I(datage15m$year-0.1), pch = 16, col = 'goldenrod3', cex = 0.8)
        points(I(datage25m$est*100) ~ datage25m$year, pch = 15, col = 'goldenrod3', cex = 0.8)
        points(I(datage35m$est*100) ~ I(datage25m$year+0.1), pch = 17, col = 'goldenrod3', cex = 0.8)
        segments(x0=datage15m$year-0.1, y0=I(datage15m$ci_l*100), 
                 x1=datage15m$year-0.1, y1=I(datage15m$ci_u*100), col = 'goldenrod3', lwd = 1)
        segments(x0=datage25m$year, y0=I(datage25m$ci_l*100), 
                 x1=datage25m$year, y1=I(datage25m$ci_u*100), col = 'goldenrod3', lwd = 1)
        segments(x0=datage35m$year+0.1, y0=I(datage35m$ci_l*100), 
                 x1=datage35m$year+0.1, y1=I(datage25m$ci_u*100), col = 'goldenrod3', lwd = 1)
    }
}


# ---- Single Function for Plots ----
plot_out <- function(mod, fp, likdat, cnt, survey_hts, out_evertest, simul = NULL, plot_legend = TRUE, yr_pred = 2018) {
par(mfrow = c(3,2), mar = c(4,4,2,2))
  plot_out_nbtest(mod, fp, likdat, cnt, simul, yr_pred)
  plot_out_nbpostest(mod, fp, likdat, cnt, simul, yr_pred)  
  plot_out_evertestneg(mod, fp, likdat, cnt, survey_hts, out_evertest, simul, plot_legend, yr_pred)
  plot_out_evertestpos(mod, fp, likdat, cnt, survey_hts, out_evertest, simul, plot_legend, yr_pred)
  plot_out_evertest(mod, fp, likdat, cnt, survey_hts, out_evertest, simul, plot_legend, yr_pred)
  plot_out_90s(mod, fp, likdat, cnt, out_evertest, survey_hts, simul, plot_legend, yr_pred) 
}

#' @export
plot_out_strat <- function(mod, fp, likdat, cnt, survey_hts, out_evertest, simul = NULL, plot_legend = TRUE, yr_pred = 2018) {
  par(mfrow = c(1,2), mar = c(4,4,2,2))
  plot_out_evertest_mbyage(mod, fp, likdat, cnt, survey_hts, out_evertest, simul, plot_legend, yr_pred) 
  plot_out_evertest_fbyage(mod, fp, likdat, cnt, survey_hts, out_evertest, simul, plot_legend, yr_pred)    
}

# ---- Output Tabular Format ----
#' @export
end_of_year <- function(year, value){
  if (length(unique(year)) != length(year)) { print('non unique years'); break }
  new_x <- year + 0.5
  new_value <- approx(year, value, new_x, method = 'linear', rule = 2)$y
  return(new_value)
} 


#' @export
tab_out_evertest <- function(mod, fp, age_grp = '15-49', gender = 'both', hiv = 'all', year_range = c(2010, 2018), simul = NULL, end_year = TRUE) {
  if (length(year_range) == 1) { year_range <- c(year_range, year_range) }
  if (is.null(simul)) {
    out <- get_out_evertest(mod, fp, age_grp, gender, hiv)
    if (end_year == TRUE) { out$value <- end_of_year(out$year, out$value) }
    out$value <- round(out$value * 100, 1)
    tab_evertest <- subset(out, year >= year_range[1] & year <= year_range[2])
  } else {
    out <- get_out_evertest(mod, fp, age_grp, gender, hiv)
    if (end_year == TRUE) { out$value <- end_of_year(out$year, out$value) }
    out$value <- round(out$value * 100, 1)
    outci <- getCI(simul$ever.test)
    outci <- subset(outci, agegr == age_grp & sex == gender & hivstatus == hiv)
    if (end_year == TRUE) {
      outci$lower <- end_of_year(outci$year, outci$lower)
      outci$upper <- end_of_year(outci$year, outci$upper) }
    outci$lower <- round(outci$lower * 100, 1)
    outci$upper <- round(outci$upper * 100, 1)
    outci <- subset(outci, year >= year_range[1] & year <= year_range[2]) 
    tab_evertest <- merge(out, outci)
  }
  row.names(tab_evertest) <- NULL
  tab_evertest
}

#' @export
tab_out_aware <- function(mod, fp, age_grp = '15-49', gender = 'both', year_range = c(2010, 2018), simul = NULL, end_year = TRUE) {
  if (length(year_range) == 1) {
    year_range <- c(year_range, year_range)
  }

  if (is.null(simul)) {
    out <- get_out_aware(mod, fp, age_grp, gender)

    if (end_year == TRUE) {
      out$value <- end_of_year(out$year, out$value)
    }
    out$value <- round(out$value * 100, 1)
    tab_aware <- subset(out, year >= year_range[1] & year <= year_range[2])

  } else {

    out <- get_out_aware(mod, fp, age_grp, gender)
    if (end_year == TRUE) {
      out$value <- end_of_year(out$year, out$value)
    }
    
    out$value <- round(out$value * 100, 1)
    outci <- getCI(simul$diagnoses)
    outci <- subset(outci, agegr == age_grp & sex == gender) 

    if (end_year == TRUE) { 
      outci$lower <- end_of_year(outci$year, outci$lower)
      outci$upper <- end_of_year(outci$year, outci$upper)
    }

    outci$lower <- round(outci$lower * 100, 1)
    outci$upper <- round(outci$upper * 100, 1)
    outci <- subset(outci, year >= year_range[1] & year <= year_range[2]) 

    tab_aware <- merge(out, outci)
  }
  
  row.names(tab_aware) <- NULL
  tab_aware
}

#' @export
tab_out_nbaware <- function(mod, fp, age_grp = '15-49', gender = 'both', year_range = c(2010, 2018), end_year = TRUE) {
  if (length(year_range) == 1) { year_range <- c(year_range, year_range) }
    out <- get_out_nbaware(mod, fp, age_grp, gender)
    if (end_year == TRUE) { out$value <- end_of_year(out$year, out$value) }
    out$value <- round(out$value, 0)
    tab_nbaware <- subset(out, year >= year_range[1] & year <= year_range[2])

  row.names(tab_nbaware) <- NULL
  tab_nbaware
}


#' @export
tab_out_artcov <- function(mod, fp, gender = 'both', year_range = c(2010, 2018)) {
  ## ART coverage is already end-of-year, no need to adjust
  
  if (length(year_range) == 1) {
    year_range <- c(year_range, year_range)
  }

  artcov_m <- data.frame(year = fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L,
                         outcome = "artcov",
                         agegr = "15+",
                         sex = "male",
                         hivstatus = "positive",
                         value = fp$art15plus_num[1,] / colSums(mod[,1,2,]))

  artcov_f <- data.frame(year = fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L,
                         outcome = "artcov",
                         agegr = "15+",
                         sex = "female",
                         hivstatus = "positive",
                         value = fp$art15plus_num[2,] / colSums(mod[,2,2,]))

  artcov_b <- data.frame(year = fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L,
                         outcome = "artcov",
                         agegr = "15+",
                         sex = "both",
                         hivstatus = "positive",
                         value = colSums(fp$art15plus_num) / colSums(mod[,,2,],,2))

  out <- rbind(artcov_b, artcov_m, artcov_f)
  out$value <- round(out$value * 100, 1)
  tab_artcov <- subset(out, year >= year_range[1] & year <= year_range[2] &
                            sex %in% gender)
  row.names(tab_artcov) <- NULL
  tab_artcov
}

#' @export
tab_out_pregprev <- function(mod, fp, year_range = c(2010, 2018), end_year = TRUE) {
  if (length(year_range) == 1) { year_range <- c(year_range, year_range) }
    out <- get_out_pregprev(mod, fp)
  out$prev <- round(out$prev * 100, 1)
  out$aware <- round(out$aware * 100, 1)
  out$artcov <- round(out$artcov * 100, 1)
  tab_pregprev <- subset(out, year >= year_range[1] & year <= year_range[2])
  row.names(tab_pregprev) <- NULL
  tab_pregprev
}
