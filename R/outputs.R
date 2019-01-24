#' Calculate number of tests conducted or number tested by year
#' 
#' @param mod model output of class 'eppasm'
#' @param fp parameter inputs (class 'specfp')
#' @param df a data.frame with indices for prediction. See [evertest()] for more information.
#'
#' @return a data.frame consisting of the number of tests, number tested in the past 12 months, and population size corresponding to rows of df
#'
#' @details
#' 
#' Number of tests (or number tested in past 12 months) are approximated by mid-year
#' counts and annual testing rates within each stratum.
#'
#' @useDynLib first90 number_testsC
#' @export
number_tests <- function(mod, fp, df, tested12m = FALSE, VERSION = "C") {

  if(tested12m){
    hts_prop <- 1 - exp(-fp$hts_rate)
    diagn_prop <- 1 - exp(-fp$diagn_rate)
  }
  
  if(VERSION != "R"){
    
    if(tested12m)
      stop("tested12m option is not implemented in C version")
    
    val <- .Call(number_testsC,
                 mod,
                 attr(mod, "testnegpop"),
                 attr(mod, "hivpop"),
                 attr(mod, "diagnpop"),
                 attr(mod, "artpop"),
                 fp$hts_rate,
                 fp$diagn_rate,
                 as.integer(df$haidx),
                 as.integer(df$hagspan),
                 as.integer(df$sidx),
                 as.integer(df$hvidx),
                 as.integer(df$yidx),
                 as.integer(fp$ss$agfirst.idx),
                 as.integer(fp$ss$h.ag.span),
                 # Mathieu added
                 attr(mod, "late_diagnoses"))
    
    names(val) <- c("pop", "tests")
    val <- as.data.frame(val)
    
  } else {

    tests <- pop <- numeric(length(df$haidx))

    if(tested12m)
      ntested12m <- numeric(length(df$haidx))

    for(i in seq_along(df$haidx)) {

      haidx <- df$haidx[i] + 1:df$hagspan[i] - 1
      sidx <- if(df$sidx[i] == 0) 1:2 else df$sidx[i]

      paidx <- fp$ss$agfirst.idx[df$haidx[i]] + 1:sum(fp$ss$h.ag.span[haidx]) - 1
            
      if(df$hvidx[i] %in% c(0, 1)){ # testing among HIV-

        pop_hivn_ha <- apply(mod[paidx, sidx, fp$ss$hivn.idx, df$yidx[i], drop = FALSE], 2:3, fastmatch::ctapply, fp$ss$ag.idx[paidx], sum)
        tested_ha <- attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivn.idx, df$yidx[i], drop = FALSE]

        tests[i] <- tests[i] +
          sum((c(pop_hivn_ha) - c(tested_ha)) * fp$hts_rate[haidx, sidx, 1, df$yidx[i]]) + # among untested HIV- population
          sum(tested_ha * fp$hts_rate[haidx, sidx, 2, df$yidx[i], drop = FALSE])  # among previously tested HIV- population

        if(tested12m)
          ntested12m[i] <- ntested12m[i] +
            sum((c(pop_hivn_ha) - c(tested_ha)) * hts_prop[haidx, sidx, 1, df$yidx[i]]) + # among untested HIV- population
            sum(tested_ha * hts_prop[haidx, sidx, 2, df$yidx[i], drop = FALSE])  # among previously tested HIV- population
        
        pop[i] <- pop[i] + sum(pop_hivn_ha)
        
      }
      
      if(df$hvidx[i] %in% c(0, 2)){ # testing among HIV+

        hivpop_ha_hm <- attr(mod, "hivpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i], drop = FALSE]
        diagn_ha_hm <- attr(mod, "diagnpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i], drop = FALSE]
        artpop_ha_hm <- colSums(attr(mod, "artpop")[ , , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i], drop = FALSE])

        
        ## Calculation proportion ever tested among undiagnosed
        undiagnosed_ha_hm <- hivpop_ha_hm - diagn_ha_hm
        testneg_ha <- attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivp.idx, df$yidx[i], drop = FALSE]
        prop_testneg <- 1 / colSums(hivpop_ha_hm - diagn_ha_hm) * c(testneg_ha)

        naive_ha_hm <- sweep(undiagnosed_ha_hm, 2:4, 1 - prop_testneg, "*")
        testneg_ha_hm <- sweep(undiagnosed_ha_hm, 2:4, prop_testneg, "*")

        tests[i] <- tests[i] +
          sum(naive_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 1, df$yidx[i]])) +
          sum(testneg_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 2, df$yidx[i]])) +
          sum(diagn_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 3, df$yidx[i]])) +
          sum(artpop_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 4, df$yidx[i]])) +
          sum(attr(mod, "late_diagnoses")[, haidx, sidx, df$yidx[i]])

        if(tested12m)
          ntested12m[i] <- ntested12m[i] +
            sum(naive_ha_hm * c(diagn_prop[ , haidx, sidx, 1, df$yidx[i]])) +
            sum(testneg_ha_hm * c(diagn_prop[ , haidx, sidx, 2, df$yidx[i]])) +
            sum(diagn_ha_hm * c(diagn_prop[ , haidx, sidx, 3, df$yidx[i]])) +
            sum(artpop_ha_hm * c(diagn_prop[ , haidx, sidx, 4, df$yidx[i]]))
        
        pop[i] <- pop[i] +
          sum(hivpop_ha_hm) +
          sum(artpop_ha_hm)
      }
    }
  
    val <- data.frame(pop = pop, tests = tests)
    if(tested12m)
      val$ntested12m <- ntested12m
  }

  val
}

#' @export
total_tests <- function(mod, df) {
  all_diag <- NA
  length(all_diag) <- length(df$haidx)
  for(i in seq_along(df$haidx)) {
    hivdx <- switch(df$hivstatus[i], 
                    all = 1:6,
                    negative = 1:2,
                    positive = 3:6)
    haidx <- df$haidx[i] + 1:df$hagspan[i] - 1
    sidx <- if(df$sidx[i] == 0) 1:2 else df$sidx[i]
    
    all_diag[i] <- colSums(attr(mod, "hivtests")[haidx, sidx, hivdx, df$yidx[i], drop = FALSE],,3)
    
    # We add late diag if hivstatus is all or positive  
    if (df$hivstatus[i] != 'negative') {
      late_diag <- colSums(attr(mod, "late_diagnoses")[, haidx, sidx, df$yidx[i], drop = FALSE],,3)
      all_diag[i] <- all_diag[i] + late_diag }
  }
  return(all_diag)
}


#' Proportion ever tested by age, sex, or HIV status
#'
#' This function calculates proportion of ever tested among a population
#' stratified by age group, sex, HIV status, and year.
#'
#' @param mod     simulation model output
#' @param fp      simulation model parameter inputs
#' @param df      a data.frame with indices for prediction, see Details.
#'
#' @return a vector
#'
#' @details
#'
#' Age groups are specified in terms of aggregate HIV age groups: 15-16, 17-19, 20-24, ..., 45-49, 50+.
#' Another function could be added to handle other age groups if needed, with additional computational
#' complexity.
#'
#' The argument df must contain the following columns:
#' 
#' * haidx:   HIV age group (1 = 15-16, 2 = 17-19, 3 = 20-24, ..., 8 = 45-49, 9 = 50+)
#' * sidx:    sex (1 = male, 2 = female, 0 = both)
#' * hvidx:   HIV status (1 = negative, 2 = positive, 0 = all)
#' * yidx:    year index
#' * hagspan: number of HIV age groups to span
#'
#' @examples
#'
#' \dontrun{
#' data(survey_hts)
#' dat <- subset(survey_hts, country == "Malawi" & outcome == "evertest")
#' df <- add_ss_indices(dat, fp$ss)
#'
#' df$pred <- evertest(mod, fp, df)
#' }
#' 
#' @export
#' @useDynLib first90 evertestC
evertest <- function(mod, fp, df, VERSION = "C") {

  if(VERSION != "R")

    val <- .Call(evertestC,
                 mod,
                 attr(mod, "testnegpop"),
                 attr(mod, "diagnpop"),
                 attr(mod, "artpop"),
                 as.integer(df$haidx),
                 as.integer(df$hagspan),
                 as.integer(df$sidx),
                 as.integer(df$hvidx),
                 as.integer(df$yidx),
                 as.integer(fp$ss$agfirst.idx),
                 as.integer(fp$ss$h.ag.span))

  else { 

    val <- numeric(length(df$haidx))

    for(i in seq_along(df$haidx)) {
      
      haidx <- df$haidx[i] + 1:df$hagspan[i] - 1
      sidx <- if(df$sidx[i] == 0) 1:2 else df$sidx[i]

      paidx <- fp$ss$agfirst.idx[df$haidx[i]] + 1:sum(fp$ss$h.ag.span[haidx]) - 1
  
      tested <- 0
      pop <- 0
      
      if(df$hvidx[i] %in% c(0, 1)){
        tested <- tested + sum(attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivn.idx, df$yidx[i]])
        pop <- pop + sum(mod[paidx, sidx, fp$ss$hivn.idx, df$yidx[i]])
      }
      
      if(df$hvidx[i] %in% c(0, 2)){
        tested <- tested +
          sum(attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivp.idx, df$yidx[i]]) +
          sum(attr(mod, "diagnpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]]) +
          sum(attr(mod, "artpop")[ , , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]])
        pop <- pop + sum(mod[paidx, sidx, fp$ss$hivp.idx, df$yidx[i]])
      }

      val[i] <- tested / pop
    }
  }

  val
}




#' Proportion diagnosed among HIV+ by age, sex, or HIV status
#'
#' This function calculates proportion diagnosed among the HIV positive population
#' stratified by age group, sex, HIV status, and year.
#'
#' @param mod     simulation model output
#' @param fp      simulation model parameter inputs
#' @param df      a data.frame with indices for prediction, see Details.
#'
#' @return a vector
#'
#' @details
#'
#' Age groups are specified in terms of aggregate HIV age groups: 15-16, 17-19, 20-24, ..., 45-49, 50+.
#' Another function could be added to handle other age groups if needed, with additional computational
#' complexity.
#'
#' The argument df must contain the following columns:
#' 
#' * haidx:   HIV age group (1 = 15-16, 2 = 17-19, 3 = 20-24, ..., 8 = 45-49, 9 = 50+)
#' * sidx:    sex (1 = male, 2 = female, 0 = both)
#' * yidx:    year index
#' * hagspan: number of HIV age groups to span
#'
#' @examples
#'
#' \dontrun{
#' data(survey_hts)
#' dat <- subset(survey_hts, country == "Malawi" & outcome == "aware")
#' df <- add_ss_indices(dat, fp$ss)
#'
#' df$pred <- diagnosed(mod, fp, df)
#' }
#' 
#' @export
#' @useDynLib first90 diagnosedC
diagnosed <- function(mod, fp, df, VERSION = "C") {

  if(VERSION != "R") {

    val <- .Call(diagnosedC,
             mod,
             attr(mod, "diagnpop"),
             attr(mod, "artpop"),
             as.integer(df$haidx),
             as.integer(df$hagspan),
             as.integer(df$sidx),
             as.integer(df$yidx),
             as.integer(fp$ss$agfirst.idx),
             as.integer(fp$ss$h.ag.span))

  } else { 

    val <- numeric(length(df$haidx))

    for(i in seq_along(df$haidx)) {
      
      haidx <- df$haidx[i] + 1:df$hagspan[i] - 1
      sidx <- if(df$sidx[i] == 0) 1:2 else df$sidx[i]

      paidx <- fp$ss$agfirst.idx[df$haidx[i]] + 1:sum(fp$ss$h.ag.span[haidx]) - 1

      diagnosed <- sum(attr(mod, "diagnpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]]) +
        sum(attr(mod, "artpop")[ , , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]])
      hivpop <- sum(mod[paidx, sidx, fp$ss$hivp.idx, df$yidx[i]])
      
      val[i] <- diagnosed / hivpop
    }
  }

  val
}

nb_aware <- function(mod, fp, df) {
  
    val <- numeric(length(df$haidx))
    
    for(i in seq_along(df$haidx)) {
      
      haidx <- df$haidx[i] + 1:df$hagspan[i] - 1
      sidx <- if(df$sidx[i] == 0) 1:2 else df$sidx[i]
      
      paidx <- fp$ss$agfirst.idx[df$haidx[i]] + 1:sum(fp$ss$h.ag.span[haidx]) - 1
      
      aware <- sum(attr(mod, "diagnpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]]) +
        sum(attr(mod, "artpop")[ , , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]])
      
      val[i] <- aware
    }
  val
}

#' Number of new diagnoses by age and sex
#'
#' @return
#'
#' A data frame reporting new diagnoses in two ways:
#'
#' 1. Approximated from the testing rates and mid-year populations
#'    as in the function [number_tests()].
#' 1. Based on the number of new diagnoses recorded in each model time step.
#'
#' The column `late_diagnoses` reports the number who are initiated to ART
#' directly from the undiagnosed population.
#'
#' The two approaches are largely for debugging purposes to understand
#' the implications of the model choices.
#'
#' @export
#' 
number_diagnoses <- function(mod, fp, df, VERSION = "R") {

  if(VERSION != "R")

    stop("C version not yet implemented")

  else { 

    diagnoses <- diagnoses_mod <- late_diagnoses <- numeric(length(df$haidx))

    for(i in seq_along(df$haidx)) {

      haidx <- df$haidx[i] + 1:df$hagspan[i] - 1
      sidx <- if(df$sidx[i] == 0) 1:2 else df$sidx[i]
      paidx <- fp$ss$agfirst.idx[df$haidx[i]] + 1:sum(fp$ss$h.ag.span[haidx]) - 1
      
      hivpop_ha_hm <- attr(mod, "hivpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i], drop = FALSE]
      diagn_ha_hm <- attr(mod, "diagnpop")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i], drop = FALSE]
      
      ## Calculation proportion ever tested among undiagnosed
      undiagnosed_ha_hm <- hivpop_ha_hm - diagn_ha_hm
      testneg_ha <- attr(mod, "testnegpop")[haidx, sidx, fp$ss$hivp.idx, df$yidx[i], drop = FALSE]
      prop_testneg <- 1 / colSums(hivpop_ha_hm - diagn_ha_hm) * c(testneg_ha)
      
      naive_ha_hm <- sweep(undiagnosed_ha_hm, 2:4, 1 - prop_testneg, "*")
      testneg_ha_hm <- sweep(undiagnosed_ha_hm, 2:4, prop_testneg, "*")
      
      diagnoses[i] <- sum(naive_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 1, df$yidx[i]])) +
        sum(testneg_ha_hm * c(fp$diagn_rate[ , haidx, sidx, 2, df$yidx[i]]))
      diagnoses_mod[i] <- sum(attr(mod, "diagnoses")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]])
      late_diagnoses[i] <- sum(attr(mod, "late_diagnoses")[ , df$haidx[i] + 1:df$hagspan[i] - 1, sidx, df$yidx[i]])
    }
  }
  
  data.frame(diagnoses = diagnoses,
             diagnoses_mod = diagnoses_mod,
             late_diagnoses = late_diagnoses)
}

#' ART coverage among age 15-49
#' @export
artcov15to49 <- function (mod, sex = 'both') {
    if (sex == 'both') { sexid <- 1:2 }
    if (sex == 'male') { sexid <- 1 }
    if (sex == 'female') { sexid <- 2 }  
    n_art <- colSums(attr(mod, "artpop")[, , 1:8, sexid, , drop = FALSE], 
        , 4)
    n_hiv <- colSums(attr(mod, "hivpop")[, 1:8, sexid, , drop = FALSE], 
        , 3)
    return(n_art/(n_hiv + n_art))
}
