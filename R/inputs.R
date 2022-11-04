#' Process programmatic data on number of tests
#' 
#' @export
select_prgmdata <- function(prgm_dat, cnt, age_group) {
  # any country with HTS data? if so, we select them
  if (any(prgm_dat$country == cnt)) {
    prg_dat <- subset(prgm_dat, country == cnt)
  # for plotting purposes, we calculate the total form years with sex-disaggregated data
  if (any(prg_dat$sex != "both")) {
      prg_dat_both <- subset(prg_dat, sex == "both")
      prg_dat_sex <- subset(prg_dat, sex == "male" | sex == "female")
      yr_both <- unique(prg_dat_both$year)
      yr_sex <- unique(prg_dat_sex$year)
      prg_dat_sex_both <- NULL
      for (i in 1:length(yr_sex)) {
        prg.i <- subset(prg_dat_sex, year == yr_sex[i])
        tot.i <- data.frame(country = cnt, year = yr_sex[i],
                            agegr = '15-99',
                            sex = 'both',
                            tot = sum(prg.i$tot),
                            totpos = sum(prg.i$totpos),
                            vct = sum(prg.i$vct),
                            vctpos = sum(prg.i$vctpos),
                            anc = prg.i$anc[prg.i$sex == 'female'],
                            ancpos = prg.i$ancpos[prg.i$sex == 'female'])
        value_verif <- prg_dat_both$totpos[prg_dat_both$year == yr_sex[i]]
        if (length(value_verif) > 0 & is.na(tot.i$totpos)) { 
              tot.i$totpos <- prg_dat_both$totpos[prg_dat_both$year == yr_sex[i]] }
        prg_dat_sex_both <- rbind(prg_dat_sex_both, tot.i)
      }
      prg_dat <- rbind(prg_dat_both[!(prg_dat_both$year %in% yr_sex), ], prg_dat_sex_both, prg_dat_sex)
    }  } 
  
  if (!any(prgm_dat$country == cnt)) {

    ## -- UPDATE HERE --
    ## * year vector needs to be extended to output results to current year
    
    prg_dat <- data.frame(country = cnt, 
                          year = 2010:2022,
                          agegr = '15-99', sex = 'both',
                          tot = NA, totpos = NA,
                          vct = NA, vctpos = NA, anc = NA, ancpos = NA)
    ## -- UPDATE ABOVE --
  }
  prg_dat <- prg_dat[order(prg_dat$year), ]
  return(prg_dat)
}

#' Process survey data on hiv testing behaviors
#' 
#' @export
select_hts <- function(survey_hts, cnt, age_group) {
  require(first90)

  # We first select stratified by HIV status for females
  dat_f <- survey_hts[survey_hts$country == cnt & (survey_hts$hivstatus != "all") & 
                      (survey_hts$agegr %in% age_group) & survey_hts$sex == "female" & 
                      survey_hts$outcome == "evertest", ]
  # We first select stratified by HIV status for males
  dat_m <- survey_hts[survey_hts$country == cnt & (survey_hts$hivstatus != "all") & 
                      (survey_hts$agegr %in% age_group) & survey_hts$sex == "male" & 
                      survey_hts$outcome == "evertest", ]

  if (cnt == 'Equatorial Guinea') {
    # special case of GE b/c data for HIV entered manually (not age disaggregated)
    dat_f <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & 
                          (survey_hts$agegr %in% age_group)  & survey_hts$sex == "female" & survey_hts$outcome == "evertest", ]
    
    dat_m <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & 
                          (survey_hts$agegr %in% age_group)  & survey_hts$sex == "male" & survey_hts$outcome == "evertest", ]
    
    dat_f_hiv <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "positive" & 
                          survey_hts$agegr == '15-49'  & survey_hts$sex == "female" & survey_hts$outcome == "evertest", ]
    
    dat_m_hiv <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "positive" & 
                          survey_hts$agegr == '15-49'  & survey_hts$sex == "male" & survey_hts$outcome == "evertest", ]
    
    dat_f <- rbind(dat_f, dat_f_hiv)
    dat_m <- rbind(dat_m, dat_m_hiv)
  }
  
  # We find surveys that collected information on HIV status but for different age groups
  all_svy_hiv <- unique(survey_hts$surveyid[survey_hts$country == cnt & (survey_hts$hivstatus != "all") &
                        survey_hts$outcome == "evertest"])
  other_svy_hiv <- all_svy_hiv[!(all_svy_hiv %in% unique(rbind(dat_f, dat_m)$surveyid))]
  
  dat_both <- NULL
  if (length(other_svy_hiv) > 0) {
    dat_f2 <- survey_hts[(survey_hts$surveyid %in% other_svy_hiv) & (survey_hts$hivstatus != "all") & 
                        survey_hts$agegr == '15-49' & survey_hts$sex == "female" & 
                        survey_hts$outcome == "evertest", ]
    dat_m2 <- survey_hts[(survey_hts$surveyid %in% other_svy_hiv) & (survey_hts$hivstatus != "all") & 
                 survey_hts$agegr == '15-49' & survey_hts$sex == "male" & 
                 survey_hts$outcome == "evertest", ]
    if (dim(dat_f2)[1] > 0 | dim(dat_m2)[1] > 0) {
      dat_f <- rbind(dat_f, dat_f2)
      dat_m <- rbind(dat_m, dat_m2)     
    } else {
      dat_both <- survey_hts[(survey_hts$surveyid %in% other_svy_hiv) & (survey_hts$hivstatus != "all") & 
                     survey_hts$agegr == '15-49' & survey_hts$sex == "both" & 
                     survey_hts$outcome == "evertest", ]
    }
  }
  
  # We find surveys that did not collect information on HIV status
  all_survey <- unique(survey_hts$surveyid[survey_hts$country == cnt & survey_hts$outcome == "evertest"])
  other_survey <-  all_survey[!(all_survey %in% unique(rbind(dat_f, dat_m, dat_both)$surveyid))]
  
  dat_fall <- NULL
  dat_mall <- NULL

  if(length(other_survey)){
    for (i in 1:length(other_survey)){
      age_strat_f <- unique(survey_hts$agegr[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                                             survey_hts$sex == "female" & survey_hts$outcome == "evertest"])
      age_strat_m <- unique(survey_hts$agegr[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                                             survey_hts$sex == "male" & survey_hts$outcome == "evertest"])    
      ## Female - If we have the correct age disagregation, we stick with it
      if (all(age_group %in% age_strat_f)) {
        dat_f_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% age_group) & survey_hts$sex == "female" & survey_hts$outcome == "evertest",]
      } else if (all(c('15-24','25-49') %in% age_strat_f)) {
        dat_f_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% c('15-24','25-49')) & survey_hts$sex == "female" & survey_hts$outcome == "evertest",]
      } else if ("15-49" %in% age_strat_f) {
        dat_f_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% c('15-49')) & survey_hts$sex == "female" & survey_hts$outcome == "evertest",]
      } else {
                dat_f_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% c('15-99')) & survey_hts$sex == "female" & survey_hts$outcome == "evertest",]
      } 
      ## Female - If we have the correct age disagregation, we stick with it
      if (all(age_group %in% age_strat_m)) {
        dat_m_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% age_group) & survey_hts$sex == "male" & survey_hts$outcome == "evertest",]
      } else if (all(c('15-24','25-49','15-99') %in% age_strat_m)) {
        dat_m_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% c('15-24','25-49')) & survey_hts$sex == "male" & survey_hts$outcome == "evertest",]
      } else if ("15-49" %in% age_strat_m) {
        dat_m_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% c('15-49')) & survey_hts$sex == "male" & survey_hts$outcome == "evertest",]
      } else {
        dat_m_i <- survey_hts[survey_hts$country == cnt & survey_hts$hivstatus == "all" & survey_hts$surveyid == other_survey[i] &
                              (survey_hts$agegr %in% c('15-99')) & survey_hts$sex == "male" & survey_hts$outcome == "evertest",]
      }
      dat_fall <- rbind(dat_fall, dat_f_i)
      dat_mall <- rbind(dat_mall, dat_m_i)
    }
  }

  dat <- rbind(dat_f, dat_m, dat_fall, dat_mall, dat_both)
  return(dat)
}

#' Function to present the survey data on shinny interface
#' 
#' @export
svy_hts_interface <- function(survey_hts, cnt, age_group = c('15-24','25-34','35-49')) {
  # Data used in the fitting
  svy_dat <- select_hts(survey_hts, cnt, age_group)
  svy_dat <- svy_dat[order(svy_dat$year, svy_dat$sex, svy_dat$hivstatus,
                    svy_dat$agegr), ]
  svy_dat$color <- "Fitting"
  # Data used for plotting
  svy_plot1549 <- subset(survey_hts, country == cnt & agegr == '15-49' &
                      sex != "both" & outcome == 'evertest')
  svy_plotage <- subset(survey_hts, country == cnt & agegr %in% age_group &
                                    sex != "both" & outcome == 'evertest' & hivstatus == 'all')
  svy_plot90s <- subset(survey_hts, country == cnt &
                                    agegr == '15-49' &
                                    sex == 'both' &
                                    (outcome %in% c('aware', 'art_both', 'art_self') |
                                     hivstatus == 'positive' & outcome == 'evertest'))
    
  svy_plot <- rbind(svy_plot1549, svy_plotage, svy_plot90s)
  svy_plot <- svy_plot[order(svy_plot$year, svy_plot$sex, svy_plot$hivstatus, svy_plot$agegr), ]
  svy_plot$color <- "Plotting"
  
  svy <- rbind(svy_dat, svy_plot)
  
  return(svy)
}


#' Add incidences corresponding to state space
#'
#' @param dat a data.frame with particular columns coded in particular ways
#' @param ss  state space definition (fp$ss)
#'
#' @examples
#' \dontrun{
#' data(survey_hts)
#' dat <- subset(survey_hts, survey_hts$country == "Malawi" & outcome == "evertest")
#' df <- add_ss_indices(dat, fp$ss)
#' }
#' @export
add_ss_indices <- function(dat, ss) {

  df <- dat

  ## Convert 15+ to 15-99 for parsing below
  agegr_tmp <- sub("15\\+", "15-99", df$agegr)
  
  df$agestart <- type.convert(sub("([0-9]+)-([0-9]+)", "\\1", agegr_tmp), as.is = TRUE)
  ageend <- type.convert(sub("([0-9]+)-([0-9]+)", "\\2", agegr_tmp), as.is = TRUE)+1L

  df$aidx <- df$agestart - ss$AGE_START + 1L
  df$agspan <- ageend - df$agestart

  ## # HIV age group indices
  hag_breaks <- ss$agfirst.idx + ss$AGE_START - 1L
  df$haidx <- match(df$agestart, hag_breaks)
  end_haidx <- match(ageend, hag_breaks)

  ## if max age is 65 or greater, assume age 50+
  end_haidx[ageend >= 65] <- length(hag_breaks) + 1L
  
  df$hagspan <- end_haidx - df$haidx
  
  ## testing history is homogenous among oldest age group
  df$haidx[df$agestart >= max(hag_breaks) & ageend >= max(hag_breaks)] <- length(hag_breaks)
  df$hagspan[df$agestart >= max(hag_breaks) & ageend >= max(hag_breaks)] <- 1L

  
 df$sidx <- match(df$sex, c("both", "male", "female")) - 1L
  df$hvidx <- match(df$hivstatus, c("all", "negative", "positive")) - 1L
  df$yidx <- df$year - ss$proj_start + 1L

  if(any(is.na(df$hagspan)))
    warning(paste(sum(is.na(df$hagspan)), "records with non-matching HIV age groups"))

  df
}


#' Create model inputs from aggregated Spectrum PJNZ
#'
#' @param pjnz_in a list of outputs from extract_pjnz
#'
#' @export
prepare_inputs_from_extracts <- function(pjnz_in){

  pjnz_aggr <- combine_inputs(pjnz_in)
  
  demp <- list(basepop = pjnz_aggr$totpop,
               Sx = pjnz_aggr$Sx,
               asfr = pjnz_aggr$asfr,
               srb = pjnz_aggr$srb,
               netmigr = pjnz_aggr$netmigr)
  
  projp <- pjnz_aggr
  projp$age14totpop <- projp$totpop["14",,]
  dimnames(projp$fert_rat)[[1]] <- gsub("(.*)-.*", "\\1", dimnames(projp$fert_rat)[[1]])
  
  fp <- create_fp(projp, demp)

  ## Set pop adjust
  fp$popadjust <- pjnz_aggr$popadjust
  
  ## Calculate incidence input
  fp$infections <- pjnz_aggr$infections[fp$ss$AGE_START + fp$ss$p.age15plus.idx, , ]
  fp$eppmod <- "directinfections_hts"
  
  ## initialise to no testing
  fp$hts_rate <- array(0.0, c(fp$ss$hAG, fp$ss$NG, fp$ss$pDS, fp$ss$PROJ_YEARS))
  fp$diagn_rate <- array(0.0, c(fp$ss$hDS, fp$ss$hAG, fp$ss$NG, 4, fp$ss$PROJ_YEARS))
  
  fp$t_hts_start <- fp$ss$PROJ_YEARS+1L
  
  fp
}


#' Create model inputs from Spectrum PJNZ
#'
#' @param pjnzlist a vector of PJNZ file names to aggregate
#'
#' @details
#'
#' The aggregation makes a number of assumptions:
#'
#' * Progression parameters are the same in all files, and values frome the first file are used.
#' * Special populations ART eligibility is the same in all files.
#'
#' @examples
#'
#' pjnzlist <- list.files("~/Documents/Data/Spectrum files/2018 final/SSA/", "CotedIvoire.*PJNZ$", full.names=TRUE, ignore.case=TRUE)
#' pjnzlist <- "~/Documents/Data/Spectrum files/2018 final/SSA/Malawi_2018_version_8.PJNZ"
#' 
#' @export
prepare_inputs <- function(pjnzlist){

  pjnz_in <- lapply(pjnzlist, extract_pjnz)
  prepare_inputs_from_extracts(pjnz_in)
}


add_fp <- function(v, fp){
  fp[names(v)] <- v
  fp
}
  

create_paedsurv_inputs <- function(projp, ss, artcd4elig_idx) {

  proj_start <- ss$proj_start
  proj_end <- proj_start + ss$PROJ_YEARS - 1L
  
  val <- list()

  ## HIV prevalence and ART coverage among age 15 entrants
  hivpop14 <- projp$age14hivpop[,,,as.character(proj_start:(proj_end-1))]
  pop14 <- projp$age14totpop[ , as.character(proj_start:(proj_end-1))]
  hiv14 <- colSums(hivpop14,,2)
  art14 <- colSums(hivpop14[5:7,,,],,2)

  val$entrantprev <- cbind(0, hiv14/pop14) # 1 year offset because age 15 population is age 14 in previous year
  val$entrantartcov <- cbind(0, art14/hiv14)
  val$entrantartcov[is.na(val$entrantartcov)] <- 0
  colnames(val$entrantprev) <- colnames(val$entrantartcov) <- as.character(proj_start:proj_end)

  hiv_noart14 <- colSums(hivpop14[1:4,,,])
  artpop14 <- hivpop14[5:7,,,]

  val$paedsurv_cd4dist <- array(0, c(ss$hDS, ss$NG, ss$PROJ_YEARS))
  val$paedsurv_artcd4dist <- array(0, c(ss$hTS, ss$hDS, ss$NG, ss$PROJ_YEARS))

  cd4convert <- rbind(c(1, 0, 0, 0, 0, 0, 0),
                      c(1, 0, 0, 0, 0, 0, 0),
                      c(1, 0, 0, 0, 0, 0, 0),
                      c(0, 1, 0, 0, 0, 0, 0),
                      c(0, 0, 0.67, 0.33, 0, 0, 0),
                      c(0, 0, 0, 0, 0.35, 0.21, 0.44))

  ## Convert age 5-14 CD4 distribution to adult CD4 distribution and normalize to
  ## sum to 1 in each sex and year.
  for(g in 1:ss$NG)
    for(i in 2:ss$PROJ_YEARS){
      
      if((hiv14[g,i-1] - art14[g,i-1]) > 0)
        val$paedsurv_cd4dist[,g,i] <- hiv_noart14[,g,i-1] %*% cd4convert / (hiv14[g,i-1] - art14[g,i-1])
      if(art14[g,i-1]){
        val$paedsurv_artcd4dist[,,g,i] <- artpop14[,,g,i-1] %*% cd4convert / art14[g,i-1]

        ## if age 14 has ART population in CD4 above adult eligibilty, assign to highest adult
        ## ART eligibility category.
        idx <- artcd4elig_idx[i]
        if(idx > 1){
          val$paedsurv_artcd4dist[,idx,g,i] <- val$paedsurv_artcd4dist[,idx,g,i] + rowSums(val$paedsurv_artcd4dist[,1:(idx-1),g,i, drop=FALSE])
          val$paedsurv_artcd4dist[,1:(idx-1),g,i] <- 0
        }
      }
    }

  val
}

.fp_targetpop <- function(totpop, ss){

  v <- list()
  
  proj_start <- ss$proj_start
  proj_end <- proj_start + ss$PROJ_YEARS - 1L

  if(!length(setdiff(proj_start:proj_end, dimnames(totpop)[[3]]))){
    v$entrantpop <- totpop[ss$AGE_START,,as.character(proj_start:proj_end)]
    v$targetpop <- totpop[ss$AGE_START+1:ss$pAG,,as.character(proj_start:proj_end)]
  } else
    stop("targetpop does not span proj_start:proj_end")

  v
}
