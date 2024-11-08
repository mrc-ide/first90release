create_fp <- function(projp,
                      demp,
                      hiv_steps_per_year = 10L,
                      proj_start = projp$yr_start,
                      proj_end = projp$yr_end,
                      AGE_START = 15L,
                      popadjust = TRUE,
                      artelig200adj=TRUE,
                      who34percelig=0,
                      projection_period = NULL,
                      art_dropout_recover_cd4 = NULL) {


  ## ########################## ##
  ##  Define model state space  ##
  ## ########################## ##

  ## Parameters defining the model projection period and state-space
  ss <- list(proj_start = proj_start,
             PROJ_YEARS = as.integer(proj_end - proj_start + 1L),
             AGE_START  = as.integer(AGE_START),
             hiv_steps_per_year = as.integer(hiv_steps_per_year))

  ## populuation projection state-space
  ss$NG <- 2
  ss$pDS <- 2               # Disease stratification for population projection (HIV-, and HIV+)

  ## macros
  ss$m.idx <- 1
  ss$f.idx <- 2

  ss$hivn.idx <- 1
  ss$hivp.idx <- 2

  ss$pAG <- 81 - AGE_START
  ss$ag.rate <- 1
  ss$p.fert.idx <- 16:50 - AGE_START
  ss$p.age15to49.idx <- 16:50 - AGE_START
  ss$p.age15plus.idx <- (16-AGE_START):ss$pAG


  ## HIV model state-space
  ss$h.ag.span <- as.integer(c(2,3, rep(5, 6), 31))   # Number of population age groups spanned by each HIV age group [sum(h.ag.span) = pAG]
  ss$hAG <- length(ss$h.ag.span)          # Number of age groups
  ss$hDS <- 7                             # Number of CD4 stages (Disease Stages)
  ss$hTS <- 3                             # number of treatment stages (including untreated)

  ss$ag.idx <- rep(1:ss$hAG, ss$h.ag.span)
  ss$agfirst.idx <- which(!duplicated(ss$ag.idx))
  ss$aglast.idx <- which(!duplicated(ss$ag.idx, fromLast=TRUE))


  ss$h.fert.idx <- which((AGE_START-1 + cumsum(ss$h.ag.span)) %in% 15:49)
  ss$h.age15to49.idx <- which((AGE_START-1 + cumsum(ss$h.ag.span)) %in% 15:49)
  ss$h.age15plus.idx <- which((AGE_START-1 + cumsum(ss$h.ag.span)) >= 15)

  invisible(list2env(ss, environment())) # put ss variables in environment for convenience

  fp <- list(ss=ss)
  fp$SIM_YEARS <- ss$PROJ_YEARS
  fp$proj.steps <- proj_start + 0.5 + 0:(ss$hiv_steps_per_year * (fp$SIM_YEARS-1)) / ss$hiv_steps_per_year

  ## Replace 6,14 with 6.14 for files saved with french locale
  projp$spectrum_version <- sub("^([0-9]+),(.*)$", "\\1.\\2", projp$spectrum_version)
  
  if (!grepl("^[4-6]\\.[0-9]", projp$spectrum_version)) {
    stop(paste0("Spectrum version not recognized: ", projp$spectrum_version))
  }

  if (is.null(projection_period)) {
    fp$projection_period <- if (projp$spectrum_version >= "6.2") {"calendar"} else {"midyear"}
  } else {
    stopifnot(projection_period %in% c("calendar", "midyear"))
    fp$projection_period <- projection_period
  }
  
  
  ## ######################## ##
  ##  Demographic parameters  ##
  ## ######################## ##

  fp$basepop <- demp$basepop[(AGE_START+1):81, , as.character(proj_start)]
  fp$Sx <- demp$Sx[(AGE_START+1):81,,as.character(proj_start:proj_end)]

  fp$asfr <- demp$asfr[,as.character(proj_start:proj_end)] # NOTE: assumes 15-49 is within projection age range
  ## Note: Spectrum averages ASFRs from the UPD file over 5-year age groups.
  ##       Prefer to use single-year of age ASFRs as provided. The below line will
  ##       convert to 5-year average ASFRs to exactly match Spectrum.
  ## fp$asfr <- apply(apply(fp$asfr, 2, tapply, rep(3:9*5, each=5), mean), 2, rep, each=5)

  fp$srb <- sapply(demp$srb[as.character(proj_start:proj_end)], function(x) c(x,100)/(x+100))

  netmigr.adj <- demp$netmigr

  if (fp$projection_period == "midyear") {
    
    ## Spectrum mid-year projection (v5.19 and earlier) adjusts net-migration to occur
    ## half in current age group and half in next age group
    
    netmigr.adj[-1,,] <- (demp$netmigr[-1,,] + demp$netmigr[-81,,])/2
    netmigr.adj[1,,] <- demp$netmigr[1,,]/2
    netmigr.adj[81,,] <- netmigr.adj[81,,] + demp$netmigr[81,,]/2
  }
    
  fp$netmigr <- netmigr.adj[(AGE_START+1):81,,as.character(proj_start:proj_end)]

  fp$entrantpop <- projp$totpop[AGE_START,,as.character(proj_start:proj_end)]
  
  ## set population adjustment
  fp$popadjust <- popadjust

  hivnpop <- projp$totpop - projp$hivpop
  fp$target_hivn_pop <- hivnpop[(AGE_START+1):81,,as.character(proj_start:proj_end)]
  fp$target_hivp_pop <- projp$hivpop[(AGE_START+1):81,,as.character(proj_start:proj_end)]

  fp$target_hivn_ha <- apply(fp$target_hivn_pop, 2:3, fastmatch::ctapply, fp$ss$ag.idx, sum)
  fp$target_artpop_ha <- apply(projp$artpop[(AGE_START+1):81,,as.character(proj_start:proj_end)],
                               2:3, fastmatch::ctapply, fp$ss$ag.idx, sum)
  hivp_ha <- apply(fp$target_hivp_pop, 2:3, fastmatch::ctapply, fp$ss$ag.idx, sum)
  fp$target_hivpop_ha <- hivp_ha - fp$target_artpop_ha
  
  
  ## ###################### ##
  ##  HIV model parameters  ##
  ## ###################### ##

  projp.h.ag <- findInterval(AGE_START + cumsum(h.ag.span) - h.ag.span, c(15, 25, 35, 45))  # NOTE: Will not handle AGE_START < 15 presently
  fp$cd4_initdist <- projp$cd4_initdist[,projp.h.ag,]
  fp$cd4_prog <- (1-exp(-projp$cd4_prog[,projp.h.ag,] / hiv_steps_per_year)) * hiv_steps_per_year
  fp$cd4_mort <- projp$cd4_mort[,projp.h.ag,]
  fp$art_mort <- projp$art_mort[,,projp.h.ag,]
  fp$artmx_timerr <- projp$artmx_timerr

  fp$cd4_nonaids_excess_mort <- projp$cd4_nonaids_excess_mort[,projp.h.ag,]
  fp$art_nonaids_excess_mort <- array(0.0, dim(fp$art_mort))
  fp$art_nonaids_excess_mort[] <- rep(projp$art_nonaids_excess_mort[,projp.h.ag,], each = hTS)

  frr_agecat <- as.integer(rownames(projp$fert_rat))
  frr_agecat[frr_agecat == 18] <- 17
  fert_rat.h.ag <- findInterval(AGE_START + cumsum(h.ag.span[h.fert.idx]) - h.ag.span[h.fert.idx], frr_agecat)

  fp$frr_cd4 <- array(1, c(hDS, length(h.fert.idx), PROJ_YEARS))
  fp$frr_cd4[,,] <- rep(projp$fert_rat[fert_rat.h.ag, as.character(proj_start:proj_end)], each=hDS)
  fp$frr_cd4 <- sweep(fp$frr_cd4, 1, projp$cd4fert_rat, "*")
  fp$frr_cd4 <- fp$frr_cd4 * projp$frr_scalar

  fp$frr_art <- array(1.0, c(hTS, hDS, length(h.fert.idx), PROJ_YEARS))
  fp$frr_art[1,,,] <- fp$frr_cd4 # 0-6 months
  fp$frr_art[2:3, , , ] <- sweep(fp$frr_art[2:3, , , ], 3, projp$frr_art6mos[fert_rat.h.ag] * projp$frr_scalar, "*") # 6-12mos, >1 years

  ## ART eligibility and numbers on treatment

  fp$art15plus_num <- projp$art15plus_num[,as.character(proj_start:proj_end)]
  fp$art15plus_isperc <- projp$art15plus_numperc[, as.character(proj_start:proj_end)] == 1

  ## convert percentage to proportion
  fp$art15plus_num[fp$art15plus_isperc] <- fp$art15plus_num[fp$art15plus_isperc] / 100

  ## eligibility starts in projection year idx
  fp$specpop_percelig <- rowSums(with(projp$artelig_specpop[-1,], mapply(function(elig, percent, year) rep(c(0, percent*as.numeric(elig)), c(year - proj_start, proj_end - year + 1)), elig, percent, year)))
  fp$artcd4elig_idx <- findInterval(-projp$art15plus_eligthresh[as.character(proj_start:proj_end)], -c(999, 500, 350, 250, 200, 100, 50))

  ## Update eligibility threshold from CD4 <200 to <250 to account for additional
  ## proportion eligible with WHO Stage 3/4.
  if(artelig200adj)
    fp$artcd4elig_idx <- replace(fp$artcd4elig_idx, fp$artcd4elig_idx==5L, 4L)

  fp$pw_artelig <- with(projp$artelig_specpop["PW",], rep(c(0, elig), c(year - proj_start, proj_end - year + 1)))  # are pregnant women eligible (0/1)

  ## percentage of those with CD4 <350 who are based on WHO Stage III/IV infection
  fp$who34percelig <- who34percelig

  if (is.null(art_dropout_recover_cd4)) {
    fp$art_dropout_recover_cd4 <- projp$spectrum_version >= "6.14"
  } else {
    stopifnot(is.logical(art_dropout_recover_cd4))
    fp$art_dropout_recover_cd4 <- art_dropout_recover_cd4
  }

  
  fp$art_dropout <- projp$art_dropout[as.character(proj_start:proj_end)]/100
  fp$median_cd4init <- projp$median_cd4init[as.character(proj_start:proj_end)]
  fp$med_cd4init_input <- as.integer(fp$median_cd4init > 0)
  fp$med_cd4init_cat <- replace(findInterval(-fp$median_cd4init, - c(1000, 500, 350, 250, 200, 100, 50)),
                                !fp$med_cd4init_input, 0L)

  fp$tARTstart <- min(unlist(apply(fp$art15plus_num > 0, 1, which)), PROJ_YEARS)

  ## New ART patient allocation options
  fp$art_alloc_method <- projp$art_alloc_method
  fp$art_alloc_mxweight <- projp$art_prop_alloc[1]

  ## Scale mortality among untreated population by ART coverage
  fp$scale_cd4_mort <- projp$scale_cd4_mort

  ## HIV prevalence and ART coverage among age 15 entrants
  hivpop14 <- projp$age14hivpop[,,,as.character(proj_start:(proj_end-1))]
  pop14 <- projp$age14totpop[ , as.character(proj_start:(proj_end-1))]
  hiv14 <- colSums(hivpop14,,2)
  art14 <- colSums(hivpop14[5:7,,,],,2)

  fp$entrantprev <- cbind(0, hiv14/pop14) # 1 year offset because age 15 population is age 14 in previous year
  fp$entrantartcov <- cbind(0, art14/hiv14)
  fp$entrantartcov[is.na(fp$entrantartcov)] <- 0
  colnames(fp$entrantprev) <- colnames(fp$entrantartcov) <- as.character(proj_start:proj_end)

  hiv_noart14 <- colSums(hivpop14[1:4,,,])
  artpop14 <- hivpop14[5:7,,,]

  fp$paedsurv_cd4dist <- array(0, c(hDS, NG, PROJ_YEARS))
  fp$paedsurv_artcd4dist <- array(0, c(hTS, hDS, NG, PROJ_YEARS))

  cd4convert <- rbind(c(1, 0, 0, 0, 0, 0, 0),
                      c(1, 0, 0, 0, 0, 0, 0),
                      c(1, 0, 0, 0, 0, 0, 0),
                      c(0, 1, 0, 0, 0, 0, 0),
                      c(0, 0, 0.67, 0.33, 0, 0, 0),
                      c(0, 0, 0, 0, 0.35, 0.21, 0.44))

  ## Convert age 5-14 CD4 distribution to adult CD4 distribution and normalize to
  ## sum to 1 in each sex and year.
  for(g in 1:NG)
    for(i in 2:PROJ_YEARS){

      if((hiv14[g,i-1] - art14[g,i-1]) > 0) {
        fp$paedsurv_cd4dist[,g,i] <- hiv_noart14[,g,i-1] %*% cd4convert / (hiv14[g,i-1] - art14[g,i-1])
      }
      if(art14[g,i-1]) {
        fp$paedsurv_artcd4dist[,,g,i] <- artpop14[,,g,i-1] %*% cd4convert / art14[g,i-1]

        ## if age 14 has ART population in CD4 above adult eligibilty, assign to highest adult
        ## ART eligibility category.
        idx <- fp$artcd4elig_idx[i]
        if(idx > 1) {
          fp$paedsurv_artcd4dist[,idx,g,i] <- fp$paedsurv_artcd4dist[,idx,g,i] + rowSums(fp$paedsurv_artcd4dist[,1:(idx-1),g,i, drop=FALSE])
          fp$paedsurv_artcd4dist[,1:(idx-1),g,i] <- 0
        }
      }
    }

  class(fp) <- "specfp"

  return(fp)
}


#'  Beers coefficients to distribute from 5-year to single-year of age
#'
create_beers <- function(n5yr) {

  ## Beer's coefficients for disaggregating 5 year age groups into
  ## single-year age groups (from John Stover)
  Afirst <- rbind(c(0.3333, -0.1636, -0.0210,  0.0796, -0.0283),
                  c(0.2595, -0.0780,  0.0130,  0.0100, -0.0045),
                  c(0.1924,  0.0064,  0.0184, -0.0256,  0.0084),
                  c(0.1329,  0.0844,  0.0054, -0.0356,  0.0129),
                  c(0.0819,  0.1508, -0.0158, -0.0284,  0.0115))
  Asecond <- rbind(c( 0.0404,  0.2000, -0.0344, -0.0128,  0.0068),
                   c( 0.0093,  0.2268, -0.0402,  0.0028,  0.0013),
                   c(-0.0108,  0.2272, -0.0248,  0.0112, -0.0028),
                   c(-0.0198,  0.1992,  0.0172,  0.0072, -0.0038),
                   c(-0.0191,  0.1468,  0.0822, -0.0084, -0.0015))
  Amid <- rbind(c(-0.0117,  0.0804,  0.1570, -0.0284,  0.0027),
                c(-0.0020,  0.0160,  0.2200, -0.0400,  0.0060),
                c( 0.0050, -0.0280,  0.2460, -0.0280,  0.0050),
                c( 0.0060, -0.0400,  0.2200,  0.0160, -0.0020),
                c( 0.0027, -0.0284,  0.1570,  0.0804, -0.0117))
  Apenult <- rbind(c(-0.0015, -0.0084,  0.0822,  0.1468, -0.0191),
                   c(-0.0038,  0.0072,  0.0172,  0.1992, -0.0198),
                   c(-0.0028,  0.0112, -0.0248,  0.2272, -0.0108),
                   c( 0.0013,  0.0028, -0.0402,  0.2268,  0.0093),
                   c( 0.0068, -0.0128, -0.0344,  0.2000,  0.0404))
  Aultim <- rbind(c( 0.0115, -0.0284, -0.0158,  0.1508,  0.0819),
                  c( 0.0129, -0.0356,  0.0054,  0.0844,  0.1329),
                  c( 0.0084, -0.0256,  0.0184,  0.0064,  0.1924),
                  c(-0.0045,  0.0100,  0.0130, -0.0780,  0.2595),
                  c(-0.0283,  0.0796, -0.0210, -0.1636,  0.3333))

  A <- do.call(rbind,
               c(list(cbind(Afirst, matrix(0, 5, n5yr-5)),
                      cbind(Asecond, matrix(0, 5, n5yr-5))),
                 lapply(0:(n5yr-6), function(i) cbind(matrix(0, 5, i), Amid, matrix(0, 5, (n5yr-5)-i))),
                 list(cbind(matrix(0, 5, n5yr-6), Apenult, matrix(0, 5, 1)),
                      cbind(matrix(0, 5, n5yr-6), Aultim, matrix(0, 5, 1)),
                      c(rep(0, n5yr-1), 1))))
  return(round(A, 4))
}
