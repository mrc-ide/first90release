#' Simulate EPP-ASM model
#'
#' @useDynLib first90 eppasmC
#' @export
simmod <- function(fp, VERSION = "C") {

  if(VERSION != "R") {
    fp$eppmodInt <- match(fp$eppmod, c("rtrend", "directincid", "directinfections", "directinfections_hts"), nomatch = 0) # 0: r-spline;

    ## projection_period codes:
    ## 0: mid-year (<= Spectrum 5.19)
    ## 1: calendar year (>= Spectrum 5.2)
    fp$projection_period_int <- match(fp$projection_period, c("midyear", "calendar")) - 1L # -1 for 0-based indexing

    mod <- .Call(eppasmC, fp)
    class(mod) <- "spec"
    return(mod)
  }

  ##################################################################################

  if(requireNamespace("fastmatch", quietly = TRUE))
    ctapply <- fastmatch::ctapply
  else
    ctapply <- tapply

  fp$ss$DT <- 1/fp$ss$hiv_steps_per_year

  ## Attach state space variables
  invisible(list2env(fp$ss, environment())) # put ss variables in environment for convenience

  ## initialize projection
  pop <- array(0, c(pAG, NG, pDS, PROJ_YEARS))
  pop[,,1,1] <- fp$basepop
  hivpop <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  diagnpop <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  artpop <- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))

  testnegpop <- array(0, c(hAG, NG, pDS, PROJ_YEARS))

  ## initialize output
  prev15to49 <- numeric(PROJ_YEARS)
  incid15to49 <- numeric(PROJ_YEARS)
  sexinc15to49out <- array(NA, c(NG, PROJ_YEARS))
  paedsurvout <- rep(NA, PROJ_YEARS)

  infections <- array(0, c(pAG, NG, PROJ_YEARS))
  hivdeaths <- array(0, c(pAG, NG, PROJ_YEARS))
  natdeaths <- array(0, c(pAG, NG, PROJ_YEARS))

  excessnonaidsdeaths <- array(0.0, c(pAG, NG, PROJ_YEARS))

  aidsdeaths_noart <- array(0.0, c(hDS, hAG, NG, PROJ_YEARS))
  natdeaths_noart <- array(0.0, c(hDS, hAG, NG, PROJ_YEARS))
  excessnonaidsdeaths_noart <- array(0.0, c(hDS, hAG, NG, PROJ_YEARS))

  aidsdeaths_art <- array(0.0, c(hTS, hDS, hAG, NG, PROJ_YEARS))
  natdeaths_art <- array(0.0, c(hTS, hDS, hAG, NG, PROJ_YEARS))
  excessnonaidsdeaths_art <- array(0.0, c(hTS, hDS, hAG, NG, PROJ_YEARS))


  hivpopdeaths <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  artpopdeaths <- array(0, c(hTS, hDS, hAG, NG, PROJ_YEARS))

  hivtests <- array(0, c(hAG, NG, 6, PROJ_YEARS))
  diagnoses <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  late_diagnoses <- array(0, c(hDS, hAG, NG, PROJ_YEARS))
  artinits <- array(0, c(hDS, hAG, NG, PROJ_YEARS))

  popadj.prob <- array(0, c(pAG, NG, PROJ_YEARS))

  entrant_prev_out <- numeric(PROJ_YEARS)
  hivp_entrants_out <- array(0, c(NG, PROJ_YEARS))

  if(fp$eppmod != "directincid") {
    ## outputs by timestep
    incrate15to49.ts.out <- rep(NA, length(fp$rvec))
    rvec <- if(fp$eppmod == "rtrend") rep(NA, length(fp$proj.steps)) else fp$rvec

    prev15to49.ts.out <- rep(NA, length(fp$rvec))
  }

  ## store last prevalence value (for r-trend model)
  prevlast <- 0


  for(i in 2:fp$SIM_YEARS){

    ## ################################### ##
    ##  Single-year population projection  ##
    ## ################################### ##

    ## age the population
    pop[-c(1,pAG),,,i] <- pop[-(pAG-1:0),,,i-1]
    pop[pAG,,,i] <- pop[pAG,,,i-1] + pop[pAG-1,,,i-1] # open age group

    ## Add lagged births into youngest age group
    hivn_entrants <- fp$entrantpop[,i-1] * (1 - fp$entrantprev[ , i])
    hivp_entrants <- fp$entrantpop[,i-1] * fp$entrantprev[ , i]

    entrant_prev_out[i] <- sum(hivp_entrants) / sum(hivn_entrants + hivp_entrants)
    hivp_entrants_out[,i] <- sum(hivp_entrants)

    pop[1,,hivn.idx,i] <- hivn_entrants
    pop[1,,hivp.idx,i] <- hivp_entrants

    ## if(i > fp$t_hts_start){
    ##   ## !! DO SOMETHING ABOUT HTS AMONG NEW ENTRANTS
    ## }

    hiv.ag.prob <- pop[aglast.idx,,hivp.idx,i-1] / apply(pop[,,hivp.idx,i-1], 2, ctapply, ag.idx, sum)
    hiv.ag.prob[is.nan(hiv.ag.prob)] <- 0

    hivpop[,,,i] <- hivpop[,,,i-1]
    hivpop[,-hAG,,i] <- hivpop[,-hAG,,i] - sweep(hivpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
    hivpop[,-1,,i] <- hivpop[,-1,,i] + sweep(hivpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
    hivpop[,1,,i] <- hivpop[,1,,i] + sweep(fp$paedsurv_cd4dist[,,i], 2, hivp_entrants * (1-fp$entrantartcov[,i]), "*")

    if(i > fp$t_hts_start) {

      hivn.ag.prob <- pop[aglast.idx,,hivn.idx,i-1] / apply(pop[,,hivn.idx,i-1], 2, ctapply, ag.idx, sum)
      hivn.ag.prob[is.nan(hivn.ag.prob)] <- 0

      testnegpop[,,,i] <- testnegpop[,,,i-1]
      testnegpop[-hAG,,hivn.idx,i] <- testnegpop[-hAG,,hivn.idx,i] - testnegpop[-hAG,,hivn.idx,i-1] * hivn.ag.prob[-hAG,]
      testnegpop[-1,,hivn.idx,i] <- testnegpop[-1,,hivn.idx,i] + testnegpop[-hAG,,hivn.idx,i-1] * hivn.ag.prob[-hAG,]

      testnegpop[-hAG,,hivp.idx,i] <- testnegpop[-hAG,,hivp.idx,i] - testnegpop[-hAG,,hivp.idx,i-1] * hiv.ag.prob[-hAG,]
      testnegpop[-1,,hivp.idx,i] <- testnegpop[-1,,hivp.idx,i] + testnegpop[-hAG,,hivp.idx,i-1] * hiv.ag.prob[-hAG,]

      ## ageing diagnosed population
      diagnpop[,,,i] <- diagnpop[,,,i-1]
      diagnpop[,-hAG,,i] <- diagnpop[,-hAG,,i] - sweep(diagnpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")
      diagnpop[,-1,,i] <- diagnpop[,-1,,i] + sweep(diagnpop[,-hAG,,i-1], 2:3, hiv.ag.prob[-hAG,], "*")

      ## !! Currently assume that age 15 entrants are not previously tested negative or diagnosed and untreated
    }

    if(i > fp$tARTstart) {
      artpop[,,,,i] <- artpop[,,,,i-1]
      artpop[,,-hAG,,i] <- artpop[,,-hAG,,i] - sweep(artpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
      artpop[,,-1,,i] <- artpop[,,-1,,i] + sweep(artpop[,,-hAG,,i-1], 3:4, hiv.ag.prob[-hAG,], "*")
      artpop[,,1,,i] <- artpop[,,1,,i] + sweep(fp$paedsurv_artcd4dist[,,,i], 3, hivp_entrants * fp$entrantartcov[,i], "*")
    }

    ## survive the population
    deaths <- sweep(pop[,,,i], 1:2, (1-fp$Sx[,,i]), "*")
    hiv.sx.prob <- 1-apply(deaths[,,hivp.idx], 2, ctapply, ag.idx, sum) / apply(pop[,,hivp.idx,i], 2, ctapply, ag.idx, sum)
    hiv.sx.prob[is.nan(hiv.sx.prob)] <- 0

    if(i > fp$t_hts_start) {
      hivn.sx.prob <- 1-apply(deaths[,,hivn.idx], 2, ctapply, ag.idx, sum) / apply(pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum)
      hivn.sx.prob[is.nan(hivn.sx.prob)] <- 0
    }

    pop[,,,i] <- pop[,,,i] - deaths
    natdeaths[,,i] <- rowSums(deaths,,2)

    hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.sx.prob, "*")

    if(i > fp$t_hts_start) {
      testnegpop[,, hivn.idx,i] <- testnegpop[,,hivn.idx,i] * hivn.sx.prob
      testnegpop[,, hivp.idx,i] <- testnegpop[,,hivp.idx,i] * hiv.sx.prob

      diagnpop[,,,i] <- sweep(diagnpop[,,,i], 2:3, hiv.sx.prob, "*")
    }
    if(i > fp$tARTstart)
      artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.sx.prob, "*")

    if (fp$projection_period == "midyear") {
      
      ## net migration
      netmigsurv <- fp$netmigr[,,i]*(1+fp$Sx[,,i])/2
      mr.prob <- 1+netmigsurv / rowSums(pop[,,,i],,2)
      hiv.mr.prob <- apply(mr.prob * pop[,,hivp.idx,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,hivp.idx,i], 2, ctapply, ag.idx, sum)
      hiv.mr.prob[is.nan(hiv.mr.prob)] <- 0

      if(i > fp$t_hts_start) {
        hivn.mr.prob <- apply(mr.prob * pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum)
        hivn.mr.prob[is.nan(hivn.mr.prob)] <- 0
      }
      
      pop[,,,i] <- sweep(pop[,,,i], 1:2, mr.prob, "*")
      
      hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.mr.prob, "*")
      if(i > fp$t_hts_start) {
        testnegpop[,, hivn.idx,i] <- testnegpop[,,hivn.idx,i] * hivn.mr.prob
        testnegpop[,, hivp.idx,i] <- testnegpop[,,hivp.idx,i] * hiv.mr.prob
        
        diagnpop[,,,i] <- sweep(diagnpop[,,,i], 2:3, hiv.mr.prob, "*")
      }
      if(i > fp$tARTstart)
        artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.mr.prob, "*")
    }
    
    ## fertility
    births.by.age <- rowSums(pop[p.fert.idx, f.idx,,i-1:0])/2 * fp$asfr[,i]
    births.by.h.age <- ctapply(births.by.age, ag.idx[p.fert.idx], sum)
    births <- fp$srb[,i] * sum(births.by.h.age)


    ## ################################ ##
    # ----  Disease model simulation  ---- 
    ## ################################ ##

    ## events at dt timestep
    for(ii in seq_len(hiv_steps_per_year)) {

      ts <- (i-2)/DT + ii

      grad <- array(0, c(hDS, hAG, NG))

      ## disease progression and mortality
      grad[-hDS,,] <- grad[-hDS,,] - fp$cd4_prog * hivpop[-hDS,,,i]  # remove cd4 stage progression (untreated)
      grad[-1,,] <- grad[-1,,] + fp$cd4_prog * hivpop[-hDS,,,i]      # add cd4 stage progression (untreated)

      if(fp$scale_cd4_mort == 1) {
        cd4mx_scale <- hivpop[,,,i] / (hivpop[,,,i] + colSums(artpop[,,,,i]))
        cd4mx_scale[!is.finite(cd4mx_scale)] <- 1.0
        cd4_mort_ts <- fp$cd4_mort * cd4mx_scale
      } else {
        cd4_mort_ts <- fp$cd4_mort
      }

      hivdeaths.ts <- cd4_mort_ts * hivpop[,,,i]  # HIV mortality, untreated
      grad <- grad - hivdeaths.ts              
      hivdeaths_hAG.ts <- colSums(hivdeaths.ts)
      aidsdeaths_noart[,,,i] <- aidsdeaths_noart[,,,i] + DT * hivdeaths.ts
      
      ## Non-AIDS excess mortality
      nonaids_excess.ts <- fp$cd4_nonaids_excess_mort * hivpop[,,,i]
      grad <- grad - nonaids_excess.ts
      nonaids_excess_hAG.ts <- DT * colSums(nonaids_excess.ts)
      excessnonaidsdeaths_noart[,,,i] <- excessnonaidsdeaths_noart[,,,i] + DT * nonaids_excess.ts


      ## ---- Distributing new infections in disease model ----
      if(fp$eppmod == "directinfections_hts") {

        ## Calculate annualised new infections by HIV age groups
        infections_ha <- apply(fp$infections[,, i], 2, ctapply, ag.idx, sum)
        grad <- grad + sweep(fp$cd4_initdist, 2:3, infections_ha, "*")

        ## move new infections per DT from negative to positive compartment
        pop[,, hivn.idx, i] <- pop[,, hivn.idx, i] - DT * fp$infections[,, i]
        pop[,, hivp.idx, i] <- pop[,, hivp.idx, i] + DT * fp$infections[,, i]
        infections[,,i] <- infections[,,i] + DT * fp$infections[,, i]
      }      

      ## HIV testing and diagnosis
      if(i >= fp$t_hts_start) {

        grad_tn <- array(0, c(hAG, NG, pDS))
        grad_diagn <- array(0, c(hDS, hAG, NG))

        ## Number of tests among never tested HIV negative persons
        hivn_pop_ha <- apply(pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum)

        hivtests[,,1,i] <- hivtests[,,1,i] + DT * (fp$hts_rate[ , , 1, i] * (hivn_pop_ha - testnegpop[ , , hivn.idx , i]))
        hivtests[,,2,i] <- hivtests[,,2,i] + DT * fp$hts_rate[ , , 2, i] * testnegpop[ , , hivn.idx , i]
        grad_tn[ , , hivn.idx] <- grad_tn[ , , hivn.idx] + fp$hts_rate[ , , 1, i] * (hivn_pop_ha - testnegpop[ , , hivn.idx , i])

        ## Add new infections to testneg population
        if (fp$eppmod == "directinfections_hts") {

          infections_ha <- apply(fp$infections[,, i], 2, ctapply, ag.idx, sum)
          
          ## number of new infections among never tested (proportional to the HIV-neg pop)
          hivn_pop_ha <- apply(pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum)
          testneg_infections_ha <- infections_ha * (testnegpop[,,hivn.idx,i] / hivn_pop_ha)

          grad_tn[ , , hivn.idx] <- grad_tn[ , , hivn.idx] - testneg_infections_ha
          grad_tn[ , , hivp.idx] <- grad_tn[ , , hivp.idx] + testneg_infections_ha
        }

        ## Do new diagnoses

        ## Remove HIV deaths among tested negative pop
        prop_tn_hivp <- testnegpop[,,hivp.idx,i] / colSums(hivpop[,,,i])
        prop_tn_hivp[!is.finite(prop_tn_hivp)] <- 0.0
        grad_tn[ , , hivp.idx] <- grad_tn[ , , hivp.idx] - hivdeaths_hAG.ts * prop_tn_hivp

        ## Remove non-AIDS excess deaths among tested negative pop
        grad_tn[ , , hivp.idx] <- grad_tn[ , , hivp.idx] - excess_nonaids_hAG.ts * prop_tn_hivp

        undiagnosed_i <- (hivpop[,,,i] - diagnpop[,,,i])
        prop_testneg <- testnegpop[ , , hivp.idx, i] / colSums(undiagnosed_i)
        prop_testneg[is.na(prop_testneg) | prop_testneg > 1 | prop_testneg < 0] <- 0

        ## Annualized new diagnoses among never tested population
        diagn_naive <- fp$diagn_rate[,,,1,i] * sweep(undiagnosed_i, 2:3, 1 - prop_testneg, "*")

        ## Annualized new diagnoses among previously tested negative population
        diagn_testneg <- fp$diagn_rate[,,,2,i] * sweep(undiagnosed_i, 2:3, prop_testneg, "*")

        hivtests[,,3,i] <- hivtests[,,3,i] + DT * colSums(diagn_naive)
        hivtests[,,4,i] <- hivtests[,,4,i] + DT * colSums(diagn_testneg)
        hivtests[,,5,i] <- hivtests[,,5,i] + DT * colSums(fp$diagn_rate[,,,3,i] * diagnpop[,,,i])
        hivtests[,,6,i] <- hivtests[,,6,i] + DT * colSums(fp$diagn_rate[,,,4,i] * colSums(artpop[,,,,i]))

        grad_tn[,,hivp.idx] <- grad_tn[,,hivp.idx] - colSums(diagn_testneg)
        grad_diagn <- grad_diagn + diagn_naive + diagn_testneg

        ## Disease progression and mortality among diagnosed (untreated) population
        grad_diagn[-hDS,,] <- grad_diagn[-hDS,,] - fp$cd4_prog * diagnpop[-hDS,,,i]  # remove cd4 stage progression (untreated)
        grad_diagn[-1,,] <- grad_diagn[-1,,] + fp$cd4_prog * diagnpop[-hDS,,,i]      # add cd4 stage progression (untreated)

        grad_diagn <- grad_diagn - fp$cd4_mort * diagnpop[,,,i]                      # HIV mortality, untreated
        grad_diagn <- grad_diagn - fp$cd4_nonaids_excess_mort * diagnpop[,,,i]       # Non-AIDS excess mortality

        diagnoses[,,,i] <- diagnoses[,,,i] + DT * (diagn_naive + diagn_testneg)

        testnegpop[ , , , i] <- testnegpop[ , , , i] + DT * grad_tn
        diagnpop[,,,i] <- diagnpop[,,,i] + DT * grad_diagn
      }

      hivpop[,,,i] <- hivpop[,,,i] + DT*grad
      hivpopdeaths[,,, i] <- hivpopdeaths[,,, i] + DT * (hivdeaths.ts + nonaids_excess.ts)


      ## ART population
      if(i >= fp$tARTstart) {

        gradART <- array(0, c(hTS, hDS, hAG, NG))

        ## progression and mortality
        gradART[1:2,,,] <- gradART[1:2,,,] - 2.0 * artpop[1:2,,,, i]      # remove ART duration progression (HARD CODED 6 months duration)
        gradART[2:3,,,] <- gradART[2:3,,,] + 2.0 * artpop[1:2,,,, i]      # add ART duration progression (HARD CODED 6 months duration)

        artdeaths.ts <- fp$art_mort * fp$artmx_timerr[ , i] * artpop[,,,,i]
        gradART <- gradART - artdeaths.ts                  # ART mortality

        hivdeaths_hAG.ts <- hivdeaths_hAG.ts + colSums(artdeaths.ts,,2)

        nonaids_excess_onart.ts <- fp$art_nonaids_excess_mort * artpop[,,,,i]
        gradART <- gradART - nonaids_excess_onart.ts
        nonaids_excess_hAG.ts <- nonaids_excess_hAG.ts +
          DT * colSums(nonaids_excess_onart.ts,,2)
        excessnonaidsdeaths_art[,,,,i] <- excessnonaidsdeaths_art[,,,,i] +
          DT * nonaids_excess_onart.ts
        
        artpop[,,,, i] <- artpop[,,,, i] + DT * gradART
        artpopdeaths[,,,, i] <- artpopdeaths[,,,, i] + DT * (artdeaths.ts + nonaids_excess_onart.ts)

        ## ART dropout
        ## remove proportion from all adult ART groups back to untreated pop

        art_dropout_ii <- fp$art_dropout[i]*colSums(artpop[1:2,,,,i])
        if (fp$art_dropout_recover_cd4) {
          art_dropout_ii[1,,] <- art_dropout_ii[1,,] +
            fp$art_dropout[i] * artpop[3:fp$ss$hTS,1,,,i]
          art_dropout_ii[-fp$ss$hDS,,] <- art_dropout_ii[-fp$ss$hDS,,] +
            fp$art_dropout[i] * artpop[3:fp$ss$hTS,-1,,,i]
        } else {
          art_dropout_ii <- art_dropout_ii +
            fp$art_dropout[i] * artpop[3:fp$ss$hTS,,,,i]
        }

        hivpop[,,,i] <- hivpop[,,,i] + DT * art_dropout_ii
        if(i >= fp$t_hts_start) { 
          diagnpop[,,,i] <- diagnpop[,,,i] + DT * art_dropout_ii
        }

        artpop[,,,,i] <- artpop[,,,,i] - DT * fp$art_dropout[i]*artpop[,,,,i]

        ## calculate number eligible for ART
        artcd4_percelig <- 1 - (1-rep(0:1, times=c(fp$artcd4elig_idx[i]-1, hDS - fp$artcd4elig_idx[i]+1))) *
          (1-rep(c(0, fp$who34percelig), c(2, hDS-2))) *
          (1-rep(fp$specpop_percelig[i], hDS))

        art15plus.elig <- sweep(hivpop[,h.age15plus.idx,,i], 1, artcd4_percelig, "*")

        ## calculate pregnant women
        if(fp$pw_artelig[i]) {
          births.dist <- sweep(fp$frr_cd4[,,i] * hivpop[,h.fert.idx,f.idx,i], 2,
                               births.by.h.age / (ctapply(pop[p.fert.idx, f.idx, hivn.idx, i], ag.idx[p.fert.idx], sum) + colSums(fp$frr_cd4[,,i] * hivpop[,h.fert.idx,f.idx,i]) + colSums(fp$frr_art[,,,i] * artpop[ ,,h.fert.idx,f.idx,i],,2)), "*")
          if(fp$artcd4elig_idx[i] > 1)
            art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] <- art15plus.elig[1:(fp$artcd4elig_idx[i]-1),h.fert.idx-min(h.age15plus.idx)+1,f.idx] + births.dist[1:(fp$artcd4elig_idx[i]-1),] # multiply by DT to account for proportion of annual births occurring during this time step
        }


        ## calculate number to initiate ART based on number or percentage

        artnum.ii <- c(0,0) # number on ART this ts
        if (fp$projection_period == "midyear" && DT*ii < 0.5) {
          for(g in 1:2) {
            if(!any(fp$art15plus_isperc[g,i-2:1])) {  # both number
              artnum.ii[g] <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
            } else if(all(fp$art15plus_isperc[g,i-2:1])) {  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-2:1] %*% c(1-(DT*ii+0.5), DT*ii+0.5))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            } else if(!fp$art15plus_isperc[g,i-2] & fp$art15plus_isperc[g,i-1]) { # transition number to percentage
              curr_coverage <- sum(artpop[,,h.age15plus.idx,g,i]) / (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i-1] - curr_coverage) * DT/(0.5-DT*(ii-1))
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            }
          }
        } else {
          for(g in 1:2){

            art_interp_w <- DT*ii
             if (fp$projection_period == "midyear") {
               art_interp_w <- art_interp_w - 0.5
             }
            
            if(!any(fp$art15plus_isperc[g,i-1:0])) {  # both number
              artnum.ii[g] <- c(fp$art15plus_num[g,i-1:0] %*% c(1-art_interp_w, art_interp_w))
            } else if(all(fp$art15plus_isperc[g,i-1:0])) {  # both percentage
              artcov.ii <- c(fp$art15plus_num[g,i-1:0] %*% c(1-art_interp_w, art_interp_w))                
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            } else if(!fp$art15plus_isperc[g,i-1] & fp$art15plus_isperc[g,i]) {  # transition number to percentage
              curr_coverage <- sum(artpop[,,h.age15plus.idx,g,i]) / (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
              artcov.ii <- curr_coverage + (fp$art15plus_num[g,i] - curr_coverage) * DT/(1.0 - art_interp_w + DT)
              artnum.ii[g] <- artcov.ii * (sum(art15plus.elig[,,g]) + sum(artpop[,,h.age15plus.idx,g,i]))
            }
          }
        }

        art15plus.inits <- pmax(artnum.ii - colSums(artpop[,,h.age15plus.idx,,i],,3), 0)

        ## calculate ART initiation distribution
        if(!fp$med_cd4init_input[i]) {

          if(fp$art_alloc_method == 4L) { ## by lowest CD4

            ## Calculate proportion to be initiated in each CD4 category
            artinit <- array(0, dim(art15plus.elig))
            remain_artalloc <- art15plus.inits
            for(m in hDS:1) {
              elig_hm <- colSums(art15plus.elig[m,,])
              init_prop <- ifelse(elig_hm == 0, elig_hm, pmin(1.0, remain_artalloc / elig_hm, na.rm=TRUE))
              artinit[m , , ] <- sweep(art15plus.elig[m,,], 2, init_prop, "*")
              remain_artalloc <- remain_artalloc - init_prop * elig_hm
            }

          } else {

            ## ## Old EPP-ASM implementation
            ##
            ## Applied 'expected mortality' weight within each cd4-age mortality strata.
            ## Different from Spectrum.
            ##             
            ## expect.mort.weight <- sweep(fp$cd4_mort[, h.age15plus.idx,], 3,
            ##                            colSums(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,],,2), "/")       ## 
            ## artinit.weight <- sweep(fp$art_alloc_mxweight * expect.mort.weight, 3, (1 - fp$art_alloc_mxweight)/colSums(art15plus.elig,,2), "+")
            ## artinit <- pmin(sweep(artinit.weight * art15plus.elig, 3, art15plus.inits, "*"),
            ##                 art15plus.elig, na.rm=TRUE)

            ## Spectrum ART initiation is 2-step process
            ## 1. Allocate by CD4 category (weighted by 'eligible' and 'expected mortality')
            ## 2. Allocate by age groups (weighted only by eligibility)

            ## First step: allocate initiation by CD4 category (_hm)
            artelig_hm <- apply(art15plus.elig, c(1, 3), sum)
            expected_deaths_hm <- apply(art15plus.elig * fp$cd4_mort[, h.age15plus.idx,], c(1, 3), sum)
            expected.mort.weight_hm <- sweep(expected_deaths_hm, 2, colSums(expected_deaths_hm), "/")
            artelig.weight_hm <- sweep(artelig_hm, 2, colSums(artelig_hm), "/")
            
            artinit.weight_hm <- fp$art_alloc_mxweight * expected.mort.weight_hm +
              (1.0 - fp$art_alloc_mxweight) * artelig.weight_hm

            artinit_hm <- sweep(artinit.weight_hm, 2, art15plus.inits, "*")

            ## Second step: within each CD4 category, allocate initiation
            ## proportionally by age
            
            ## Proportion initiating in each sex x CD4 category
            artinit_prob <- artinit_hm / artelig_hm
            artinit <- sweep(art15plus.elig, c(1, 3), artinit_prob, "*")

            artinit <- pmin(artinit, art15plus.elig, na.rm=TRUE)
          }

        } else {

          CD4_LOW_LIM <- c(500, 350, 250, 200, 100, 50, 0)
          CD4_UPP_LIM <- c(1000, 500, 350, 250, 200, 100, 50)

          medcd4_idx <- fp$med_cd4init_cat[i]

          medcat_propbelow <- (fp$median_cd4init[i] - CD4_LOW_LIM[medcd4_idx]) / (CD4_UPP_LIM[medcd4_idx] - CD4_LOW_LIM[medcd4_idx])

          elig_below <- colSums(art15plus.elig[medcd4_idx,,,drop=FALSE],,2) * medcat_propbelow
          if(medcd4_idx < hDS)
            elig_below <- elig_below + colSums(art15plus.elig[(medcd4_idx+1):hDS,,,drop=FALSE],,2)

          elig_above <- colSums(art15plus.elig[medcd4_idx,,,drop=FALSE],,2) * (1.0-medcat_propbelow)
          if(medcd4_idx > 1)
            elig_above <- elig_above + colSums(art15plus.elig[1:(medcd4_idx-1),,,drop=FALSE],,2)

          initprob_below <- pmin(art15plus.inits * 0.5 / elig_below, 1.0, na.rm=TRUE)
          initprob_above <- pmin(art15plus.inits * 0.5 / elig_above, 1.0, na.rm=TRUE)
          initprob_medcat <- initprob_below * medcat_propbelow + initprob_above * (1-medcat_propbelow)

          artinit <- array(0, dim=c(hDS, hAG, NG))

          if(medcd4_idx < hDS)
            artinit[(medcd4_idx+1):hDS,,] <- sweep(art15plus.elig[(medcd4_idx+1):hDS,,,drop=FALSE], 3, initprob_below, "*")
          artinit[medcd4_idx,,] <- sweep(art15plus.elig[medcd4_idx,,,drop=FALSE], 3, initprob_medcat, "*")
          if(medcd4_idx > 0)
            artinit[1:(medcd4_idx-1),,] <- sweep(art15plus.elig[1:(medcd4_idx-1),,,drop=FALSE], 3, initprob_above, "*")
        }

        if(i >= fp$t_hts_start) {

          ## 'newdiagn' is the number of new diagnoses are in deficit in
          ## the diagnosed population to match ART initiations
          newdiagn <- pmax(artinit - diagnpop[,,,i], 0)
          diagn_surplus <- pmax(diagnpop[,,,i] - artinit, 0)

          ## To not artificially inflate awareness, we substitute back those
          ## late diagnoses from those in available CD4 cell count categories.
          frac_exc <- colSums(newdiagn) / colSums(diagn_surplus)

          ## We only substitute back if there is enough people
          frac_exc[!is.finite(frac_exc)] <- 0
          frac_exc[frac_exc > 1] <- 1
          to_put_back <- sweep(diagn_surplus, 2:3, frac_exc, "*")

          ## 'prop_testneg' calculates the proportion of the undiagnosed HIV+
          ## population who have previously tested negative
          prop_testneg <- testnegpop[ , , hivp.idx, i] / colSums(hivpop[,,,i] - diagnpop[,,,i])
          prop_testneg[is.na(prop_testneg) | prop_testneg > 1 | prop_testneg < 0] <- 0

          late_diagnoses[,,,i] <- late_diagnoses[,,,i] + newdiagn - to_put_back
          diagnoses[,,,i] <- diagnoses[,,,i] + newdiagn - to_put_back

          # Here, we remove from the diagnpop the artinitiation (minus the late diagnoses)
          diagnpop[,,,i] <- diagnpop[,,,i] - (artinit - newdiagn + to_put_back)

          ## here, the 'testnegpop' now becomes aware according to their relative proportion.
          newdiagn_ha <- colSums(newdiagn) - colSums(to_put_back)

          hivtests[ , , 3, i] <- hivtests[ , , 3, i] + (1 - prop_testneg) * newdiagn_ha
          hivtests[ , , 4, i] <- hivtests[ , , 4, i] + prop_testneg * newdiagn_ha
          testnegpop[ , , hivp.idx, i] <- testnegpop[ , , hivp.idx, i] - prop_testneg * newdiagn_ha
        }

        hivpop[, h.age15plus.idx,, i] <- hivpop[, h.age15plus.idx,, i] - artinit
        artpop[1,, h.age15plus.idx,, i] <- artpop[1,, h.age15plus.idx,, i] + artinit
        artinits[,,,i] <- artinits[,,,i] + artinit
      }

      ## Remove hivdeaths from pop
      calc.agdist <- function(x) {d <- x/rep(ctapply(x, ag.idx, sum), h.ag.span); d[is.na(d)] <- 0; d}
      hivdeaths_p.ts <- apply(DT*hivdeaths_hAG.ts, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)  # HIV deaths by single-year age

      pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] - hivdeaths_p.ts
      hivdeaths[,,i] <- hivdeaths[,,i] + hivdeaths_p.ts

      nonaids_excess_p.ts <- apply(nonaids_excess_hAG.ts, 2, rep, h.ag.span) * apply(pop[,,hivp.idx,i], 2, calc.agdist)  # Non-AIDS excess deaths by single-year age
      pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] - nonaids_excess_p.ts
      excessnonaidsdeaths[,,i] <- excessnonaidsdeaths[,,i] + nonaids_excess_p.ts      

    }
    # ---- End Disease Model ----

    ## ## Code for calculating new infections once per year to match prevalence (like Spectrum)
    ## ## incidence
    ## prev.i <- sum(pop[p.age15to49.idx,,2,i]) / sum(pop[p.age15to49.idx,,,i]) # prevalence age 15 to 49
    ## incrate15to49.i <- (fp$prev15to49[i] - prev.i)/(1-prev.i)

    ## Direct incidence input
     if(fp$eppmod %in% c("directincid", "directinfections")) {
       if(fp$eppmod == "directincid") {
         if(fp$incidpopage == 0L) # incidence for 15-49 population
           p.incidpop.idx <- p.age15to49.idx
         else if(fp$incidpopage == 1L) # incidence for 15+ population
           p.incidpop.idx <- p.age15plus.idx
         incrate.i <- fp$incidinput[i]
         
         sexinc <- incrate.i*c(1, fp$incrr_sex[i])*sum(pop[p.incidpop.idx,,hivn.idx,i-1])/(sum(pop[p.incidpop.idx,m.idx,hivn.idx,i-1]) + fp$incrr_sex[i]*sum(pop[p.incidpop.idx, f.idx,hivn.idx,i-1]))
         agesex.inc <- sweep(fp$incrr_age[,,i], 2, sexinc/(colSums(pop[p.incidpop.idx,,hivn.idx,i-1] * fp$incrr_age[p.incidpop.idx,,i])/colSums(pop[p.incidpop.idx,,hivn.idx,i-1])), "*")
         infections[,,i] <- agesex.inc * pop[,,hivn.idx,i-1]
      } else if(fp$eppmod == "directinfections") {
         infections[,,i] <- fp$infections[,,i]
       }

       infections_ha <- apply(infections[,,i], 2, ctapply, ag.idx, sum)

       if(i >= fp$t_hts_start) {
         hivn_pop_ha <- apply(pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum)
         testneg_infections_ha <- infections_ha * testnegpop[,,hivn.idx,i] / hivn_pop_ha
         testnegpop[ , , hivn.idx, i] <- testnegpop[ , , hivn.idx, i] - testneg_infections_ha
         testnegpop[ , , hivp.idx, i] <- testnegpop[ , , hivp.idx, i] + testneg_infections_ha
       }
     
       pop[,,hivn.idx,i] <- pop[,,hivn.idx,i] - infections[,,i]
       pop[,,hivp.idx,i] <- pop[,,hivp.idx,i] + infections[,,i]
     
       hivpop[,,,i] <- hivpop[,,,i] + sweep(fp$cd4_initdist, 2:3, infections_ha, "*")
     } # if(fp$eppmod %in% c("directincid", "directinfections"))

    if (fp$projection_period == "calendar") {
      
      ## net migration
      netmigsurv <- fp$netmigr[,,i]
      mr.prob <- 1+netmigsurv / rowSums(pop[,,,i],,2)
      hiv.mr.prob <- apply(mr.prob * pop[,,hivp.idx,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,hivp.idx,i], 2, ctapply, ag.idx, sum)
      hiv.mr.prob[is.nan(hiv.mr.prob)] <- 0

      if(i >= fp$t_hts_start) {
        hivn.mr.prob <- apply(mr.prob * pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum) /  apply(pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum)
        hivn.mr.prob[is.nan(hivn.mr.prob)] <- 0
      }
      
      pop[,,,i] <- sweep(pop[,,,i], 1:2, mr.prob, "*")
      
      hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hiv.mr.prob, "*")
      if(i >= fp$t_hts_start) {
        testnegpop[,, hivn.idx,i] <- testnegpop[,,hivn.idx,i] * hivn.mr.prob
        testnegpop[,, hivp.idx,i] <- testnegpop[,,hivp.idx,i] * hiv.mr.prob
        
        diagnpop[,,,i] <- sweep(diagnpop[,,,i], 2:3, hiv.mr.prob, "*")
      }
      if(i >= fp$tARTstart)
        artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, hiv.mr.prob, "*")
    }

    
    ## adjust HIV population to match target population size
    if(fp$popadjust) {
      
      hivpop_ha <- colSums(hivpop[,,,i])
      hivpop_ha_adj <- fp$target_hivpop_ha[,,i] / hivpop_ha
      hivpop_ha_adj[!is.finite(hivpop_ha_adj)] <- 1.0
      hivpop[,,,i] <- sweep(hivpop[,,,i], 2:3, hivpop_ha_adj, "*")
      
      if(i >= fp$t_hts_start) {
        
        hivn_ha <- apply(pop[,,hivn.idx,i], 2, ctapply, ag.idx, sum)
        hivn_ha_adj <- fp$target_hivn_ha[,,i] / hivn_ha
        hivn_ha_adj[!is.finite(hivn_ha_adj)] <- 1.0
        testnegpop[,,hivn.idx,i] <- testnegpop[,,hivn.idx,i] * hivn_ha_adj
        
        testnegpop[,,hivp.idx,i] <- testnegpop[,,hivp.idx,i] * hivpop_ha_adj
        diagnpop[,,,i] <- sweep(diagnpop[,,,i], 2:3, hivpop_ha_adj, "*")
      }
      
      if(i >= fp$tARTstart) {
        artpop_ha <- colSums(artpop[,,,,i], dims=2)
        artpop_ha_adj <- fp$target_artpop_ha[,,i] / artpop_ha
        artpop_ha_adj[!is.finite(artpop_ha_adj)] <- 1.0
        artpop[,,,,i] <- sweep(artpop[,,,,i], 3:4, artpop_ha_adj, "*")
      }

      pop[,,hivn.idx,i] <- fp$target_hivn_pop[,,i]
      pop[,,hivp.idx,i] <- fp$target_hivp_pop[,,i]
    }

    ## prevalence and incidence 15 to 49
    prev15to49[i] <- sum(pop[p.age15to49.idx,,hivp.idx,i]) / sum(pop[p.age15to49.idx,,,i])

    if (fp$projection_period == "calendar") {
      ## incidence: interpolated denominator
      incid15to49_denom <- 0.5 * (sum(pop[p.age15to49.idx,,hivn.idx,i-1]) + sum(pop[p.age15to49.idx,,hivn.idx,i]))      
    } else {
      incid15to49_denom <- sum(pop[p.age15to49.idx,,hivn.idx,i-1])
    }
    incid15to49[i] <- sum(infections[p.age15to49.idx,,i]) / incid15to49_denom

  }


  attr(pop, "prev15to49") <- prev15to49
  attr(pop, "incid15to49") <- incid15to49
  attr(pop, "sexinc") <- sexinc15to49out
  attr(pop, "hivpop") <- hivpop
  attr(pop, "diagnpop") <- diagnpop
  attr(pop, "artpop") <- artpop

  attr(pop, "testnegpop") <- testnegpop

  attr(pop, "infections") <- infections
  attr(pop, "hivdeaths") <- hivdeaths
  attr(pop, "natdeaths") <- natdeaths

  attr(pop, "hivpopdeaths") <- hivpopdeaths
  attr(pop, "artpopdeaths") <- artpopdeaths

  attr(pop, "excessnonaidsdeaths") <- excessnonaidsdeaths

  attr(pop, "aidsdeaths_noart") <- aidsdeaths_noart
  attr(pop, "natdeaths_noart") <- natdeaths_noart
  attr(pop, "excessnonaidsdeaths_noart") <- excessnonaidsdeaths_noart
  
  attr(pop, "aidsdeaths_art") <- aidsdeaths_art
  attr(pop, "natdeaths_art") <- natdeaths_art
  attr(pop, "excessnonaidsdeaths_art") <- excessnonaidsdeaths_art
  

  attr(pop, "hivtests") <- hivtests
  attr(pop, "diagnoses") <- diagnoses
  attr(pop, "late_diagnoses") <- late_diagnoses
  attr(pop, "artinits") <- artinits

  attr(pop, "popadjust") <- popadj.prob

  if(fp$eppmod != "directincid") {
    attr(pop, "incrate15to49_ts") <- incrate15to49.ts.out
    attr(pop, "prev15to49_ts") <- prev15to49.ts.out
  }

  attr(pop, "entrantprev") <- entrant_prev_out
  attr(pop, "hivp_entrants") <- hivp_entrants_out
  class(pop) <- "spec"
  return(pop)
}
