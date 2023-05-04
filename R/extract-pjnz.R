
first90_read_csv_character <- function(file) {
  v <- vroom::vroom(file, delim = ",",
                    col_types = vroom::cols(.default = vroom::col_character()),
                    .name_repair = "minimal", progress = FALSE)
  v <- as.data.frame(v)
  v
}

#' Extract outputs from PJNZ needed for first90 model
#'
#' @param pjnz filepath to PJNZ file
#' @param dp_file filepath to a .DP file
#' @param pjn_file filepath to a .PJN file
#'
#' @return a list
#'
#' @export
extract_pjnz <- function(pjnz = NULL, dp_file= NULL, pjn_file = NULL){

  if(!is.null(pjnz) && is.null(dp_file) && is.null(pjn_file)) {
    dp_file <- file_in_zip(pjnz, ".DP$")
    pjn_file <- file_in_zip(pjnz, ".PJN$")
  } else if(!is.null(pjnz) && (!is.null(dp_file) || !is.null(pjn_file))) 
    stop("Must provide either a Spectrum .PJNZ file or a .DP and .PJN file. Do not provide both.")
  else if(is.null(pjnz) && (is.null(dp_file) || is.null(pjn_file)))
    stop("Must provide either a Spectrum .PJNZ file or both a .DP and .PJN file")

  dp <- first90_read_csv_character(dp_file)
  pjn <- first90_read_csv_character(pjn_file)


  year_start <- as.integer(dpsub(dp, "<FirstYear MV2>",2,4))
  year_end <- as.integer(dpsub(dp, "<FinalYear MV2>",2,4))
  proj_years <- year_start:year_end

  v <- list()

  ## Information
  v$country <- get_pjn_country(pjn)
  v$region <- get_pjn_region(pjn)
  v$iso3 <- get_pjn_iso3(pjn)

  v$projection_name <- pjn[which(pjn[,1] == "<Projection Name>")+2, 4]
  v$spectrum_version <- pjn[which(pjn[,1] == "<Projection General>")+4, 4]

  ## Replace comma decimal separator save on Francophone locale computers
  ## (e.g. replace 6,14 by 6.14
  v$spectrum_version <- sub("^([0-9]+),(.*)$", "\\1.\\2", v$spectrum_version)


  v$valid_date <- dpsub(dp, "<ValidDate MV>", 2, 3)

  ## Spectrum custom population used
  if (exists_dptag(dp, "<RegionalAdjustPopCBState MV>")) {
    v$popadjust <- dpsub(dp, "<RegionalAdjustPopCBState MV>", 2, 4)
    v$popadjust <- as.integer(v$popadjust) == 1
  } else {
    v$popadjust <- FALSE
  }
    
  
  ## Projection inputs
  v$yr_start <- year_start
  v$yr_end <- year_end
  
  ## Demographic inputs
  v$totpop <- get_dp_totpop(dp, proj_years)
  v$Sx <- get_dp_Sx(dp, proj_years)
  v$asfr <- calc_asfr(get_dp_tfr(dp, proj_years),
                      get_dp_asfd(dp, proj_years))
  v$srb <- get_dp_srb(dp, proj_years)
  v$births <- get_dp_births(dp, proj_years)
  v$netmigr <- calc_netmigr(get_dp_totnetmig(dp, proj_years),
                            get_dp_netmigagedist(dp, proj_years),
                            v$Sx)

  ## Epidemic results
  v$hivpop <- get_dp_hivpop(dp, proj_years)
  v$infections <- get_dp_infections(dp, proj_years)
  v$artpop <- get_dp_artpop(dp, proj_years)

  ## AIM parameters

  dn <- list(sex = c("male", "female"), year = proj_years)
  col_idx <- 3+seq_along(proj_years)

  v$relinfectART <- 1.0 - as.numeric(dpsub(dp, "<AdultInfectReduc MV>", 2, 4))
  v$incidpopage <- as.integer(dpsub(dp, "<EPPPopulationAges MV>", 2, 4))  # Adults 15-49 = 0; Adults 15+ = 1

  v$incrr_sex <- get_dp_incrr_sex(dp, proj_years)
  v$incrr_age <- get_dp_incrr_age(dp, proj_years)

  v_frr <- get_dp_frr(dp, proj_years)
  v[names(v_frr)] <- v_frr

  ## # Natural history

  v$cd4_initdist <- get_dp_cd4_initdist(dp)
  v$cd4_prog <- get_dp_cd4_prog(dp)
  v$cd4_mort <- get_dp_cd4_mort(dp)
  v$art_mort <- get_dp_art_mort(dp)
  v$artmx_timerr <- get_dp_artmx_timerr(dp, proj_years)
  
  ## # ART programme data
  v$art15plus_numperc <- array(as.numeric(unlist(dpsub(dp, "<HAARTBySexPerNum MV>", 4:5, col_idx))), lengths(dn), dn)
  v$art15plus_num <- array(as.numeric(unlist(dpsub(dp, "<HAARTBySex MV>", 4:5, col_idx))), lengths(dn), dn)

  male_15plus_needart <- dpsub(dp, "<NeedARTDec31 MV>", 4:17*3 + 3, col_idx)
  male_15plus_needart <- vapply(lapply(male_15plus_needart, as.numeric), sum, numeric(1))

  female_15plus_needart <- dpsub(dp, "<NeedARTDec31 MV>", 4:17*3 + 4, col_idx)
  female_15plus_needart <- vapply(lapply(female_15plus_needart, as.numeric), sum, numeric(1))

  v$art15plus_needart <- rbind(male_15plus_needart, female_15plus_needart)
  dimnames(v$art15plus_needart) <- dn

  ## # Adult on ART adjustment factor
  ## 
  ## * Implemented from around Spectrum 6.2 (a few versions before)
  ## * Allows user to specify scalar to reduce number on ART in each year ("<AdultARTAdjFactor>")
  ## * Enabled / disabled by checkbox flag ("<AdultARTAdjFactorFlag>")
  ## * Scaling factor only applies to number inputs, not percentages (John Stover email, 20 Feb 2023)
  ##   -> Even if scaling factor specified in a year with percentage input, ignore it.

  if (exists_dptag(dp, "<AdultARTAdjFactorFlag>") &&
        dpsub(dp, "<AdultARTAdjFactorFlag>", 2, 4) == 1) {

    adult_artadj_factor <- array(as.numeric(unlist(dpsub(dp, "<AdultARTAdjFactor>", 3:4, col_idx))), lengths(dn), dn)

    ## Only apply if is number (! is percentage)
    adult_artadj_factor <- adult_artadj_factor ^ as.numeric(!v$art15plus_numperc)

    v$art15plus_num <- v$art15plus_num * adult_artadj_factor
  }

  v$art15plus_eligthresh <- setNames(as.numeric(dpsub(dp, "<CD4ThreshHoldAdults MV>", 2, col_idx)), proj_years)

  artelig_specpop <- setNames(dpsub(dp, "<PopsEligTreat MV>", 3:9, 2:6),
                              c("description", "pop", "elig", "percent", "year"))
  artelig_specpop$pop <- c("PW", "TBHIV", "DC", "FSW", "MSM", "IDU", "OTHER")
  artelig_specpop$elig <- as.logical(as.integer(artelig_specpop$elig))
  artelig_specpop$percent <- as.numeric(artelig_specpop$percent)/100
  artelig_specpop$year <- as.integer(artelig_specpop$year)
  artelig_specpop$idx <- match(as.integer(artelig_specpop$year), proj_years)
  rownames(artelig_specpop) <- artelig_specpop$pop
  v$artelig_specpop <- artelig_specpop

  v$median_cd4init <- setNames(as.numeric(dpsub(dp, "<MedCD4CountInit MV>", 2, col_idx)), proj_years)
  v$art_dropout <- setNames(as.numeric(dpsub(dp, "<PercLostFollowup MV>", 2, col_idx)), proj_years)

  v$art_alloc_method <- get_dp_art_alloc_method(dp)
  v$art_prop_alloc <- get_dp_art_prop_alloc(dp)
  v$scale_cd4_mort <- get_dp_scale_cd4_mort(dp)

  v$age14hivpop <- get_dp_age14hivpop(dp, proj_years)

  v
}

#' Combine PJNZ inputs
#'
#' @param lst a list of inputs, each returned from extract_pjnz()
#'
#' @export

combine_inputs <- function(lst){

  stopifnot(vapply(lst, "[[", integer(1), "yr_start") == lst[[1]]$yr_start,
            vapply(lst, "[[", integer(1), "yr_end") == lst[[1]]$yr_end)

  ## # Convert any ART percentage inputs to counts
  lst <- lapply(lst, 
                function(x){
                  isperc <- x$art15plus_numperc == 1
                  x$art15plus_num[isperc] <- x$art15plus_needart[isperc] * x$art15plus_num[isperc] / 100
                  x$art15plus_numperc[] <- 0
                  x
                })

  if(length(lst) == 1)
    return(lst[[1]])

  
  v <- list()
  
  ## Inputs provided as a concatenation of all regions
  nm <- c("country", "region", "iso3", "projection_name", "spectrum_version", "valid_date")
  vv <- do.call(Map, c(f = c, lapply(lst, "[", nm)))
  v[names(vv)] <- vv
  
  ## Inputs assumed not to change across regions, take first in list
  ## Note: FRR parameters don't matter for first90 model at present. Might
  ##       make a difference in future if start modelling ANC testing separately,
  ##       in which case more care needs to be taken.
  nm <- c("yr_start", "yr_end", "relinfectART",
          "fert_rat", "cd4fert_rat", "frr_art6mos", "frr_scalar",  
          "cd4_initdist", "cd4_prog", "cd4_mort", "art_mort",
          "art15plus_eligthresh", "artelig_specpop", "art15plus_numperc",
          "art_alloc_method", "art_prop_alloc", "scale_cd4_mort")
          
  v[nm] <- lst[[1]][nm]
    

  ## Inputs calculated as functions of all inputs

  ## # Sums over regions
  nm <- c("totpop", "births", "netmigr",  "hivpop", "infections", "artpop",
          "art15plus_num", "art15plus_needart", "age14hivpop")
  vv <- do.call(Map, c(f = list, lapply(lst, "[", nm)))
  vv <- lapply(vv, Reduce, f="+")
  v[names(vv)] <- vv


  ## # Popadjust: use if any region used
  v$popadjust <- any(vapply(lst, "[[", logical(1), "popadjust"))

  ## # Demographic inputs

  totpop <- lapply(lst, "[[", "totpop")
  v$Sx <- 1 - Reduce("+", Map("*", lapply(lst, function(x) 1 - x$Sx), totpop)) / Reduce("+", totpop)

  fert_pop <- lapply(lst, function(x) x$totpop[16:50, "female", ])
  v$asfr <- Reduce("+", Map("*", lapply(lst, "[[", "asfr"), fert_pop)) / Reduce("+", fert_pop)
  
  m_births <- Map("*", lapply(lst, "[[", "births"),
                  lapply(lapply(lst, "[[", "srb"), function(x) x/(1+x)))
  f_births <- Map("*", lapply(lst, "[[", "births"),
                  lapply(lapply(lst, "[[", "srb"), function(x) 1/(1+x)))
  v$srb <- Reduce("+", m_births) / Reduce("+", f_births)
  
  ## # Incidence rate ratios

  inf15to49 <- colSums(v$infections[16:50,,])
  susc15to49 <- colSums(v$totpop[16:50,,]) - colSums(v$hivpop[16:50,,])
  incid15to49 <- inf15to49 / cbind(0, susc15to49[ , -ncol(susc15to49)])
  v$incrr_sex <- incid15to49[2, ] / incid15to49[1, ]
  v$incrr_sex[!is.finite(v$incrr_sex)] <- 1.0

  agegr <- rep(factor(agegr_labels(), agegr_labels()), times=c(rep(5, 16), 1))
  inf5yr <- apply(v$infections, 2:3, tapply, agegr, sum)
  susc5yr <- apply(v$totpop, 2:3, tapply, agegr, sum) - apply(v$hivpop, 2:3, tapply, agegr, sum) + inf5yr
  incid5yr <- inf5yr / susc5yr
  v$incrr_age <- sweep(incid5yr, 2:3, incid5yr["25-29",,], "/")
  v$incrr_age[!is.finite(v$incrr_age)] <- 1.0

  
  ## # artmx_timerr, median CD4 and ART dropout: weighted mean by number on ART
  ## # Note: median CD4 could go really wrong if entered for some regions and not
  ##         for others. Set NA if not entered and use weighted mean for only
  ##         regions where it is entered.

  art15plus <- lapply(lapply(lst, "[[", "art15plus_num") , colSums)

  art15plus_aggr <- Reduce("+", art15plus)

  artmx_timerr <- lapply(lst, "[[", "artmx_timerr")
  artmx_timerr <- sweep(Reduce("+", Map("sweep", artmx_timerr, 2, art15plus, "*")), 2, Reduce("+", art15plus), "/")
  artmx_timerr[is.na(artmx_timerr) | is.infinite(artmx_timerr)] <- 1.0
  v$artmx_timerr <- artmx_timerr

  medcd4 <- lapply(lst, "[[", "median_cd4init")
  medcd4_aggr <- Reduce("+", Map("*", medcd4, art15plus)) / Reduce("+", Map("*", lapply(medcd4, "!=", 0), art15plus))
  medcd4_aggr[is.na(medcd4_aggr) | is.infinite(medcd4_aggr)] <- 0
  v$median_cd4init <- medcd4_aggr

  art_dropout <- lapply(lst, "[[", "art_dropout")
  art_dropout <- Reduce("+", Map("*", art_dropout, art15plus)) / Reduce("+", art15plus)
  art_dropout[is.na(art_dropout) | is.infinite(art_dropout)] <- 0
  v$art_dropout <- art_dropout


  v <- v[names(lst[[1]])] # reorder to match input list
  v
}


#' Internal helper functions
#' 
exists_dptag <- function(dp, tag, tagcol=1){tag %in% dp[,tagcol]}

dpsub <- function(dp, tag, rows, cols, tagcol=1){
  dp[which(dp[,tagcol]==tag)+rows, cols]
}

specpop_dimnames <- function(proj_years){
  list(age = 0:80,
       sex = c("male", "female"),
       year = proj_years )
}

agegr_labels <- function(){
  c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", 
    "35-39", "40-44", "45-49", "50-54", "55-59", "60-64", "65-69", 
    "70-74", "75-79", "80+")
}

NG <- 2
AG <- 17
DS <- 7
TS <- 3
fAG <- 7  # fertile age groups


#' Extract arrays from Spectrum DP file
#'
#' @export
get_dp_totpop <- function(dp, proj_years){

  dn <- specpop_dimnames(proj_years)  
  col_idx <- 3 + seq_along(proj_years)
  
  if(exists_dptag(dp, "<BigPop3>"))
    totpop <- dpsub(dp, "<BigPop3>", 2:163, col_idx)
  else if(exists_dptag(dp, "<BigPop MV>"))
    totpop <- dpsub(dp, "<BigPop MV>", 3:164, col_idx)
  else if(exists_dptag(dp, "<BigPop MV2>"))
    totpop <- dpsub(dp, "<BigPop MV2>", c(3+0:80, 246+0:80), col_idx)
  else if(exists_dptag(dp, "<BigPop MV3>"))
    totpop <- dpsub(dp, "<BigPop MV3>", 3:164, col_idx)

  totpop <- sapply(totpop, as.numeric)
  totpop <- array(totpop, lengths(dn), dn)
  
  totpop
}

get_dp_hivpop <- function(dp, proj_years){

  dn <- specpop_dimnames(proj_years)
  col_idx <- 3 + seq_along(proj_years)

  if(exists_dptag(dp, "<HIVBySingleAge MV>"))
    hivpop <- dpsub(dp, "<HIVBySingleAge MV>", c(3:83, 85:165), col_idx)
  else if(exists_dptag(dp, "<HIVBySingleAge MV2>"))
    hivpop <- dpsub(dp, "<HIVBySingleAge MV2>", 3:164, col_idx)
  
  hivpop <- array(as.numeric(unlist(hivpop)), lengths(dn), dn)
  
  hivpop
}

get_dp_artpop <- function(dp, proj_years){

  dn <- specpop_dimnames(proj_years)
  col_idx <- 3 + seq_along(proj_years)

  artpop_m <- dpsub(dp, "<OnARTBySingleAge MV>", 0:80*3 + 3, col_idx)
  artpop_f <- dpsub(dp, "<OnARTBySingleAge MV>", 0:80*3 + 4, col_idx)
  artpop <- rbind(artpop_m, artpop_f)

  artpop <- array(as.numeric(unlist(artpop)), lengths(dn), dn)
  
  artpop
}

get_dp_infections5yr <- function(dp, proj_years){

  dn <- list(agegr = agegr_labels(), sex = c("male", "female"), year = proj_years)

  v <- dpsub(dp, "<NewInfections MV>", 2:55, 3+seq_along(proj_years))
  v <- rbind(v[1:17*3 + 2,], v[1:17*3 + 3,])
  v <- array(as.numeric(unlist(v)), lengths(dn), dn)

  v
}

get_dp_infections <- function(dp, proj_years){

  dn <- specpop_dimnames(proj_years)
  v <- dpsub(dp, "<NewInfectionsBySingleAge MV>", 2:244, 3+seq_along(proj_years))
  v <- rbind(v[0:80 * 3 + 2, ], v[0:80 * 3 + 3, ])
  v <- array(as.numeric(unlist(v)), lengths(dn), dn)
  
  v
}
      

#' Get mortality probability (Sx) from Spectrum DP file
#' 
#' @details
#' This function extracts Sx for ages 0:79 and 80+. Spectrum calculates a
#' separate Sx for age 80. The population projection model in EPP-ASM needs
#' to be updated to handle this.
#' 
get_dp_Sx <- function(dp, proj_years){

  dn <- specpop_dimnames(proj_years)
  col_idx <- 3 + seq_along(proj_years)

  if(exists_dptag(dp, "<SurvRate MV>"))
    Sx <- dpsub(dp, "<SurvRate MV2>", 3+c(0:79,81, 83+0:79, 83+81), col_idx)
  else if(exists_dptag(dp, "<SurvRate MV2>"))
    Sx <- dpsub(dp, "<SurvRate MV2>", 3+c(0:79, 81, 82+0:79, 82+81), col_idx)

  Sx <- array(as.numeric(unlist(Sx)), lengths(dn), dn)
  Sx
}

get_dp_totnetmig <- function(dp, proj_years){

  col_idx <- 3 + seq_along(proj_years)

  if(exists_dptag(dp, "<MigrRate MV>"))
    totnetmig <- sapply(dpsub("<MigrRate MV>", c(5, 8), col_idx), as.numeric)
  else if(exists_dptag(dp, "<MigrRate MV2>"))
    totnetmig <- sapply(dpsub(dp, "<MigrRate MV2>", c(4, 6), col_idx), as.numeric)

  dimnames(totnetmig) <- list(sex = c("male", "female"), year = proj_years)

  totnetmig
}

get_dp_netmigagedist <- function(dp, proj_years){

  col_idx <- 3 + seq_along(proj_years)

  if(exists_dptag(dp, "<MigrAgeDist MV>"))
    netmigagedist <- dpsub(dp, "<MigrAgeDist MV>", 6 + c(1:17*2, 37+1:17*2), col_idx)
  else if(exists_dptag(dp, "<MigrAgeDist MV2>"))
    netmigagedist <- dpsub(dp, "<MigrAgeDist MV2>", 2+1:34, col_idx)

  dn <- list(agegr = agegr_labels(), sex = c("male", "female"), year = proj_years)

  netmigagedist <- array(as.numeric(unlist(netmigagedist))/100, lengths(dn), dn)
  netmigagedist
}



#' Get sex ratio at birth
#'
#' @export
get_dp_srb <- function(dp, proj_years){

  srb <- dpsub(dp, "<SexBirthRatio MV>", 2, 3+seq_along(proj_years))
  srb <- as.numeric(srb)
  names(srb) <- proj_years

  srb
}

#' Get age-specific fertility rate by single-year
#'
#' @export
get_dp_tfr <- function(dp, proj_years){

  tfr <- dpsub(dp, "<TFR MV>", 2, 3+seq_along(proj_years))
  tfr <- as.numeric(tfr)
  names(tfr) <- proj_years
  tfr
}

get_dp_asfd <- function(dp, proj_years){
  
  asfd <- dpsub(dp, "<ASFR MV>", 3:9, 3+seq_along(proj_years))
  asfd <- as.numeric(unlist(asfd))/100
  asfd <- array(asfd,
                c(7, length(proj_years)),
                list(agegr = agegr_labels()[4:10],
                     year = proj_years))

  asfd
}

get_dp_births <- function(dp, proj_years){
  births <- dpsub(dp, "<Births MV>", 2, 3+seq_along(proj_years))
  births <- as.numeric(births)
  names(births) <- proj_years

  births
}

#' Extract AIM module parameters
#' 

get_dp_frr <- function(dp, proj_years){

  col_idx <- 3 + seq_along(proj_years)
  
  if(exists_dptag(dp, "<HIVTFR MV>")){
    fert_rat <- vapply(dpsub(dp, "<HIVTFR MV>", 2:8, col_idx), as.numeric, numeric(fAG))
    dimnames(fert_rat) <- list(agegr = agegr_labels()[4:10], year = proj_years)
  } else if(exists_dptag(dp, "<HIVTFR MV2>")) {
    ## this version of Spectrum stratified fertility reduction by 15-17, 18-19, 20-24, ...
    fert_rat <- vapply(dpsub(dp, "<HIVTFR MV2>", 2:7, col_idx), as.numeric, numeric(6))
    dimnames(fert_rat) <- list(agegr = c("15-17", "18-19", "20-24", "25-29", "30-34", "35-49"),
                               year = proj_years)
  } else if(exists_dptag(dp, "<HIVTFR MV3>")) {
    fert_rat <- vapply(dpsub(dp, "<HIVTFR MV3>", 2:8, col_idx), as.numeric, numeric(fAG))
    dimnames(fert_rat) <- list(agegr = agegr_labels()[4:10], year = proj_years)
  } else if(exists_dptag(dp, "<HIVTFR MV4>")) {
    fert_rat <- vapply(dpsub(dp, "<HIVTFR MV4>", 2:8, col_idx), as.numeric, numeric(fAG))
    dimnames(fert_rat) <- list(agegr = agegr_labels()[4:10], year = proj_years)
  }


  if(exists_dptag(dp, "<FertCD4Discount MV>"))
    cd4fert_rat <- as.numeric(dpsub(dp, "<FertCD4Discount MV>", 2, 4+seq_len(DS)))
  else
    cd4fert_rat <- rep(1.0, DS)

  if(exists_dptag(dp, "<RatioWomenOnART MV>"))
    frr_art6mos <- rep(as.numeric(dpsub(dp, "<RatioWomenOnART MV>", 2, 4)), fAG)
  else if(exists_dptag(dp, "<RatioWomenOnART MV2>"))
    frr_art6mos <- as.numeric(dpsub(dp, "<RatioWomenOnART MV2>", 2, 3 + seq_len(fAG)))
  else
    frr_art6mos <- rep(1.0, fAG)
  names(frr_art6mos) <- agegr_labels()[4:10]

  if(exists_dptag(dp, "<FRRbyLocation MV>"))
    frr_scalar <- as.numeric(dpsub(dp, "<FRRbyLocation MV>", 2, 4))
  else
    frr_scalar <- 1.0


  list(fert_rat = fert_rat,
       cd4fert_rat = cd4fert_rat,
       frr_art6mos = frr_art6mos,
       frr_scalar = frr_scalar)
}

get_dp_incrr_sex <- function(dp, proj_years){
  incrr_sex <- as.numeric(dpsub(dp, "<HIVSexRatio MV>", 2, 3 + seq_along(proj_years)))
  names(incrr_sex) <- proj_years
  incrr_sex
}

get_dp_incrr_age <- function(dp, proj_years){

  col_idx <- 3 + seq_along(proj_years)
  
  if(exists_dptag(dp, "<DistOfHIV MV>"))
    incrr_age <- dpsub(dp, "<DistOfHIV MV>", c(4:20, 22:38), col_idx)
  else if(exists_dptag(dp, "<DistOfHIV MV2>"))
    incrr_age <- dpsub(dp, "<DistOfHIV MV2>", 3:36, col_idx)

  dn <- list(agegr = agegr_labels(), sex = c("male", "female"), year = proj_years)
  incrr_age <- array(as.numeric(unlist(incrr_age)), lengths(dn), dn)

  incrr_age
}

get_dp_cd4_initdist<- function(dp){

  cd4_initdist <- array(NA, c(DS, 4, NG),
                        list(cd4stage=1:DS,
                             agecat=c("15-24", "25-34", "35-44", "45+"),
                             sex=c("male", "female")))

  cd4_initdist[,,"male"] <- array(as.numeric(dpsub(dp, "<AdultDistNewInfectionsCD4 MV>", 3, 4:31))/100, c(DS, 4))
  cd4_initdist[,,"female"] <- array(as.numeric(dpsub(dp, "<AdultDistNewInfectionsCD4 MV>", 4, 4:31))/100, c(DS, 4))

  cd4_initdist
}

get_dp_cd4_prog<- function(dp){

  cd4_prog <- array(NA, c(DS-1, 4, NG),
                    list(cd4stage=1:(DS-1),
                         agecat=c("15-24", "25-34", "35-44", "45+"),
                         sex=c("male", "female")))
  
  cd4_prog[,,"male"] <- array(as.numeric(dpsub(dp, "<AdultAnnRateProgressLowerCD4 MV>", 3, 4:31)), c(DS, 4))[1:(DS-1),]
  cd4_prog[,,"female"] <- array(as.numeric(dpsub(dp, "<AdultAnnRateProgressLowerCD4 MV>", 4, 4:31)), c(DS, 4))[1:(DS-1),]

  cd4_prog
}

get_dp_cd4_mort<- function(dp){

  cd4_mort <- array(NA, c(DS, 4, NG),
                    list(cd4stage=1:DS,
                         agecat=c("15-24", "25-34", "35-44", "45+"),
                         sex=c("male", "female")))
  
  cd4_mort[,,"male"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4NoART MV>", 3, 4:31)), c(DS, 4))
  cd4_mort[,,"female"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4NoART MV>", 4, 4:31)), c(DS, 4))

  cd4_mort
}

get_dp_art_mort <- function(dp){

  art_mort <- array(NA, c(TS, DS, 4, NG),
                    list(artdur=c("ART0MOS", "ART6MOS", "ART1YR"),
                         cd4stage=1:DS,
                         agecat=c("15-24", "25-34", "35-44", "45+"),
                         sex=c("male", "female")))
  
  if(exists_dptag(dp, "<AdultMortByCD4WithART0to6 MV>")) {
      art_mort[1,,,"male"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART0to6 MV>", 3, 4:31)), c(DS, 4))
      art_mort[1,,,"female"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART0to6 MV>", 4, 4:31)), c(DS, 4))
      art_mort[2,,,"male"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART7to12 MV>", 3, 4:31)), c(DS, 4))
      art_mort[2,,,"female"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART7to12 MV>", 4, 4:31)), c(DS, 4))
      art_mort[3,,,"male"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithARTGt12 MV>", 3, 4:31)), c(DS, 4))
      art_mort[3,,,"female"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithARTGt12 MV>", 4, 4:31)), c(DS, 4))
  } else if(exists_dptag(dp, "<AdultMortByCD4WithART0to6 MV2>")) {
      art_mort[1,,,"male"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART0to6 MV2>", 2, 4:31)), c(DS, 4))
      art_mort[1,,,"female"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART0to6 MV2>", 3, 4:31)), c(DS, 4))
      art_mort[2,,,"male"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART7to12 MV2>", 2, 4:31)), c(DS, 4))
      art_mort[2,,,"female"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithART7to12 MV2>", 3, 4:31)), c(DS, 4))
      art_mort[3,,,"male"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithARTGt12 MV2>", 2, 4:31)), c(DS, 4))
      art_mort[3,,,"female"] <- array(as.numeric(dpsub(dp, "<AdultMortByCD4WithARTGt12 MV2>", 3, 4:31)), c(DS, 4))
  }

  art_mort
}

get_dp_age14hivpop <- function(dp, proj_years){

  PAED_DS <- 6 # number of paediatric stages of infection
  age14hivpop <- as.numeric(unlist(dpsub(dp, "<ChAged14ByCD4Cat MV>", 1+1:(NG*PAED_DS*(4+TS)), 3+seq_along(proj_years))))
  
  age14hivpop <- array(age14hivpop, c(4+TS, PAED_DS, NG, length(proj_years)),
                       list(artstage=c("PERINAT", "BF0MOS", "BF6MOS", "BF1YR", "ART0MOS", "ART6MOS", "ART1YR"),
                            cd4cat = c("CD4_1000", "CD4_750", "CD4_500", "CD4_350", "CD4_200", "CD4_0"),
                            sex = c("male", "female"), year = proj_years))

  age14hivpop
}

get_dp_art_alloc_method <- function(dp){

  if(exists_dptag(dp, "<NewARTPatAllocationMethod MV2>"))
    art_alloc_method <- as.integer(dpsub(dp, "<NewARTPatAllocationMethod MV2>", 2, 4))
  else
    art_alloc_method <- 1L

  art_alloc_method
}

get_dp_art_prop_alloc <- function(dp){

  if(exists_dptag(dp, "<NewARTPatAlloc MV>"))
    art_prop_alloc <- as.numeric(dpsub(dp, "<NewARTPatAlloc MV>", 2, 4:5))
  else
    art_prop_alloc <- c(0.5, 0.5)
  names(art_prop_alloc) <- c("mx", "elig")

  art_prop_alloc
}

get_dp_scale_cd4_mort <- function(dp){

  vers_str <- dpsub(dp, "<ValidVers MV>", 2, 4)
  version <- sub("^([0-9\\.]+).*", "\\1", vers_str)
  betav <- if(grepl("Beta", vers_str))
             as.numeric(sub(".*Beta ([0-9]+)$", "\\1", vers_str))
           else
             NA
  if(version >= "5.73" && (betav >= 15 | is.na(betav)))
    scale_cd4_mort <- 1L
  else
    scale_cd4_mort <- 0L

  scale_cd4_mort
}

get_dp_artmx_timerr <- function(dp, proj_years){

  col_idx <- 3 + seq_along(proj_years)
  
  artmx_timerr <- array(1.0, c(3, length(proj_years)), list(artdur=c("ART0MOS", "ART6MOS", "ART1YR"), year = proj_years))
  if(exists_dptag(dp, "<MortalityRates MV>")){
    val <- as.numeric(dpsub(dp, "<MortalityRates MV>", 2, col_idx))
    artmx_timerr["ART0MOS", ] <- val
    artmx_timerr["ART6MOS", ] <- val
    artmx_timerr["ART1YR", ] <- val
  } else if(exists_dptag(dp, "<MortalityRates MV2>")) {
    val <- array(as.numeric(unlist(dpsub(dp, "<MortalityRates MV2>", 2:3, col_idx))), c(2, length(proj_years)))
    artmx_timerr["ART0MOS", ] <- val[1, ]
    artmx_timerr["ART6MOS", ] <- val[1, ]
    artmx_timerr["ART1YR", ] <- val[2, ]
  }

  artmx_timerr
}


  



#' Calculate ASFR from TFR and fertility distribution
#' 
#' @param tfr vector of annual TFR values
#' @param asfd array of proportion of births by 5 year age group 15-49
#'
#' @return array of age-specific fertility rate by single-year of age 15-49.
calc_asfr <- function(tfr, asfd){

  ## Normalise the ASFD to sum exactly to 1.0
  asfd_sum <- colSums(asfd)
  asfd_sum[asfd_sum == 0.0] <- 1.0
  asfd <- sweep(asfd, 2, asfd_sum, "/")
  
  asfr <- apply(asfd / 5, 2, rep, each=5)
  asfr <- sweep(asfr, 2, tfr, "*")
  dimnames(asfr) <- list(age = 15:49, year = dimnames(asfd)[[2]])

  asfr
}

calc_netmigr <- function(totnetmig, netmigagedist, Sx){

  ## Normalise netmigagedist to sum to exactly 1.0
  netmigagedist_sum <- colSums(netmigagedist)
  netmigagedist_sum[netmigagedist_sum == 0.0] <- 1.0
  netmigagedist <- sweep(netmigagedist, 2:3, netmigagedist_sum, "/")
  
  netmigr5 <- sweep(netmigagedist, 2:3, totnetmig, "*")
  A <- create_beers(17)
  netmigr <- apply(netmigr5, 2:3, function(x) A %*% x)
  dimnames(netmigr) <- specpop_dimnames(colnames(totnetmig))

  ## For age <5 years, Spectrum disaggregates netmigration proportional to survival in
  ## the **base year**.
  u5prop <- array(dim = c(5, 2))
  u5prop[1, ] <- Sx[1, , 1] * 2
  u5prop[2, ] <- Sx[2, , 1] * u5prop[1, ]
  u5prop[3, ] <- Sx[3, , 1] * u5prop[2, ]
  u5prop[4, ] <- Sx[4, , 1] * u5prop[3, ]
  u5prop[5, ] <- Sx[5, , 1] * u5prop[4, ]
  
  u5prop <- sweep(u5prop, 2, colSums(u5prop), "/")
  
  netmigr[1:5 , 1, ] <- u5prop[ , 1, drop = FALSE] %*% netmigr5[1, 1, ]
  netmigr[1:5 , 2, ] <- u5prop[ , 2, drop = FALSE] %*% netmigr5[1, 2, ]
  
  netmigr
}

#' Get country name from parsed PJN
#' @param pjn parsed PJN file
#' @details
#' `pjn` should be via `first90_read_csv_character(pjn_file)`
#' @export
get_pjn_country <- function(pjn){
  cc <- as.integer(pjn[which(pjn[, 1] == "<Projection Parameters>") + 2, 4])
  idx <- which(spectrum5_countrylist$Code == cc)
  return(spectrum5_countrylist$Country[idx])
}

#' Get subnational region from parsed PJN
#' @param pjn parsed PJN file
#' @details
#' `pjn` should be via `first90_read_csv_character(pjn_file)`
#' @export
get_pjn_region <- function(pjn){
  region <- pjn[which(pjn[, 1] == "<Projection Parameters - Subnational Region Name2>") + 2, 4]
  if (is.na(region) || region == "") 
    return(NULL)
  else
    return(region)
}

get_pjn_iso3 <- function(pjn){
  cc <- as.integer(pjn[which(pjn[, 1] == "<Projection Parameters>") + 2, 4])
  idx <- which(spectrum5_countrylist$Code == cc)
  return(spectrum5_countrylist$Code[idx])
}

#' @export
file_in_zip <- function(zfile, ext){
  ff <- grep(ext, unzip(zfile, list = TRUE)$Name, value = TRUE)
  unz(zfile, ff)
}

#' @export 
read_country <- function(pjnz = NULL, pjn_file = NULL){

  if(!is.null(pjnz) && is.null(pjn_file))
    pjn_file <- file_in_zip(pjnz, ".PJN$")
  else if(!is.null(pjnz) && !is.null(pjn_file))
    stop("Must provide either a Spectrum .PJNZ file or a .PJN file. Do not provide both.")

  pjn <- first90_read_csv_character(pjn_file)
  get_pjn_country(pjn)
}

#' @export 
read_region <- function(pjnz = NULL, pjn_file = NULL){

  if(!is.null(pjnz) && is.null(pjn_file))
    pjn_file <- file_in_zip(pjnz, ".PJN$")
  else if(!is.null(pjnz) && !is.null(pjn_file))
    stop("Must provide either a Spectrum .PJNZ file or a .PJN file. Do not provide both.")

  pjn <- first90_read_csv_character(pjn_file)
  get_pjn_region(pjn)
}
