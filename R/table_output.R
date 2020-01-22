#' Table to write CSV outputs for Spectrum
#'
#' @return a data.frame to write to CSV file for ingestion into Spectrum
#' 
#' @details
#' Presently this returns point estimates for age 15+ population by sex:
#'
#' * Number PLHIV,
#' * Ever tested among PLHIV
#' * Aware of HIV+ status
#' * On ART
#'
#' PLHIV is mid-year estimate. All other outcomes are end of year estimate.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' write.csv(spectrum_output_table(mod, fp), 
#'           "~/Downloads/Malawi-shiny90-example-output.csv", 
#'            row.names = FALSE)
#' }
#'
spectrum_output_table <- function(mod, fp) {

  aware_m <- get_out_aware(mod, fp, "15+", "male")
  aware_f <- get_out_aware(mod, fp, "15+", "female")
  evertest_m <- get_out_evertest(mod, fp, "15+", "male", "positive")
  evertest_f <- get_out_evertest(mod, fp, "15+", "female", "positive")
  
  aware_m$value <- end_of_year(aware_m$year, aware_m$value)
  aware_f$value <- end_of_year(aware_f$year, aware_f$value)
  evertest_m$value <- end_of_year(evertest_m$year, evertest_m$value)
  evertest_f$value <- end_of_year(evertest_f$year, evertest_f$value)
  
  year_out <- aware_m$year
  year_proj <- fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L
  out_idx <- year_proj %in% year_out

  plhiv <- t(colSums(mod[,, 2, out_idx]))
  colnames(plhiv) <- paste0("plhiv_", c("m", "f"))

  onart <- t(fp$art15plus_num[, out_idx])
  colnames(onart) <- paste0("onart_", c("m", "f"))
  
  ##  If we use the numbers on ART in mod (that get capped) all is fine.
  ##  We do not have the situation where kos < art
  ##  The code below verifies this...
      # onart_f <- colSums(attr(mod, "artpop")[, , 1:9, 2, out_idx, drop = FALSE], , 4)
      # onart_m <- colSums(attr(mod, "artpop")[, , 1:9, 1, out_idx, drop = FALSE], , 4)
      # onart <- cbind(onart_m, onart_f)
      # colnames(onart) <- paste0("onart_", c("m", "f"))
  
  evertest <- cbind(evertest_m = evertest_m$value,
                    evertest_f = evertest_f$value) * plhiv
  aware <- cbind(aware_m = aware_m$value,
                 aware_f = aware_f$value) * plhiv
  
  # Number adults 15+ undiagnosed and infected in the past year
  prb_dx_1yr_m <- pool_prb_dx_one_yr(mod, fp, year = c(2000:2019),
                   age = c("15-24","25-34", "35-49", "50-99"),
                   sex = c("male"))
  prb_dx_1yr_f <- pool_prb_dx_one_yr(mod, fp, year = c(2000:2019),
                   age = c("15-24","25-34", "35-49", "50-99"),
                   sex = c("female"))
  prb_dx_1yr <- cbind(prb_dx_1yr_m = c(prb_dx_1yr_m$prb1yr, 
                                       rep(NA, length(year_out) - nrow(prb_dx_1yr_m))),
                      prb_dx_1yr_f = c(prb_dx_1yr_f$prb1yr, 
                                       rep(NA, length(year_out) - nrow(prb_dx_1yr_m)))) 
  new_inf_m <- apply(fp$infections[, 1, out_idx], MARGIN = 2, FUN = sum)
  new_inf_f <- apply(fp$infections[, 2, out_idx], MARGIN = 2, FUN = sum)

  notdx_hiv_one_yr_m <- new_inf_m * (1 - prb_dx_1yr[, "prb_dx_1yr_m"])
  notdx_hiv_one_yr_f <- new_inf_f * (1 - prb_dx_1yr[, "prb_dx_1yr_f"])
  
  notdx_hiv_one_yr <- cbind(notdx_hiv_one_yr_m, notdx_hiv_one_yr_f) 
  
  val <- data.frame(year = year_out,
             plhiv,
             evertest,
             aware,
             onart,
             notdx_hiv_one_yr)
  
  ## If the numbers aware are lower than those on ART, we make them equal before
  ## importing back in Spectrum. This is necessary if countries overestimate their
  ## ART numbers... if that is the case, Spectrum will cap initiations if there are 
  ## not enough people to be initiated (based on sex, age, cd4 strata). Here, we are
  ## just making the numbers equal but countries should be strongly encouraged to 
  ## revisit their ART numbers.
  val$aware_m <- ifelse(val$aware_m < val$onart_m, val$onart_m, val$aware_m)
  val$aware_f <- ifelse(val$aware_f < val$onart_f, val$onart_f, val$aware_f)
  val$evertest_m <- ifelse(val$evertest_m < val$onart_m, val$onart_m, val$evertest_m)
  val$evertest_f <- ifelse(val$evertest_f < val$onart_f, val$onart_f, val$evertest_f)  
  
  return(val)
}
