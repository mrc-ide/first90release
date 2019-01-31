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
#' @examples
#'
#' \dontrun{
#' write.csv(spectrum_output_table(mod, fp), 
#'           "~/Downloads/Malawi-shiny90-example-output.csv", 
#'            row.names = FALSE)
#' }
#' @export 
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

  plhiv <- t(colSums(mod[,,2,out_idx]))
  colnames(plhiv) <- paste0("plhiv_", c("m", "f"))

  onart <- t(fp$art15plus_num[,out_idx])
  colnames(onart) <- paste0("onart_", c("m", "f"))

  evertest <- cbind(evertest_m = evertest_m$value,
                    evertest_f = evertest_f$value) * plhiv
  aware <- cbind(aware_m = aware_m$value,
                 aware_f = aware_f$value) * plhiv

  data.frame(year = year_out,
             plhiv,
             evertest,
             aware,
             onart)
}
