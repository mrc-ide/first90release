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
#' @examples
#'
#' write.csv(spectrum_output_table(mod, fp), 
#'           "~/Downloads/Malawi-shiny90-example-output.csv", 
#'            row.names = FALSE)
#' @export 
spectrum_output_table <- function(mod, fp) {

  year <- fp$ss$proj_start + seq_len(fp$ss$PROJ_YEARS) - 1L

  plhiv <- t(colSums(mod[,,2,]))
  colnames(plhiv) <- paste0("plhiv_", c("m", "f"))
  
  onart <- t(colSums(attr(mod, "artpop"),,3))
  colnames(onart) <- paste0("onart_", c("m", "f"))

  aware <- onart + t(colSums(attr(mod, "diagnpop"),,2))
  colnames(aware) <- paste0("aware_", c("m", "f"))

  evertest <- aware + t(colSums(attr(mod, "testnegpop")[,,2,]))
  colnames(evertest) <- paste0("evertest_", c("m", "f"))

  data.frame(year, plhiv, evertest, aware, onart)
}
  
