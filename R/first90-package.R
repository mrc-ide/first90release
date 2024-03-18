#' @useDynLib first90, .registration = TRUE
NULL

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
## usethis namespace: end
NULL

## Suppress R CMD check NOTEs about "no visible binding for
## global variable". These NOTEs aren't that important as we're using
## but the thousands of lines of NOTE they print does obscure real issues.
utils::globalVariables(c(
  "DT", "PROJ_YEARS", "ag.idx", "agegr", "aglast.idx", "country", "f.idx",
  "h.ag.span", "h.age15plus.idx", "h.fert.idx", "hAG", "hDS", "hTS",
  "hiv_steps_per_year", "hivn.idx", "hivp.idx", "hivstatus", "logit", "m.idx",
  "nearPD", "outcome", "p.age15plus.idx", "p.age15to49.idx", "p.fert.idx",
  "pAG", "pDS", "quantile", "sex", "tot", "totpos", "year"))
