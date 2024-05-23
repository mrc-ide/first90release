# first90 1.6.10

* Implement recovery to next higher CD4 category following ART interruption for those on ART greater than one year.
* Qualify package names, fix some R CMD warnings and notes

# first90 1.6.9

* Bug fix: account for end-year net migration in the ART population in the first year of ART start (implemented in v1.6.0).

# first90 1.6.8

* Don't interpolate number aware outputs in Spectrum output table when using calendar year interpolation.

# first90 1.6.7

* Update PJNZ extraction for adult ART need Dec 31 for 2023 PJNZ files. Previously child ART was
  note recorded in the .DP file tag `<NeedARTDec31 MV>`, and so it was fine to extract the total
  value. Now child ART is recorded, and so need to sum the adult age groups only.
  
# first90 1.6.6

* Implement Spectrum Adult ART scalar adjustment. This is a user input that 
  allows the input number on ART to be adjusted by a scalar to account for 
  over/under-reporting of treatment numbers.

# first90 1.6.5

* Exclude NA values to set max `ylim` in `plot_prv_pos_yld()`. This addresses issue with this plot not showing in interface and preventing plot downloads when data NA values are simulated.

# first90 1.6.4

* Extend output plots and tables through 2022.

# first90 1.6.3

* Correct initial theta0 value for 2023 extended parameters.

# first90 1.6.2

* Use `vroom::vroom()` to read PJNZ files; much faster for reading `.DP` file.
* Handle Spectrum version numbers saved on Francophone locale devices: e.g. 6,13 instead of 6.13.
* In `read_specdp_demog_param()`, ensure no zero totals when normalising age-specific fertility 
  distribution and net-migration age distribution.
* Fix calculation for target number on treatment during transition from number to percent input.

# first90 1.6.1

Updates for 2023 UNAIDS estimates:

* Extend output arrays, plots, and tables to output 2022 results.
* Extend theta0 initial values for additional year random-walk parameters.

Other:

* Remove `data.table` dependency.

# first90 1.6.0

Updates for Spectrum transition from mid-year projection period to calendar year
projection period. Implemented in Spectrum 6.2 in November 2022.

Previous mid-year interpolation is still available in code by using argument
`create_fp(..., projection_period = "midyear")`. The default argument 
`projection_period = NULL` will choose the projection period based on the 
Spectrum inputs version number. For versions 6.19 and earlier,
`projection_period = "midyear"`. For versions 6.2 and later, `projection_period =
"calendar"`.


* Normalise `asfd` and `netmigagedist` to sum to exactly 1.0 before distributing
  to age groups.
* Disaggregate under-5 net migration proportional to paediatric survival probabilities
  in the base year.

* Net-migrations added at end of projection step, consistent with WPP 2022. No longer 
  (1) adjust net migration for half-period survival, nor (2) adjust net migration to 
  be half in current age group and half in next age group.
  
* End-year ART input interpolation adjusted to reflect calendar year projection period,
  such that ART input aligns to end of projection year. 
  _Note:_ since net migration is applied _after_ ART initiations, the number on ART will
  be scaled by the net migration proportions. Therefore the model output for number on
  ART will not exactly match the ART inputs.
  
* In likelihood calculations for household survey proportion ever tested, interopolate
  model outputs to mid-year. This affects `ll_evertest()`. Likelihood for annual 
  HIV tests and diagnoses (`ll_prgdat()`) is unaffected.
  
  
# first90 1.5.5

* Patch to `eppmod = "directinfections_hts"` option (v1.5.0) to avoid referencing 
  uninitialised memory.
* Fix error in ART deaths calculation in R version of simulation code. ART mortality time trend
  was omitted from deaths removed from total HIV population.

# first90 1.5.4

* Update time to diagnosis functions to account for CD4 at seroconversion in categories below 350.

# first90 1.5.3

* Add missing @export tags to ensure README.Rmd runs

# first90 1.5.2

* Update initial values to median of 2021 .shiny90 file final parameter estimates.

# first90 1.5.1

* Update output tables and plots show results through 2021 (missed in v1.5.0 update)

# first90 1.5.0

Updates for 2022 UNAIDS estimates:

* Extend output arrays, plots, and tables to output 2021 results.
* Extend theta0 initial values for additional year random-walk parameters.

* Intercalate new infections in each HIV time step. Equal number of infections are added in each time step. This is controlled by the option `eppmod = "directinfections_hts"`.
* Option `eppmod = "directinfections_hts"` is specified as the default option in 
  `prepare_inputs(pjnz)`. Therefore, new package version will not reproduce simulations
  from a previous model fit. Simulating from a previous fit requires manually specifying
  `fp$eppmod <- "directinfections"` after calling `prepare_inputs()`.
  
* `extract_pjnz()` reads whether custom population adjustment was used in the Spectrum file from the tag `"<RegionalAdjustPopCBState MV>"`. This is used to set `popadjust = TRUE` automatically if a custom population was used in Spectrum.

# first90 1.4.3

* Patch: in function `add_ss_indices()`, add argument `type.convert(..., as.is = TRUE)` 
  to suppress R 4.0 warning.
* Bugfix: remove duplicate declaration of `double incrate_g[NG];` in `EPP_DIRECTINCID` incidence
  option. This did not affect any results because this option is not used in Shiny90 application; 
  new infections by sex and age group were directly specified (option `EPP_DIRECTINFECTIONS`).

# first90 1.4.2

* Implement backward compatibility for simulating previous .shiny90 outputs.
  * The updates for each year add an additional knot, which changes the length
    of the parameter vector. This update parses the knot definition based on 
    the number of parameters.
  * Function `create_hts_param()` will throw and error if an unexpected parameter
    vector length is identified. 
  
# first90 1.4.1

Updates for 2021 UNAIDS estimates:

* Extended arrays to output 2020 results.

* Changed prior for baseline testing rate to normal(log(0.005), sd = 0.25). (likelihood.R#138)
  
* Revise HIV retesting rate to increase log-linearly betwee 2005 and 2010. The 
  implication of this is to reduce the retesting rate during the 2000s compared
  to the 2010s, which reduces implausibly high predictions for number of HIV
  tests conducted during the 2000s in the abscence of HTS programme data (Giguere et al 2020).

# first90 1.3.5

* Updated the calculation of 'retests' among HIV positive adults to include HIV positive
  tests conducted amongst adults who had previously tested HIV negative. Previously new
  HIV diagnoses were excluded from the count of retests and only HIV positive adults who 
  were previously diagnosed or on ART were counted as 'retests'.  This should not change
  any of the standard displayed outputs because by convention the 'retest' figures have
  focused on retests among HIV negative adults while separate results have been reported
  for retesting among HIV positive results.


# first90 1.3.4

* Add adjustment in Spectrum output CSV that the number aware of status is
  greater than the number on ART. This occurs in some cases when the 
  'internal' number on ART in Spectrum is greater than the input end of
  year numbers due to capping the number on ART below number of PLHIV 
  within age/sex strata.  This also means that if the number on ART 
  in Spectrum is greater than the number of PLHIV, then the number aware
  will also be greater than the number of PLHIV. (This occurs seldomly,
  but can occur.)

# first90 1.3.3

* Output the proportion of the undiagnosed population who were infected 
  in the past year in the Spectrum outputs download.

# first90 0.0.0.9000

* Added a `NEWS.md` file to track changes to the package.
