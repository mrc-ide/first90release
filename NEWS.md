# first90 1.5.2

* Update initial values to median of 2021 .shiny90 file final parameter estiamtes.

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
