# first90 1.3.5 

* Remove function `prb_dx_one_yr_cpp()` (commented out). This is creating error with windows compilation.

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
