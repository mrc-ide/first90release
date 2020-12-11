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
