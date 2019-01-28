## Guidance on program data

This should be data sourced from national testing programs. For each year provide either sex 
aggregated data (1 row with sex = "both") or sex disaggregated data (2 rows with sex="male" and 
sex="female" respectively). This should be for age-group 15-99 only.

This needs to be a data-frame with the following columns. If these columns names are wrong the model
will either not run, or the results will not be accurate.

* country: Must match exactly the value extracted from your Spectrum file:
```
 pjnz <- "~/path/to/yourfile.PJNZ"
 country <- read_country(pjnz)
```

* year: Year in which the survey was conducted; year of survey fieldwork midpoint if survey spanned multiple years.
* sex: For each year EITHER provide disaggregated values for "male" and "female", or an aggregated value for "both"
* agegr: This should be "15-99" for all rows.
* tot: This is the annual number of tests performed at the national level among the population aged 15+ years of age. 
This number should be equal to the total number of tests administered as part of HIV Testing and Counseling (HTC) and 
during antenatal care (ANC), and for which the clients received the results.
* totpos: Out of the total annual number of tests, how many were found to be HIV positive. This number should be equal 
to the number of positive tests found during HTC (in non-pregnant population) and during ANC among pregnant women.
* vct: Total annual number of tests performed in the population aged 15+ years outside of ANC services, and for which
 clients received the results. If only the overall total is available, please input NA.
* vctpos: Annual number of tests that were found to be positive for HIV outside of ANC services.
* anc: Annual number of pregnant women tested for HIV (and that received their results) as part of ANC services.
* ancpos: Annual number of pregnant women found to be HIV positive during ANC services.
