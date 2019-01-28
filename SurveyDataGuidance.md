## Guidance on survey data

Please provide survey data on the proportion of people ever tested by sex and age group.
This needs to be a data-frame with the following columns. If these columns names are wrong the model
will either not run, or the results will not be accurate.

* country: Must exactly match the value extracted from your Spectrum file:

```
 pjnz <- "~/path/to/yourfile.PJNZ"
 country <- read_country(pjnz)
```

* surveyId
* year: Year in which survey was conducted; year of survey fieldwork midpoint if survey spanned multiple years.
* agegr: 15-24, 25-34, 35-49, 50+, 15-49 or 15+
* sex: male, female or both.
* hivstatus: positive, negative or all.
* est: Estimate for proportion ever tested for HIV; as a proportion (e.g. 0.876 rather than 87.6%).
* se: As a proportion.
* ci_l: As a proportion.
* ci_u: As a proportion.
* counts: Unweighted counts of survey respondents in stratification group.
* outcome: Only data with outcome "evertest" will be used in the model fitting.
