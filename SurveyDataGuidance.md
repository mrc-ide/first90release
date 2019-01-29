## Guidance on survey data

Please provide survey data on the proportion of people ever tested by sex, age group (15-24, 25-34, and 35-49 years), 
and HIV serostatus (if available). This needs to be a data-frame with the following columns. 
If these columns names are wrong the model will either not run, or the results will not be accurate. 
In order to visualize model fits, it is useful to also provide information on the proportion of people ever tested 
among the 15-49 years old overall, and stratified by sex and HIV serostatus (if available).

* country: Must exactly match the value extracted from your Spectrum file:

```
 pjnz <- "~/path/to/yourfile.PJNZ"
 country <- read_country(pjnz)
```

* surveyId: Provide a unique identifier for your survey (each survey must have a unique name).
* year: Year in which survey was conducted; year of survey fieldwork midpoint if survey spanned multiple years.
* agegr: "15-24", "25-34", "35-49", (and "15-49 for plotting model fits). 
* sex: "male", "female" or "both".
* hivstatus: "positive", "negative" or "all".
* est: Estimate for proportion ever tested for HIV who received the result of their last HIV test; as a proportion (e.g. 0.876 rather than 87.6%).
* se: Standard Error of the estimate (as a proportion). Should take into account survey design effects.
* ci_l: Lower limit of the 95% confidence interval of the survey estimate (as a proportion).
* ci_u: Uower limit of the 95% confidence interval of the survey estimate (as a proportion).
* counts: Unweighted counts of the number of survey respondents included in the stratification group.
* outcome: Only data with outcome "evertest" will be used in the model fitting (i.e., include "evertest" for all rows).
