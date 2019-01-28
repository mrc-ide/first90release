---
output: github_document
---

<!-- README.md is generated from README.Rmd.
Please edit that file and run `knitr::knit("README.Rmd") from the root of this directory` -->


# first90

UNAIDS put forward the ambitious 90-90-90 target to end the AIDS epidemic by 2030. This target aims for 90% of people living with HIV (PLHIV) to be aware of their HIV-positive status, 90% of those diagnosed to receive antiretroviral therapy, and 90% of those on treatment to have a suppressed viral load by 2020 (each reaching 95% by 2030). HIV testing remains an important bottleneck in this cascade, however, and obtaining reliable epidemiological data on the proportion of PLHIV aware of their status is difficult. Such information is nevertheless crucial to effectively monitor HIV prevention efforts. Tracking progress towards achievement of this “first 90” target could be improved by combining population-based surveys and programmatic data on the number of HIV tests performed (and yield) in a coherent deterministic/statistical model. This type of integrative systems modelling is especially useful to fully consider HIV incidence, mortality, testing behaviours, as well as to coherently combine different sources of data. 

The goal of the first90 package is to provide annual estimates of the proportion of PLHIV that are aware of their status, by combining estimates of PLHIV from EPP/Spectrum, annual programmatic data on the number of HIV tests performed (and yield), and nationally-representative survey of HIV testing behaviors. 

## Installation

Install via Github using `devtools`:

``` r
devtools::install_github("mrc-ide/first90release")
```

## Example: Malawi

This example demonstrates basic model steps.


```r
# Read PJNZ file(s)
pjnz <- "~/Downloads/Malawi_2018_version_8.PJNZ"
cnt <- first90::read_country(pjnz)

fp <- first90::prepare_inputs(pjnz)
# first90::prepare_inputs can also take a list of files, if using regional files
# e.g. fp <- first90::prepare_inputs(list.files("~/Documents/Data/", "CotedIvoire.*PJNZ$", full.names=TRUE, ignore.case=TRUE))

# We visualize the PJNZ data
first90::plot_pjnz(fp)
```

![plot of chunk example](man/figures/README-example-1.png)

The following functions enable users to produce invidual plots.
```
pjnz_summary <- first90::get_pjnz_summary_data(fp)
first90::plot_pjnz_pop(pjnz_summary)
first90::plot_pjnz_plhiv(pjnz_summary)
first90::plot_pjnz_prv(pjnz_summary)
first90::plot_pjnz_inc(pjnz_summary)
```


```r
age_group <- c('15-24','25-34','35-49')
# Import and prepare your survey data. See [guidance](SurveyDataGuidance.md)
survey_hts <- data.frame(country="Malawi",
                                  surveyid="Survey1",
                                  year=2000,
                                  agegr="15-99",
                                  sex="both",
                                  outcome="evertest",
                                  hivstatus="positive",
                                  est=0.553,
                                  se=0.0159,
                                  ci_l=2.5652e-12,
                                  ci_u=8958e-12,
                                  counts=16168)

dat <- first90::select_hts(survey_hts, cnt, age_group)

# Import and prepare your programmatic data. See [guidance](ProgramDataGuidance.md)
prgm_dat <- data.frame(country = "Malawi",
                            year = 2010,
                            sex = 'both',
                            tot = 215269,
                            totpos = 50115,
                            vct = NA,
                            vctpos = NA,
                            anc = NA,
                            ancpos = NA)
prg_dat <- first90::select_prgmdata(prgm_dat, cnt, age_group)

# We visualize the program data
first90::plot_inputdata(prg_dat, fp)
```

![plot of chunk unnamed-chunk-1](man/figures/README-unnamed-chunk-1-1.png)

The following functions enable users to produce invidual plots.
```
first90::plot_input_tot(prgm_dat, fp)
first90::plot_input_totpos(prgm_dat, fp)
first90::plot_input_anctot(prgm_dat, fp)
first90::plot_input_ancpos(prgm_dat, fp)
```

```r
# ---- Enter parameters here ----
# We create the likelihood data
likdat <- first90::prepare_hts_likdat(dat, prg_dat, fp)

# Starting parameters
data("theta0", package="first90")
first90::ll_hts(theta0, fp, likdat)
#> [1] -3428.327

opt <- optim(theta0, ll_hts, fp = fp, likdat = likdat, method = "BFGS", 
             control = list(fnscale = -1, trace = 4, REPORT = 1, maxit = 250), hessian = TRUE)
#> initial  value 3428.327287 
#> iter   2 value 370.374792
#> iter   3 value 212.081672
#> iter   4 value 186.126159
#> iter   5 value 141.537213
#> iter   6 value 101.059665
#> iter   7 value 92.621121
#> iter   8 value 87.475615
#> iter   9 value 84.788073
#> iter  10 value 75.022805
#> iter  11 value 61.528578
#> iter  12 value 55.627036
#> iter  13 value 54.013246
#> iter  14 value 52.348387
#> iter  15 value 49.595268
#> iter  16 value 47.570854
#> iter  17 value 44.622040
#> iter  18 value 41.998845
#> iter  19 value 40.912793
#> iter  20 value 39.263882
#> iter  21 value 38.740845
#> iter  22 value 38.238546
#> iter  23 value 37.284149
#> iter  24 value 36.847265
#> iter  25 value 36.578046
#> iter  26 value 36.162281
#> iter  27 value 35.842953
#> iter  28 value 35.687427
#> iter  29 value 35.506092
#> iter  30 value 35.282857
#> iter  31 value 35.224500
#> iter  32 value 34.995619
#> iter  33 value 34.940506
#> iter  34 value 34.890671
#> iter  35 value 34.857324
#> iter  36 value 34.801386
#> iter  37 value 34.748816
#> iter  38 value 34.655957
#> iter  39 value 34.565003
#> iter  40 value 34.434560
#> iter  41 value 34.106836
#> iter  42 value 33.511850
#> iter  43 value 33.273834
#> iter  44 value 33.010539
#> iter  45 value 32.636078
#> iter  46 value 32.408091
#> iter  47 value 32.117210
#> iter  48 value 31.766076
#> iter  49 value 31.202110
#> iter  50 value 30.658590
#> iter  51 value 29.788505
#> iter  52 value 29.136107
#> iter  53 value 28.328386
#> iter  54 value 28.210781
#> iter  55 value 28.009187
#> iter  56 value 27.279623
#> iter  57 value 26.819745
#> iter  58 value 26.147700
#> iter  59 value 25.922341
#> iter  60 value 25.811214
#> iter  61 value 25.795109
#> iter  62 value 25.782505
#> iter  63 value 25.743289
#> iter  64 value 25.736549
#> iter  65 value 25.732321
#> iter  66 value 25.728007
#> iter  67 value 25.726405
#> iter  68 value 25.709005
#> iter  69 value 25.707754
#> iter  70 value 25.699064
#> iter  71 value 25.678734
#> iter  72 value 25.651630
#> iter  73 value 25.629433
#> iter  74 value 25.609648
#> iter  75 value 25.596536
#> iter  76 value 25.570177
#> iter  77 value 25.563056
#> iter  78 value 25.548120
#> iter  79 value 25.546800
#> iter  80 value 25.546575
#> iter  81 value 25.543256
#> iter  82 value 25.541494
#> iter  83 value 25.540472
#> iter  84 value 25.540166
#> iter  85 value 25.539892
#> iter  86 value 25.539629
#> iter  87 value 25.539454
#> iter  88 value 25.539222
#> iter  89 value 25.539208
#> iter  90 value 25.538520
#> iter  91 value 25.538021
#> iter  92 value 25.537742
#> iter  93 value 25.537733
#> iter  94 value 25.537686
#> iter  95 value 25.537428
#> iter  96 value 25.536062
#> iter  97 value 25.535666
#> iter  98 value 25.535663
#> iter  99 value 25.535657
#> iter 100 value 25.535655
#> iter 100 value 25.535655
#> iter 101 value 25.535481
#> iter 102 value 25.535227
#> iter 103 value 25.535013
#> iter 104 value 25.534804
#> iter 105 value 25.534688
#> iter 105 value 25.534688
#> final  value 25.534688 
#> converged

simul <- first90::simul.test(opt, fp, sim = 400)

# ---- Plots for FITS ----
fp <- first90::create_hts_param(opt$par, fp)
mod <- first90::simmod(fp)

# ---- The Fitted Parameters ----
first90::optimized_par(opt)
#> Loading required package: Matrix
#>                               Parameter_Name Estimate  LCI  UCI
#> 1                    RR testing: men in 2005     0.85 0.59 1.00
#> 2                    RR testing: men in 2012     0.94 0.81 1.02
#> 3                         RR re-testing 2010     1.94 1.04 5.61
#> 4                         RR re-testing 2015     1.94 1.04 5.68
#> 5                  RR testing: PLHIV unaware     1.56 1.50 1.62
#> 6  RR re-testing: PLHIV aware (not ART) 2010     1.02 0.11 4.86
#> 7  RR re-testing: PLHIV aware (not ART) 2017     1.02 0.08 5.46
#> 8  RR re-testing: PLHIV on ART (*RR not ART)     0.14 0.01 0.74
#> 9                         RR among 25-34 men     1.55 1.37 1.75
#> 10                          RR among 35+ men     3.17 2.91 3.42
#> 11                      RR among 25-34 women     2.19 1.99 2.39
#> 12                        RR among 35+ women     3.89 3.46 4.28
#> 13                        RR OI Dx (ART Cov)     1.60 1.50 1.60

# ---- Functions for individuals model fits ----
out_evertest <- first90::get_out_evertest(mod, fp)

first90::plot_out(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
```

![plot of chunk unnamed-chunk-2](man/figures/README-unnamed-chunk-2-1.png)
  
The model fits by age and sex.

```r
first90::plot_out_strat(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
```

![plot of chunk unnamed-chunk-3](man/figures/README-unnamed-chunk-3-1.png)
The following functions enable users to produce invidual plots.
```
first90::plot_out_nbtest(mod, fp, likdat, cnt, simul)
first90::plot_out_nbpostest(mod, fp, likdat, cnt, simul)
first90::plot_out_evertestneg(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
first90::plot_out_evertestpos(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
first90::plot_out_evertest(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
first90::plot_out_90s(mod, fp, likdat, cnt, out_evertest, survey_hts, simul)
first90::plot_out_evertest_fbyage(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
first90::plot_out_evertest_mbyage(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
```
We can compare HIV tests' positivity through time, the estimated *true* yield of new HIV diagnoses, and compare those to population-level HIV prevalence.


```r
par(mfrow = c(1,1))
first90::plot_prv_pos_yld(mod, fp, likdat, cnt, yr_pred = 2018)
```

![plot of chunk unnamed-chunk-4](man/figures/README-unnamed-chunk-4-1.png)

We can also examine some ouptuts related to the distribution of HIV tests performed in those susceptibles to HIV and PLHIV by different awareness and treatment status. (First sets of graphs is on the absolute scale, second one on the relative scale.)


```r
par(mfrow = c(1,2))
first90::plot_retest_test_neg(mod, fp, likdat, cnt)
first90::plot_retest_test_pos(mod, fp, likdat, cnt)
```

![plot of chunk unnamed-chunk-5](man/figures/README-unnamed-chunk-5-1.png)

```r

par(mfrow = c(1,2))
first90::plot_retest_test_neg(mod, fp, likdat, cnt, relative = TRUE)
first90::plot_retest_test_pos(mod, fp, likdat, cnt, relative = TRUE)
```

![plot of chunk unnamed-chunk-5](man/figures/README-unnamed-chunk-5-2.png)

Finally, tabular outputs can be obtained by using the following functions.

```r
# ---- Tabular outputs ----
first90::tab_out_evertest(mod, fp, simul = simul)
#>   year  outcome agegr  sex hivstatus value lower upper
#> 1 2010 evertest 15-49 both       all   9.8   8.4  12.4
#> 2 2011 evertest 15-49 both       all  11.6   9.9  14.0
#> 3 2012 evertest 15-49 both       all  13.2  11.3  16.2
#> 4 2013 evertest 15-49 both       all  14.6  12.4  18.9
#> 5 2014 evertest 15-49 both       all  15.9  13.2  21.4
#> 6 2015 evertest 15-49 both       all  17.1  13.7  23.9
#> 7 2016 evertest 15-49 both       all  18.2  14.2  26.6
#> 8 2017 evertest 15-49 both       all  19.1  14.7  28.9
#> 9 2018 evertest 15-49 both       all  20.0  15.0  30.7
first90::tab_out_aware(mod, fp, simul = simul)
#>   year outcome agegr  sex hivstatus value lower upper
#> 1 2010   aware 15-49 both  positive  27.2  27.2  29.2
#> 2 2011   aware 15-49 both  positive  33.8  33.8  34.6
#> 3 2012   aware 15-49 both  positive  40.5  40.5  41.2
#> 4 2013   aware 15-49 both  positive  46.5  46.5  47.7
#> 5 2014   aware 15-49 both  positive  51.7  51.6  53.3
#> 6 2015   aware 15-49 both  positive  56.8  56.8  58.6
#> 7 2016   aware 15-49 both  positive  62.5  62.5  64.3
#> 8 2017   aware 15-49 both  positive  67.8  67.7  69.6
#> 9 2018   aware 15-49 both  positive  72.3  72.2  74.1
first90::tab_out_nbaware(mod, fp)
#>   year      outcome agegr  sex hivstatus  value
#> 1 2010 number aware 15-49 both  positive 195845
#> 2 2011 number aware 15-49 both  positive 249283
#> 3 2012 number aware 15-49 both  positive 306953
#> 4 2013 number aware 15-49 both  positive 360294
#> 5 2014 number aware 15-49 both  positive 407764
#> 6 2015 number aware 15-49 both  positive 454584
#> 7 2016 number aware 15-49 both  positive 504249
#> 8 2017 number aware 15-49 both  positive 549282
#> 9 2018 number aware 15-49 both  positive 587187
first90::tab_out_artcov(mod, fp)
#>   year outcome agegr  sex hivstatus value
#> 1 2010  artcov   15+ both  positive  29.4
#> 2 2011  artcov   15+ both  positive  36.7
#> 3 2012  artcov   15+ both  positive  44.6
#> 4 2013  artcov   15+ both  positive  50.3
#> 5 2014  artcov   15+ both  positive  55.8
#> 6 2015  artcov   15+ both  positive  60.0
#> 7 2016  artcov   15+ both  positive  66.8
#> 8 2017  artcov   15+ both  positive  71.4
#> 9 2018  artcov   15+ both  positive  75.9
```
