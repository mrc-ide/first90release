
<!-- README.md is generated from README.Rmd. Please edit that file -->
first90
=======

UNAIDS put forward the ambitious 90-90-90 target to end the AIDS epidemic by 2030. This target aims for 90% of people living with HIV (PLHIV) to be aware of their HIV-positive status, 90% of those diagnosed to receive antiretroviral therapy, and 90% of those on treatment to have a suppressed viral load by 2020 (each reaching 95% by 2030). HIV testing remains an important bottleneck in this cascade, however, and obtaining reliable epidemiological data on the proportion of PLHIV aware of their status is difficult. Such information is nevertheless crucial to effectively monitor HIV prevention efforts. Tracking progress towards achievement of this “first 90” target could be improved by combining population-based surveys and programmatic data on the number of HIV tests performed (and yield) in a coherent deterministic/statistical model. This type of integrative systems modelling is especially useful to fully consider HIV incidence, mortality, testing behaviours, as well as to coherently combine different sources of data.

The goal of the first90 package is to provide annual estimates of the proportion of PLHIV that are aware of their status, by combining estimates of PLHIV from EPP/Spectrum, annual programmatic data on the number of HIV tests performed (and yield), and nationally-representative survey of HIV testing behaviors.

Installation
------------

Install via Github using `devtools`:

``` r
devtools::install_github("mrc-ide/first90")
```

This is a private repository. Installation via Github requires a Personal Access Token (PAT). To create a token, visit <https://github.com/settings/tokens/new>. Ensure that the box "repo" is ticked, and click "Generate token".

Create a file in home directory `~/.Renviron` with the following line (replacing `<YOUR TOKEN>` with token that appears):

    GITHUB_PAT = <YOUR TOKEN>

Save the file, restart R, and use the command `install_github()` above.

Example: Malawi
---------------

This example demonstrates basic model steps.

``` r
cnt <- "Malawi"
age_group <- c('15-24','25-34','35-49')

# Read PJNZ file(s)
pjnz <- "~/Google Drive/McGill/Research/Wisbech/Spectrum files/Malawi_2018_version_8.PJNZ"
fp <- prepare_inputs(pjnz)

# We visualize the PJNZ data
plot_pjnz(fp)
```

![](man/figures/README-example-1.png)

The following functions enable users to produce invidual plots.

    pjnz_summary <- get_pjnz_summary_data(fp)
    plot_pjnz_pop(pjnz_summary)
    plot_pjnz_plhiv(pjnz_summary)
    plot_pjnz_prv(pjnz_summary)
    plot_pjnz_inc(pjnz_summary)

``` r
# All surveys should be included. If serology was not collected, the model should
# be fitted overall
data(survey_hts)
dat <- dat_test <- select_hts(survey_hts, cnt, age_group)

# We prepare the program data
data(prgm_dat)
prg_dat <- select_prgmdata(prgm_dat, cnt, age_group)

# We visualize the program data
plot_inputdata(prg_dat, fp)
```

![](man/figures/README-unnamed-chunk-1-1.png)

The following functions enable users to produce invidual plots.

    plot_input_tot(prgm_dat, fp)
    plot_input_totpos(prgm_dat, fp)
    plot_input_anctot(prgm_dat, fp)
    plot_input_ancpos(prgm_dat, fp)

``` r
# ---- Enter parameters here ----
# We create the likelihood data
likdat <- prepare_hts_likdat(dat, prg_dat, fp)

# Starting parameters
data(theta0)
ll_hts(theta0, fp, likdat)
#> [1] -13801.64

opt <- optim(theta0, ll_hts, fp = fp, likdat = likdat, method = "BFGS", 
             control = list(fnscale = -1, trace = 4, REPORT = 1, maxit = 250), hessian = TRUE)
#> initial  value 13801.640245 
#> iter   2 value 8534.729033
#> iter   3 value 6012.879598
#> iter   4 value 5959.329582
#> iter   5 value 5347.251127
#> iter   6 value 4614.303482
#> iter   7 value 4411.501226
#> iter   8 value 4347.134154
#> iter   9 value 4280.783307
#> iter  10 value 4183.540230
#> iter  11 value 4124.194076
#> iter  12 value 4035.166598
#> iter  13 value 4005.548672
#> iter  14 value 3975.052745
#> iter  15 value 3938.852808
#> iter  16 value 3914.482890
#> iter  17 value 3877.065883
#> iter  18 value 3852.556175
#> iter  19 value 3816.662841
#> iter  20 value 3602.969351
#> iter  21 value 3474.652594
#> iter  22 value 3364.761592
#> iter  23 value 3186.673116
#> iter  24 value 3068.396406
#> iter  25 value 2980.821757
#> iter  26 value 2883.062528
#> iter  27 value 2752.062613
#> iter  28 value 2695.152034
#> iter  29 value 2632.177829
#> iter  30 value 2605.798911
#> iter  31 value 2569.553112
#> iter  32 value 2511.754026
#> iter  33 value 2385.772192
#> iter  34 value 2296.874753
#> iter  35 value 2233.635205
#> iter  36 value 2191.844359
#> iter  37 value 2155.164806
#> iter  38 value 2001.722851
#> iter  39 value 1982.408574
#> iter  40 value 1934.754416
#> iter  41 value 1864.043628
#> iter  42 value 1841.536140
#> iter  43 value 1764.507459
#> iter  44 value 1647.811384
#> iter  45 value 1542.302552
#> iter  46 value 1382.316787
#> iter  47 value 1272.413850
#> iter  48 value 1253.242377
#> iter  49 value 1130.062404
#> iter  50 value 1019.733511
#> iter  51 value 996.674939
#> iter  52 value 985.160317
#> iter  53 value 980.457466
#> iter  54 value 971.728088
#> iter  55 value 968.203247
#> iter  56 value 965.470428
#> iter  57 value 963.480493
#> iter  58 value 961.978989
#> iter  59 value 960.360585
#> iter  60 value 959.249424
#> iter  61 value 958.327285
#> iter  62 value 957.745761
#> iter  63 value 957.517254
#> iter  64 value 957.441684
#> iter  65 value 957.395900
#> iter  66 value 957.355779
#> iter  67 value 957.303001
#> iter  68 value 957.235538
#> iter  69 value 957.132734
#> iter  70 value 956.976147
#> iter  71 value 956.772458
#> iter  72 value 956.573834
#> iter  73 value 956.444601
#> iter  74 value 956.383929
#> iter  75 value 956.351233
#> iter  76 value 956.332351
#> iter  77 value 956.323973
#> iter  78 value 956.319430
#> iter  79 value 956.317595
#> iter  80 value 956.316623
#> iter  81 value 956.316391
#> iter  82 value 956.316336
#> iter  82 value 956.316329
#> iter  82 value 956.316328
#> final  value 956.316328 
#> converged

simul <- simul.test(opt, fp, sim = 400)
#> Loading required package: MASS
#> Loading required package: Matrix

# ---- Plots for FITS ----
fp <- create_hts_param(opt$par, fp)
mod <- simmod(fp)

# ---- The Fitted Parameters ----
optimized_par(opt)
#>                               Parameter_Name Estimate  LCI  UCI
#> 1                    RR testing: men in 2005     0.81 0.74 0.88
#> 2                    RR testing: men in 2012     0.44 0.41 0.46
#> 3                              RR re-testing     1.13 1.09 1.19
#> 4                  RR testing: PLHIV unaware     0.92 0.82 1.04
#> 5  RR re-testing: PLHIV aware (not ART) 2010     1.70 1.50 1.93
#> 6  RR re-testing: PLHIV aware (not ART) 2017     1.16 0.93 1.43
#> 7  RR re-testing: PLHIV on ART (*RR not ART)     0.03 0.00 0.17
#> 8                         RR among 25-34 men     1.53 1.38 1.68
#> 9                           RR among 35+ men     0.82 0.74 0.92
#> 10                      RR among 25-34 women     1.28 1.20 1.37
#> 11                        RR among 35+ women     0.62 0.58 0.66
#> 12            % of OI tested for HIV in 2015     0.82 0.74 0.86
#> 13      SD penatly for testing rate (female)     0.10 0.10 0.10
#> 14   SD penatly for RR PLHIV aware (not ART)     0.25 0.01 0.25

# ---- Functions for individuals model fits ----
data(survey_hts)
out_evertest <- get_out_evertest(mod, fp)

plot_out(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
```

![](man/figures/README-unnamed-chunk-2-1.png)

The model fits by age and sex.

``` r
plot_out_strat(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
```

![](man/figures/README-unnamed-chunk-3-1.png) The following functions enable users to produce invidual plots.

    plot_out_nbtest(mod, fp, likdat, cnt, simul)
    plot_out_nbpostest(mod, fp, likdat, cnt, simul)
    plot_out_evertestneg(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
    plot_out_evertestpos(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
    plot_out_evertest(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
    plot_out_90s(mod, fp, likdat, cnt, out_evertest, simul)
    plot_out_evertest_fbyage(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)
    plot_out_evertest_mbyage(mod, fp, likdat, cnt, survey_hts, out_evertest, simul)

We can compare HIV tests' positivity through time, the estimated *true* yield of new HIV diagnoses, and compare those to population-level HIV prevalence.

``` r
par(mfrow = c(1,1))
plot_prv_pos_yld(mod, fp, likdat, cnt, yr_pred = 2018) 
```

![](man/figures/README-unnamed-chunk-4-1.png)

We can also examine some ouptuts related to the distribution of HIV tests performed in those susceptibles to HIV and PLHIV by different awareness and treatment status. (First sets of graphs is on the absolute scale, second one on the relative scale.)

``` r
par(mfrow = c(1,2))
plot_retest_test_neg(mod, fp, likdat, cnt)
plot_retest_test_pos(mod, fp, likdat, cnt)
```

![](man/figures/README-unnamed-chunk-5-1.png)

``` r

par(mfrow = c(1,2))
plot_retest_test_neg(mod, fp, likdat, cnt, relative = TRUE)
plot_retest_test_pos(mod, fp, likdat, cnt, relative = TRUE)
```

![](man/figures/README-unnamed-chunk-5-2.png)

Finally, tabular outputs can be obtained by using the following functions.

``` r
# ---- Tabular outputs ----
tab_out_evertest(mod, fp, simul = simul)
#>   year  outcome agegr  sex hivstatus value lower upper
#> 1 2010 evertest 15-49 both       all  58.6  58.3  59.3
#> 2 2011 evertest 15-49 both       all  63.6  63.3  64.2
#> 3 2012 evertest 15-49 both       all  67.3  67.1  67.9
#> 4 2013 evertest 15-49 both       all  70.0  69.8  70.5
#> 5 2014 evertest 15-49 both       all  72.2  72.1  72.7
#> 6 2015 evertest 15-49 both       all  74.5  74.3  75.0
#> 7 2016 evertest 15-49 both       all  77.6  77.4  78.1
#> 8 2017 evertest 15-49 both       all  81.2  81.0  81.8
#> 9 2018 evertest 15-49 both       all  83.6  83.3  84.7
tab_out_aware(mod, fp, simul = simul)
#>   year outcome agegr  sex hivstatus value lower upper
#> 1 2010   aware 15-49 both  positive  64.2  63.3  66.3
#> 2 2011   aware 15-49 both  positive  68.5  67.7  70.4
#> 3 2012   aware 15-49 both  positive  72.0  71.2  73.7
#> 4 2013   aware 15-49 both  positive  74.8  74.2  76.3
#> 5 2014   aware 15-49 both  positive  77.2  76.6  78.6
#> 6 2015   aware 15-49 both  positive  79.6  79.1  80.9
#> 7 2016   aware 15-49 both  positive  82.5  82.0  83.7
#> 8 2017   aware 15-49 both  positive  85.6  85.1  86.8
#> 9 2018   aware 15-49 both  positive  87.6  87.2  88.8
tab_out_artcov(mod, fp)
#>   year outcome agegr  sex hivstatus value
#> 1 2010  artcov 15-49 both  positive  24.0
#> 2 2011  artcov 15-49 both  positive  29.9
#> 3 2012  artcov 15-49 both  positive  36.9
#> 4 2013  artcov 15-49 both  positive  43.3
#> 5 2014  artcov 15-49 both  positive  48.9
#> 6 2015  artcov 15-49 both  positive  53.7
#> 7 2016  artcov 15-49 both  positive  59.4
#> 8 2017  artcov 15-49 both  positive  65.3
#> 9 2018  artcov 15-49 both  positive  70.1
```
