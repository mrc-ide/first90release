context("test-plot-input-anc.R")

pjnz <- "sample_files/Malawi_2018_version_8.PJNZ"
filePath <- normalizePath(pjnz)

fp <- prepare_inputs(filePath)

test_that("plot_input_anctot does not error given NAs", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = 2010,
                                sex = 'both',
                                tot = NA,
                                totpos = NA,
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    expect_silent(plot_input_anctot(test_prgm_dat, fp))
})

test_that("plot_input_anctot does not error given only sex disaggregated data", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2010, 2011),
                                sex = c('female', 'male', 'female'),
                                tot = c(107634, 107634, 107634),
                                totpos = 50115,
                                vct = 107634,
                                vctpos = 25057,
                                anc = 107635,
                                ancpos = 25058)

    expect_silent(plot_input_anctot(test_prgm_dat, fp))
})

test_that("plot_input_anctot does not error given sex aggregated data", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = 2010,
                                sex = 'both',
                                tot = 215269,
                                totpos = 50115,
                                vct = 107634,
                                vctpos = 25057,
                                anc = 107635,
                                ancpos = 25058)

    expect_silent(plot_input_anctot(test_prgm_dat, fp))
})

test_that("plot_input_anctot does not error given mixed sex aggregated and sex dis-aggregated data", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2011, 2011),
                                sex = c('both', 'female', 'male'),
                                tot = 215269,
                                totpos = 50115,
                                vct = 107634,
                                vctpos = 25057,
                                anc = 107635,
                                ancpos = 25058)

    expect_silent(plot_input_anctot(test_prgm_dat, fp))
})
