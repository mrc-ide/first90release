context("test-input-plots.R")

pjnz <- "sample_files/Malawi_2018_version_8.PJNZ"
filePath <- normalizePath(pjnz)

fp <- prepare_inputs(filePath)

test_that("combine_rows can combine sex disaggregated rows to give totals", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2010, 2011, 2012, 2012),
                                sex = c('female', 'male', 'female', 'female', 'male'),
                                tot = c(100, 200, 400, 700, 100),
                                totpos = NA,
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    result <- combine_rows(test_prgm_dat)

    expect_identical(result$tot, c(300, NA, 800))
    expect_identical(result$year, c(2010, 2011, 2012))
})

test_that("combine_rows can combine sex disaggregated and aggregated rows to give totals", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2010, 2011, 2012),
                                sex = c('female', 'male', 'female', 'both'),
                                tot = c(100, 100, 100, 500),
                                totpos = c(100, 100, 100, 500),
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    result <- combine_rows(test_prgm_dat)

    expect_identical(result$tot, c(200, NA, 500))
    expect_identical(result$year, c(2010, 2011, 2012))
})

test_that("combine_rows discards aggregated rows if dis-aggregated are provided", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2010, 2010),
                                sex = c('both','female', 'male'),
                                tot = c(100, 100, 100),
                                totpos = NA,
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    result <- combine_rows(test_prgm_dat)

    expect_identical(result$tot, c(200))
    expect_identical(result$year, c(2010))
})

test_that("combine_rows aggregates totpos", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2010, 2011, 2011),
                                sex = c('female','male','female', 'male'),
                                tot = c(200, 200, 200, 200),
                                totpos = c(100, 100,100, 100),
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    result <- combine_rows(test_prgm_dat)

    expect_identical(result$tot, c(400, 400))
    expect_identical(result$totpos, c(200, 200))
    expect_identical(result$year, c(2010, 2011))
})

test_that("plot_input_tot does not error given NAs", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = 2010,
                                sex = 'both',
                                tot = NA,
                                totpos = NA,
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    plot_input_tot(test_prgm_dat, fp)
})

test_that("plot_input_tot does not error given only sex disaggregated data", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2010, 2011),
                                sex = c('female', 'male', 'female'),
                                tot = c(107634, 107634, 107634),
                                totpos = 50115,
                                vct = 107634,
                                vctpos = 25057,
                                anc = 107635,
                                ancpos = 25058)

    plot_input_tot(test_prgm_dat, fp)
})

test_that("plot_input_tot does not error given sex aggregated data", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = 2010,
                                sex = 'both',
                                tot = 215269,
                                totpos = 50115,
                                vct = 107634,
                                vctpos = 25057,
                                anc = 107635,
                                ancpos = 25058)

    plot_input_tot(test_prgm_dat, fp)
})

test_that("plot_input_tot does not error given mixed sex aggregated and sex dis-aggregated data", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = c(2010, 2011, 2011),
                                sex = c('both', 'female', 'male'),
                                tot = 215269,
                                totpos = 50115,
                                vct = 107634,
                                vctpos = 25057,
                                anc = 107635,
                                ancpos = 25058)

    plot_input_tot(test_prgm_dat, fp)
})

test_that("plot_input_tot does not error given incomplete data", {

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = 2010,
                                sex = 'female',
                                tot = 100,
                                totpos = 50,
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    plot_input_tot(test_prgm_dat, fp)

    test_prgm_dat <- data.frame(country = "Malawi",
                                year = 2010,
                                sex = 'male',
                                tot = 100,
                                totpos = 50,
                                vct = NA,
                                vctpos = NA,
                                anc = NA,
                                ancpos = NA)

    plot_input_tot(test_prgm_dat, fp)
})
