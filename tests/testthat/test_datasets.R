context("test_datasets")

data(survey_hts)

test_that("survey dataset has correct factor levels", {
  expect_equal(levels(survey_hts$sex), c("both", "female", "male"))
})

test_that("no non-ASCII characters in factor levels", {
  expect_equal(length(tools::showNonASCII(levels(survey_hts$country))), 0)
})
