context("pjn-country")

test_that("pjn country is ready from country code correctly", {
  ## Added as a regression check for reported issue due to 
  ## errors in name of São Tomé and Príncipe
  expect_equal(get_pjn_country(pjn_data), "São Tomé and Príncipe")
})

test_that("iso3 cane be retrieved from pjn file", {
  expect_equal(get_pjn_iso3(pjn_data) , 678)
})