context("test-simmod.R")

pjnz <- "sample_files/Malawi_2018_version_8.PJNZ"
filePath <- normalizePath(pjnz)

fp <- prepare_inputs(filePath)
mod <- simmod(fp)
modR <- simmod(fp, "R")

test_that("population outputs match target population", {
  expect_equal(c(fp$target_hivn_pop), c(mod[,,1,]))
  expect_equal(c(fp$target_hivp_pop), c(mod[,,2,]))
  expect_equal(c(fp$target_hivpop_ha), c(colSums(attr(mod, "hivpop"), dims=1)))
  expect_equal(c(fp$target_artpop_ha), c(colSums(attr(mod, "artpop"), dims=2)))
})

test_prev <- c(0, 1e-05, 2e-05, 4e-05, 6e-05, 0.00012, 0.00019, 3e-04, 0.00047,
               0.00073, 0.00113, 0.00175, 0.00273, 0.00428, 0.00671, 0.01046,
               0.01566, 0.02302, 0.03285, 0.04535, 0.06027, 0.07676, 0.09359,
               0.1094, 0.12308, 0.134, 0.14168, 0.14668, 0.14895, 0.14887, 0.1469,
               0.14342, 0.13876, 0.13335, 0.12769, 0.1223, 0.1175, 0.11383,
               0.11128, 0.10968, 0.1085, 0.10737, 0.1061, 0.10472, 0.10304,
               0.101, 0.09852, 0.0957, 0.09275, 0.0897, 0.08664, 0.08338, 0.07997)

test_incid <- c(0, 1e-05, 1e-05, 2e-05, 3e-05, 7e-05, 8e-05, 0.00013, 0.00019,
                3e-04, 0.00046, 0.00072, 0.00113, 0.00179, 0.0028, 0.00433, 0.00632,
                0.0091, 0.01241, 0.01601, 0.01944, 0.02219, 0.02392, 0.02453,
                0.02413, 0.02301, 0.02117, 0.01942, 0.01731, 0.01538, 0.01377,
                0.01255, 0.01167, 0.01107, 0.01073, 0.0105, 0.00995, 0.00971,
                0.00941, 0.00903, 0.00868, 0.00824, 0.00759, 0.00686, 0.00622,
                0.00561, 0.00501, 0.00436, 0.00396, 0.00362, 0.00339, 0.0031,
                0.00276)

test_artbycd4 <- structure(c(88242.22074, 86889.8681, 762756.74522, 106907.71164,
                             105605.23361, 1037838.36237, 106500.96704, 102268.55893, 1249894.45165,
                             105585.7832, 99297.93153, 1447553.7935, 103748.931, 97432.02137,
                             1434061.82078, 43443.70138, 40440.27621, 595785.36283, 20111.8456,
                             18603.62889, 300275.28866), .Dim = c(3L, 7L))

test_that("epidemic matches", {
  expect_equal(test_prev, round(attr(mod, "prev15to49"), 5))
  expect_equal(test_incid, round(attr(mod, "incid15to49"), 5))
  expect_equal(test_artbycd4, round(rowSums(attr(mod, "artpop"), dims=2), 5))
})
