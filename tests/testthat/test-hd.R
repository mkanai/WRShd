library(testthat)

context("Unit Tests for hd")

set.seed(2)
x <- rnorm(1e3)

test_that("Compute the Harrell-Davis estimate of the qth quantile", {
          expect_that(WRShd::hd(x), equals(WRS::hd(x)))
          expect_that(WRShd::hd(x, q = .75), equals(WRS::hd(x, q = .75)))
          expect_that(WRShd::hd(x, na.rm = TRUE), equals(WRS::hd(x, na.rm = TRUE)))
          expect_that(WRShd::hd(x, cores = 2), equals(WRS::hd(x)))
         })