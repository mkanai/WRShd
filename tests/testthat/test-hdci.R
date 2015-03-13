library(testthat)

context("Unit Tests for hdci")

set.seed(2)
x <- rnorm(1e3)

test_that("Compute a 1-alpha confidence for the Harrell-Davis estimate of the qth quantile", {
    expect_that(WRShd::hdci(x), equals(WRS::hdci(x)))
    expect_that(WRShd::hdci(x, q = .75), equals(WRS::hdci(x, q = .75)))
    expect_that(do.call(WRShd::hdci, list(x = x, alpha = .1)), throws_error("Use the function qcipb. Generally works well even when alpha is not equal to .05"))
    expect_that(WRShd::hdci(x, nboot = 200), equals(WRS::hdci(x, nboot = 200)))
    expect_that(WRShd::hdci(x, cores = 2), equals(WRS::hdci(x)))
})