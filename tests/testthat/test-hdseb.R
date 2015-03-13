library(testthat)

context("Unit Tests for hdseb")

set.seed(2)
x <- rnorm(1e3)

test_that("Compute a a bootstrap standard error of the Harrell-Davis estimate of the qth quantile", {
    expect_that(WRShd::hdseb(x), equals(WRS::hdseb(x)))
    expect_that(WRShd::hdseb(x, q = .75), equals(WRS::hdseb(x, q = .75)))
    expect_that(WRShd::hdseb(x, nboot = 200), equals(WRS::hdseb(x, nboot = 200)))
    expect_that(WRShd::hdseb(x, cores = 2), equals(WRS::hdseb(x)))
})