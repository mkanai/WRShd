library(testthat)

context("Unit Tests for hdpb")

set.seed(2)
x <- rep(rnorm(1e2), 2)

test_that("Compute a bootstrap 1-alpha confidence for the Harrell-Davis estimate of the qth quantile", {
    expect_that(WRShd::hdpb(x), equals(WRS::hdpb(x)))
    expect_that(WRShd::hdpb(x, q = .75, nboot = 200), equals(WRS::hdpb(x, q = .75, nboot = 200)))
    expect_that(WRShd::hdpb(x, alpha = .1, nboot = 200), equals(WRS::hdpb(x, alpha = .1, nboot = 200)))
    expect_that(WRShd::hdpb(x, cores = 2, nboot = 200), equals(WRS::hdpb(x, nboot = 200)))
})