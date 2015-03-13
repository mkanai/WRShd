library(testthat)

context("Unit Tests for qcomhd")

set.seed(2)
x <- rnorm(1e2)
y <- rnorm(1e2)

test_that("Compare quantiles via the Harrell-Davis estimator", {
    expect_that(WRShd::qcomhd(x, y, plotit = FALSE), equals(WRS::qcomhd(x, y), plotit = FALSE))
    expect_that(WRShd::qcomhd(x, y, q = .75, nboot = 200, plotit = FALSE), equals(WRS::qcomhd(x, y, q = .75, nboot = 200, plotit = FALSE)))
    expect_that(WRShd::qcomhd(x, y, alpha = .1, nboot = 200, plotit = FALSE), equals(WRS::qcomhd(x, y, alpha = .1, nboot = 200, plotit = FALSE)))
    expect_that(WRShd::qcomhd(x, y, cores = 2, nboot = 200, plotit = FALSE), equals(WRS::qcomhd(x, y, nboot = 200, plotit = FALSE)))
})