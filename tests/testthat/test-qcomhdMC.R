library(testthat)

context("Unit Tests for qcomhdMC")

set.seed(2)
x <- rnorm(1e2)
y <- rnorm(1e2)

test_that("Compare quantiles via the Harrell-Davis estimator", {
    expect_that(WRShd::qcomhdMC(x, y, plotit = FALSE), equals(WRS::qcomhdMC(x, y, plotit = FALSE)))
})