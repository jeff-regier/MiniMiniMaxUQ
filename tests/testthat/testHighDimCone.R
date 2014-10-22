require(testthat)
require(MiniMiniMaxUQ)

context("the 21-dimensional cone example")

p = 21

z = rep(.5, p)

f <- function(x) {
  max(abs(x - z))
}

X = matrix(0.5, nrow=2*p, ncol=p)
for (i in 1:p) {
  X[i,i] = 0
  X[p + i, i] = 1
}
X = rbind(rep(0.5, p), X)
f.X = apply(X, 1, f)

test_that("empirical Lipschitz constant is correct with various distance metrics", {
  expect_equal(find.K.hat(X, f.X, Inf), 1)
  expect_equal(find.K.hat(X, f.X, 2), 1)
  expect_equal(find.K.hat(X, f.X, 1), 1)
})

K.hat = 1

test_that("pointwise uncertainty meets established conditions", {
  expect_equal(pointwise.uncertainty(X, f.X, K.hat, z), 0)
  expect_equal(pointwise.uncertainty(X, f.X, K.hat, rep(1, p)), .25)
  expect_true(pointwise.uncertainty(X, f.X, K.hat, runif(p)) < .250001)
  expect_true(pointwise.uncertainty(X, f.X, K.hat, runif(p)) < .250001)
  expect_true(pointwise.uncertainty(X, f.X, K.hat, runif(p)) < .250001)
  expect_true(pointwise.uncertainty(X, f.X, K.hat, runif(p)) < .250001)
})

test_that("bounds on max uncertainty surround the true value", {
  expect_true(corners.uncertainty.bound(X, f.X, K.hat) < .250001)
  expect_true(lower.bound.max.uncertainty(X, f.X, K.hat) < .25001)
  expect_true(upper.bound.max.uncertainty(X, f.X, K.hat) > 0.24999)
  expect_equal(branch.and.bound.max.uncertainty(X, f.X, K.hat), 0.25)
})

