require(testthat)
require(MiniMiniMaxUQ)

context("the 21-dimensional cone example")

num.dims = 21

z = rep(.5, num.dims)

f <- function(x) {
  max(abs(x - z))
}

xs = rep(0.5, num.dims)
xs2 = matrix(0.5, nrow=2*num.dims, ncol=num.dims)
for (i in 1:num.dims) {
  xs2[i,i] = 0
  xs2[num.dims + i, i] = 1
}
xs = rbind(xs, xs2)
ys = apply(xs, 1, f)

test_that("empirical Lipschitz constant is correct with various distance metrics", {
  expect_equal(find.K.hat(xs, ys, Inf), 1)
  expect_equal(find.K.hat(xs, ys, 2), 1)
  expect_equal(find.K.hat(xs, ys, 1), 1)
})

K.hat = 1

test_that("pointwise uncertainty meets established conditions", {
  expect_equal(pointwise.uncertainty(xs, ys, K.hat, z), 0)
  expect_equal(pointwise.uncertainty(xs, ys, K.hat, rep(1, 21)), .25)
  expect_true(pointwise.uncertainty(xs, ys, K.hat, runif(21)) < .250001)
  expect_true(pointwise.uncertainty(xs, ys, K.hat, runif(21)) < .250001)
  expect_true(pointwise.uncertainty(xs, ys, K.hat, runif(21)) < .250001)
  expect_true(pointwise.uncertainty(xs, ys, K.hat, runif(21)) < .250001)
})

test_that("bounds on max uncertainty surround the true value", {
  expect_true(corners.uncertainty.bound(xs, ys, K.hat) < .250001)
  expect_true(lower.bound.max.uncertainty(xs, ys, K.hat) < .25001)
  expect_true(upper.bound.max.uncertainty(xs, ys, K.hat) > 0.24999)
  expect_equal(branch.and.bound.max.uncertainty(xs, ys, K.hat), 0.25)
})

