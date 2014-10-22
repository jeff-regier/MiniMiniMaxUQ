require(testthat)
require(MiniMiniMaxUQ)

context("the borehole function")

p = 8

borehole <- function(x)
{
  rw <- x[1]
  r  <- x[2]
  Tu <- x[3]
  Hu <- x[4]
  Tl <- x[5]
  Hl <- x[6]
  L  <- x[7]
  Kw <- x[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

ranges = matrix(c(0.05, 0.15, 100, 50000, 63070, 115600, 990, 1110, 63.1, 116, 700, 820, 1120, 1680, 9855, 12045), nrow=2)

hypercube.to.domain = function(ww) {
  (ww * (ranges[2,] - ranges[1,])) + ranges[1,]
}

borehole.hypercube <- function(ww) {
  x = hypercube.to.domain(ww)
  return(borehole(x))
}

X = t(matrix(sapply(0:(2^p-1), function(i) int.to.bits(i, p)), nrow=p))
X = rbind(X, rep(0.5, 8))
X2 = matrix(0.5, nrow=16, ncol=8)
for (i in 1:8) {
  X2[i,i] = 0
  X2[8 + i, i] = 1
}
X = rbind(X, X2)
f.X = apply(X, 1, borehole.hypercube)

K.hat = find.K.hat(X, f.X)

test_that("lower bound on computational burden is trivial---0 or less", {
  expect_equal(max(lower.bound.computational.burden(X, f.X, 0.001, Inf, K.hat), 0), 0)
})

test_that("uncertainty is 0 at all corners, when all corners are in X",{ 
          expect_equal(corners.uncertainty.bound(X, f.X, K.hat), 0.)
})

max.uncertainty = branch.and.bound.max.uncertainty(X, f.X, K.hat)
test_that("bounds on max uncertainty surround the truth", {
  expect_true(lower.bound.max.uncertainty(X, f.X, K.hat) <= max.uncertainty)
  expect_true(upper.bound.max.uncertainty(X, f.X, K.hat) >= max.uncertainty)
})

cbs = uncertainty.confidence.bounds(X, f.X, K.hat, n=5000, confidence=.999)
test_that("confidence bounds are sensible", {
  expect_true(cbs[1,1] < 124)
  expect_true(124 < cbs[1,2])
  expect_true(124 < cbs[2,1])
  expect_true(124 < cbs[2,2])
  expect_true(cbs[1,2] < cbs[1,3])
  expect_true(cbs[2,2] < cbs[2,3])
})

