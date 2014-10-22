# Copyright 2014 Jeffrey Regier

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


sup.dist <- function(a, b) {
  lp_distance(a, b, -42.)
}

#' Finds the empirical Lipschitz constant of \eqn{f},
#' a lower bound on the Lipschitz constant of \eqn{f}.
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param lp.norm the \eqn{L_p} norm for computing the distances between points in \eqn{X}. Defaults to sup-norm.
#' @param min.distance only pairs of points in \eqn{X} farther apart than this threshold are considered. Defaults to 0.
#' 
find.K.hat <- function(X, f.X, lp.norm=Inf, min.distance=0) {
  stopifnot(max(X) <= 1. && min(X) >= 0.)
  lp.norm.cpp = if (lp.norm == Inf) -42. else lp.norm
  find_K_hat(as.matrix(X), f.X, lp.norm.cpp, min.distance)
}

#' Computes pointwise uncertainty
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param K.hat the empirical Lipschitz constant of \eqn{f}
#' @param x the point for which to compute pointwise uncertainty
#' 
pointwise.uncertainty <- function(X, f.X, K.hat, x) {
  pointwise_uncertainty(X, f.X, K.hat, x)
}

#------------------------------------------------------

compute.gamma.bar <-function(f.X, p, K.hat) {
  p.dim.area = function(y) (sum(abs((2 * (f.X - y) / K.hat)^p)))
  optimize(p.dim.area, c(min(f.X), max(f.X)))$minimum
}

unit.ball.volume <- function(lp.norm, p.dims) {
  stopifnot(lp.norm == Inf || lp.norm == 2)
  
  if (lp.norm == Inf)
    2^p.dims
  else if (lp.norm == 2)
    pi^(p.dims/2) / gamma(p.dims/2 + 1)
}

#' Computes a lower bound on the minimum computational burden.
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param epsilon the maximum uncertainty allowed
#' @param lp.norm the \eqn{L_p} norm for computing the distances between points in \eqn{[0, 1]^p}. Defaults to sup-norm.
#' @param K.hat the empirical Lipschitz constant of \eqn{f}. Defaults to find.K.hat(X, f.X).
#' 
lower.bound.computational.burden <- function(X, f.X, epsilon, lp.norm=Inf, K.hat=NULL) {
  stopifnot(max(X) <= 1. && min(X) >= 0.)
  stopifnot(epsilon >= 0.)
  stopifnot(lp.norm == Inf || lp.norm %% 1 == 0)
  stopifnot(is.null(K.hat) || K.hat >= 0.)
  
  if (is.null(K.hat))
    K.hat = find.K.hat(X, f.X, lp.norm)
  
  p = ncol(X)
  gamma.bar = compute.gamma.bar(f.X, p, K.hat)
  
  C = unit.ball.volume(lp.norm, p)
  ceiling(epsilon^-p * (K.hat^p/C - sum((abs(f.X - gamma.bar))^p)))
}

#------------------------------------------------------

int.to.bits <- function(i, p.dims) {
  rev(as.integer(intToBits(i))[1:p.dims])
}

#' Computes a lower bound on the maximum uncertainty by evaluating
#' the pointwise error at corners of the domain.
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param n the number of consecutive corners to consider. Defaults to 5000.
#' @param first.corner the first corner to consider. Defaults to 0.
#' @param K.hat the empirical Lipschitz constant of \eqn{f}. Defaults to \code{find.K.hat(X, f.X)}.
#' 
corners.uncertainty.bound <- function(X, f.X, n=5000, first.corner=0, K.hat=NULL) {
  stopifnot(is.null(K.hat) || K.hat >= 0.)
  
  if (is.null(K.hat))
    K.hat = find.K.hat(X, f.X, Inf)

  p.dims = ncol(X)
  max.uncertainty.lb = 0.
  for (i in first.corner:min(first.corner + n, (2^p.dims-1))) {
    bits = int.to.bits(i, p.dims)
    cur.lb = pointwise_uncertainty(X, f.X, K.hat, bits)
    if (cur.lb > max.uncertainty.lb) {
      cat(sprintf("maximum uncertainty so far at corner # %d: %.3f\n", i, cur.lb))
      max.uncertainty.lb = cur.lb
    }
  }
  
  max.uncertainty.lb
}

#------------------------------------------------------

random.point <- function(bottom.left, top.right) {
  edge.lengths = top.right - bottom.left
  offsets = runif(length(bottom.left))
  bottom.left + edge.lengths * offsets
}

mid.point.1d <- function(a, b) {
  (a + b) / 2
}

half.point <- function(a, b, dim.id) {
  ret = a
  ret[dim.id] = mid.point.1d(a[dim.id], b[dim.id])
  ret
}

divide.box <- function(box, dim.id) {
  sub.box.1 = rbind(box[1,], half.point(box[2,], box[1,], dim.id))
  sub.box.2 = rbind(half.point(box[1,], box[2,], dim.id), box[2,])
  list(sub.box.1, sub.box.2)
}

widest.divided.box <- function(box) {
  tol = .0001
  widths = box[2,] - box[1,]
  widest.dim = which.max(widths)
  if (widths[widest.dim] > tol)
    divide.box(box, widest.dim)
  else
    list(box)
}

#' Computes a lower bound on the maximum uncertainty within the specified box,
#' by evaluating the pointwise error at the box's corners and at random points
#' in the box's interior.
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param K.hat the empirical Lipschitz constant of \eqn{f}
#' @param box the top-left corner and the bottom-right corner of a hypercube contained in the domain. Defaults to \eqn{[0,1]^p}.
#' 
lower.bound.max.uncertainty <- function(X, f.X, K.hat, box=NULL) {
  if (is.null(box)) {
    p.dims = ncol(X)
    box = rbind(rep(0, p.dims), rep(1, p.dims))
  }
  
  num.lb.samples = 10
  bottom.left = box[1,]
  top.right = box[2,]
  get.rp = function() random.point(bottom.left, top.right)
  random.pts = replicate(num.lb.samples, get.rp())
  random.pts.mat = t(matrix(random.pts, ncol=num.lb.samples))
  all.pts = rbind(random.pts.mat, bottom.left, top.right)
  bnd.pointwise.uncertainty = function(x) pointwise_uncertainty(X, f.X, K.hat, x)
  lower.bounds = apply(all.pts, 1, bnd.pointwise.uncertainty)
  max(lower.bounds)
}

#' Computes an upper bound on the maximum uncertainty within the specified box.
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param K.hat the empirical Lipschitz constant of \eqn{f}
#' @param box the top-left corner and the bottom-right corner of a hypercube contained in the domain. Defaults to \eqn{[0,1]^p}.
#' 
upper.bound.max.uncertainty <- function(X, f.X, K.hat, box=NULL) {
  if (is.null(box)) {
    p.dims = ncol(X)
    box = rbind(rep(0, p.dims), rep(1, p.dims))
  }
  
  bottom.left = box[1,]
  top.right = box[2,]
  top.deviations = f_deviations(X, K.hat, top.right)
  bottom.deviations = f_deviations(X, K.hat, bottom.left)
  paired.deviations = rbind(bottom.deviations, top.deviations)
  greatest.deviations = apply(paired.deviations, 2, max)
  do_pointwise_uncertainty(f.X, greatest.deviations)
}

#' Computes the maximum uncertainty by a branch and bound algorithm.
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param K.hat the empirical Lipschitz constant of \eqn{f}
#' @param tol solutions this close to the maximum uncertainty are good enough. Defaults to 0.001.
#' 
branch.and.bound.max.uncertainty <- function(X, f.X, K.hat=NULL, tol=1e-3) {
  stopifnot(max(X) <= 1. && min(X) >= 0.)
  stopifnot(is.null(K.hat) || K.hat >= 0.)

  if (is.null(K.hat))
    K.hat = find.K.hat(X, f.X, Inf)

  p.dims = ncol(X)
  domain = rbind(rep(0, p.dims), rep(1, p.dims))
  candidates = list(domain)
  get.box.ub = function(box) upper.bound.max.uncertainty(X, f.X, K.hat, box)
  upper.bounds = c(get.box.ub(domain))
  max.lb = corners.uncertainty.bound(X, f.X, 10, 0, K.hat)

  prev.max.lb = -42.
  prev.max.ub = upper.bounds[1] + 42.

  while (length(candidates) != 0) {
    cur.box = candidates[[1]]
    candidates = candidates[-1]

    cur.ub = upper.bounds[[1]]
    upper.bounds = upper.bounds[-1]
    upper.bound.max.uncertainty(X, f.X, K.hat, cur.box)

    if (cur.ub > (max.lb + tol)) {
      cur.lb = lower.bound.max.uncertainty(X, f.X, K.hat, cur.box)
      if (cur.lb > max.lb) max.lb = cur.lb
      children = widest.divided.box(cur.box)
      candidates = c(candidates, children)
      upper.bounds = c(upper.bounds, sapply(children, get.box.ub))
    }
    
    if (length(upper.bounds) == 0)
      break
    
    max.ub = max(upper.bounds)
    
    if (max.ub != prev.max.ub || max.lb != prev.max.lb) {
      cat(sprintf("max uncertainty is in [%.3f, %.3f]; queue length is %d\n", 
                  max.lb, max.ub, length(candidates)))
    }
    
    if (max.lb > max.ub + tol)
      break

    prev.max.lb = max.lb
    prev.max.ub = max.ub
  }
  
  max.lb
}

#------------------------------------------------------

q.CI <- function(x, p, P) {
  # x is the sample, p (0<p<1) the quantile, P the confidence level
  x <- sort(x)
  n <- length(x)
  s <- min(which(pbinom((0:n),n,p) >= 1-(1-P)/2))
  r <- max(which(pbinom((0:n),n,p) <= (1-P)/2))
  c(x[r],x[s])
  # x[r] is the lower limit, x[s] the upper limit, of the CI
}

#' Computes confidence bounds for 4 statistics (first quartile, median, third quartile, mean)
#' of the distribution of pointwise uncertainties.
#' @param X a finite subset of \eqn{[0, 1]^p} where \eqn{f} is observed
#' @param f.X \eqn{f|_X}, the values of \eqn{f} on \eqn{X}
#' @param K.hat the empirical Lipschitz constant of \eqn{f}
#' @param n the number of samples. Defaults to 1000.
#' @param confidence the size of the confidence bounds. Defaults to 0.95.
#' 
uncertainty.confidence.bounds <- function(X, f.X, K.hat, n=1000, confidence=0.95) {
  p.dims = ncol(X)
  
  rand.pts = matrix(runif(p.dims * n), ncol=p.dims)
  wrapped.pointwise.uncertainty = function(x) pointwise_uncertainty(X, f.X, K.hat, x)
  errors = apply(rand.pts, 1, wrapped.pointwise.uncertainty)
  
  mu = mean(errors)
  confidence2 = 1 - ((1 - confidence) / 2)
  mu.error = qnorm(confidence2) * sd(errors) / sqrt(n)

  q1 = q.CI(errors, 1/4, confidence)
  q2 = q.CI(errors, 2/4, confidence)
  q3 = q.CI(errors, 3/4, confidence)
  mean.CI = c(mu - mu.error, mu + mu.error)

  df = as.data.frame(cbind(q1, q2, q3, mean.CI),
                     row.names=c("lower", "upper"))
  names(df) = c("1st Qu.", "Median", "3rd Qu.", "Mean")
  df
}
