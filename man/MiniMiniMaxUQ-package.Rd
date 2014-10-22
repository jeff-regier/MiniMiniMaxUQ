\name{MiniMiniMaxUQ-package}
\alias{MiniMiniMaxUQ-package}
\alias{MiniMiniMaxUQ}
\docType{package}
\title{
Mini-Minimax Uncertainty Quantification for Emulators
}
\description{
MiniMiniMaxUQ seeks to answer the following: How
    accurately can a black-box function f be emulated from a given set of data?
    How many evaluations of f are required to guarantee that f can be emulated
    to a given level of accuracy?
}
\details{
\tabular{ll}{
Package: \tab MiniMiniMaxUQ\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-10-20\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Jeffrey Regier and Philip Stark

Maintainer: Jeffrey Regier \email{jeff@@stat.berkeley.edu}
}
\references{
Jeffrey Regier and Philip Stark. "Mini-Minimax Uncertainty Quantification for Emulators." arXiv preprint arXiv:1303.3079 (2013).
}
\keyword{ package }
\examples{
library(MiniMiniMaxUQ)

X = expand.grid(1:9, 1:9) / 10
f <- function(x) sin(x[1]) + cos(x[2])
f.X = apply(X, 1, f)

K.hat = find.K.hat(X, f.X)

pointwise.uncertainty(X, f.X, K.hat, c(.55,.33))

lower.bound.computational.burden(X, f.X, K.hat, epsilon=.1)
corners.uncertainty.bound(X, f.X, K.hat)
lower.bound.max.uncertainty(X, f.X, K.hat)
upper.bound.max.uncertainty(X, f.X, K.hat)
branch.and.bound.max.uncertainty(X, f.X, K.hat)

uncertainty.confidence.bounds(X, f.X, K.hat)
}
