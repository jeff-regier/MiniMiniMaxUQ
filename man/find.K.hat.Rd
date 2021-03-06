% Generated by roxygen2 (4.0.2): do not edit by hand
\name{find.K.hat}
\alias{find.K.hat}
\title{Finds \eqn{\hat{K},}{} the empirical Lipschitz constant of \eqn{f},
a lower bound on the Lipschitz constant of \eqn{f}.}
\usage{
find.K.hat(X, f.X, lp.norm = Inf, min.distance = 0)
}
\arguments{
\item{X}{\eqn{X,}{} a matrix whose rows are points \eqn{\textrm{in } [0, 1]^p}{} where \eqn{f} is observed}

\item{f.X}{\eqn{f|_X,}{} a vector whose values are \eqn{f} evaluated at each row of \eqn{X}}

\item{lp.norm}{an integer (or \code{Inf}) specifying the \eqn{L_p} norm for computing the distances between points in \eqn{X}. Defaults to sup-norm.}

\item{min.distance}{only pairs of points in \eqn{X} farther apart than this threshold are considered. Defaults to 0.}
}
\description{
Finds \eqn{\hat{K},}{} the empirical Lipschitz constant of \eqn{f},
a lower bound on the Lipschitz constant of \eqn{f}.
}
\examples{
X <- expand.grid(1:9, 1:9) / 10
f <- function(x) sin(x[1]) + cos(x[2])
f.X <- apply(X, 1, f)
find.K.hat(X, f.X)
}

