### Mini-minimax uncertainty quantification for emulators
=============

[![Build Status](https://travis-ci.org/jeff-regier/MiniMiniMaxUQ.svg?branch=master)](https://travis-ci.org/jeff-regier/MiniMiniMaxUQ)
[![Coverage Status](https://coveralls.io/repos/jeff-regier/MiniMiniMaxUQ/badge.svg?branch=master&service=github)](https://coveralls.io/github/jeff-regier/MiniMiniMaxUQ?branch=master)

=============

The `MiniMiniMaxUQ` R package implements the optimization and bounding procedures described in

> [Jeffrey Regier and Philip Stark. "Mini-Minimax Uncertainty Quantification for Emulators." arXiv preprint arXiv:1303.3079 (2013).](http://arxiv.org/pdf/1303.3079.pdf)

#### Installation

To get the current development version from github:

```R
# install.packages("devtools")
devtools::install_github("jeff-regier/MiniMiniMaxUQ")
```

#### Example

The code below invokes each function exported by `MiniMiniMaxUQ`:
```R
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
```
Each function is further described in the man pages. 
Additional examples appear in `tests/testthat`.


#### License

`MiniMiniMaxUQ` is free software, licensed under GPLv3.
