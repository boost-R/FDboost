# Base-learners for Functional Covariates

Base-learners that fit historical functional effects that can be used
with the tensor product, as, e.g., `hbistx(...) %X% bolsc(...)`, to form
interaction effects (Ruegamer et al., 2018). For expert use only! May
show unexpected behavior compared to other base-learners for functional
data!

## Usage

``` r
bhistx(
  x,
  limits = "s<=t",
  standard = c("no", "time", "length"),
  intFun = integrationWeightsLeft,
  inS = c("smooth", "linear", "constant"),
  inTime = c("smooth", "linear", "constant"),
  knots = 10,
  boundary.knots = NULL,
  degree = 3,
  differences = 1,
  df = 4,
  lambda = NULL,
  penalty = c("ps", "pss"),
  check.ident = FALSE
)
```

## Arguments

- x:

  object of type `hmatrix` containing time, index and functional
  covariate; note that `timeLab` in the `hmatrix`-object must be equal
  to the name of the time-variable in `timeformula` in the
  `FDboost`-call

- limits:

  defaults to `"s<=t"` for an historical effect with s\<=t; either one
  of `"s<t"` or `"s<=t"` for \[l(t), u(t)\] = \[T1, t\]; otherwise
  specify limits as a function for integration limits \[l(t), u(t)\]:
  function that takes \\s\\ as the first and `t` as the second argument
  and returns `TRUE` for combinations of values (s,t) if \\s\\ falls
  into the integration range for the given \\t\\.

- standard:

  the historical effect can be standardized with a factor. "no" means no
  standardization, "time" standardizes with the current value of time
  and "lenght" standardizes with the lenght of the integral

- intFun:

  specify the function that is used to compute integration weights in
  `s` over the functional covariate \\x(s)\\

- inS:

  historical effect can be smooth, linear or constant in s, which is the
  index of the functional covariates x(s).

- inTime:

  historical effect can be smooth, linear or constant in time, which is
  the index of the functional response y(time).

- knots:

  either the number of knots or a vector of the positions of the
  interior knots (for more details see
  [`bbs`](https://rdrr.io/pkg/mboost/man/baselearners.html)`)`.

- boundary.knots:

  boundary points at which to anchor the B-spline basis (default the
  range of the data). A vector (of length 2) for the lower and the upper
  boundary knot can be specified.

- degree:

  degree of the regression spline.

- differences:

  a non-negative integer, typically 1, 2 or 3. Defaults to 1. If
  `differences` = *k*, *k*-th-order differences are used as a penalty
  (*0*-th order differences specify a ridge penalty).

- df:

  trace of the hat matrix for the base-learner defining the base-learner
  complexity. Low values of `df` correspond to a large amount of
  smoothing and thus to "weaker" base-learners.

- lambda:

  smoothing parameter of the penalty, computed from `df` when `df` is
  specified.

- penalty:

  by default, `penalty="ps"`, the difference penalty for P-splines is
  used, for `penalty="pss"` the penalty matrix is transformed to have
  full rank, so called shrinkage approach by Marra and Wood (2011)

- check.ident:

  use checks for identifiability of the effect, based on Scheipl and
  Greven (2016); see Brockhaus et al. (2017) for identifiability checks
  that take into account the integration limits

## Value

Equally to the base-learners of package mboost:

An object of class `blg` (base-learner generator) with a `dpp` function
(dpp, data pre-processing).

The call of `dpp` returns an object of class `bl` (base-learner) with a
`fit` function. The call to `fit` finally returns an object of class
`bm` (base-model).

## Details

`bhistx` implements a base-learner for functional covariates with
flexible integration limits `l(t)`, `r(t)` and the possibility to
standardize the effect by `1/t` or the length of the integration
interval. The effect is `stand * int_{l(t)}^{r_{t}} x(s)beta(t,s) ds`.
The base-learner defaults to a historical effect of the form
\\\int\_{T1}^{t} x_i(s)beta(t,s) ds\\, where \\T1\\ is the minimal index
of \\t\\ of the response \\Y(t)\\. `bhistx` can only be used if \\Y(t)\\
and \\x(s)\\ are observd over the same domain \\s,t \in \[T1, T2\]\\.
The base-learner `bhistx` can be used to set up complex interaction
effects like factor-specific historical effects as discussed in Ruegamer
et al. (2018).

Note that the data has to be supplied as a `hmatrix` object for model
fit and predictions.

## References

Brockhaus, S., Melcher, M., Leisch, F. and Greven, S. (2017): Boosting
flexible functional regression models with a high number of functional
historical effects, Statistics and Computing, 27(4), 913-926.

Marra, G. and Wood, S.N. (2011): Practical variable selection for
generalized additive models. Computational Statistics & Data Analysis,
55, 2372-2387.

Ruegamer D., Brockhaus, S., Gentsch K., Scherer, K., Greven, S. (2018).
Boosting factor-specific functional historical models for the detection
of synchronization in bioelectrical signals. Journal of the Royal
Statistical Society: Series C (Applied Statistics), 67, 621-642.

Scheipl, F., Staicu, A.-M. and Greven, S. (2015): Functional Additive
Mixed Models, Journal of Computational and Graphical Statistics, 24(2),
477-501. <https://arxiv.org/abs/1207.5947>

Scheipl, F. and Greven, S. (2016): Identifiability in penalized
function-on-function regression models. Electronic Journal of
Statistics, 10(1), 495-526.

## See also

[`FDboost`](https://boost-R.github.io/FDboost/reference/FDboost.md) for
the model fit and
[`bhist`](https://boost-R.github.io/FDboost/reference/bsignal.md) for
simple hisotorical effects.

## Examples

``` r
if(require(refund)){
## simulate some data from a historical model
## the interaction effect is in this case not necessary
n <- 100
nygrid <- 35
data1 <- suppressWarnings(pffrSim(scenario = c("int", "ff"), limits = function(s,t){ s <= t }, 
                n = n, nygrid = nygrid))
data1$X1 <- scale(data1$X1, scale = FALSE) ## center functional covariate                  
dataList <- as.list(data1)
dataList$tvals <- attr(data1, "yindex")

## create the hmatrix-object
X1h <- with(dataList, hmatrix(time = rep(tvals, each = n), id = rep(1:n, nygrid), 
                             x = X1, argvals = attr(data1, "xindex"), 
                             timeLab = "tvals", idLab = "wideIndex", 
                             xLab = "myX", argvalsLab = "svals"))
dataList$X1h <- I(X1h)   
dataList$svals <- attr(data1, "xindex")
## add a factor variable 
dataList$zlong <- factor(gl(n = 2, k = n/2, length = n*nygrid), levels = 1:2)  
dataList$z <- factor(gl(n = 2, k = n/2, length = n), levels = 1:2)

## do the model fit with main effect of bhistx() and interaction of bhistx() and bolsc()
mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) + 
               bhistx(x = X1h, df = 5, knots = 5) %X% bolsc(zlong), 
              timeformula = ~ bbs(tvals, knots = 10), data = dataList)
              
## alternative parameterization: interaction of bhistx() and bols()
mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) %X% bols(zlong), 
              timeformula = ~ bbs(tvals, knots = 10), data = dataList)

# \donttest{
  # find the optimal mstop over 5-fold bootstrap (small example to reduce run time)
  cv <- cvrisk(mod, folds = cv(model.weights(mod), B = 5))
  mstop(cv)
  mod[mstop(cv)]
  
  appl1 <- applyFolds(mod, folds = cv(rep(1, length(unique(mod$id))), type = "bootstrap", B = 5))

 # plot(mod)
# }
}
#> Loading required package: refund
#> Use a smooth offset.
#> Warning: matrix is structurally rank deficient; using augmented matrix with additional 15 row(s) of zeros
#> Warning: matrix is structurally rank deficient; using augmented matrix with additional 15 row(s) of zeros
#> Use a smooth offset.
#> Warning: matrix is structurally rank deficient; using augmented matrix with additional 30 row(s) of zeros
#> The smooth offset is fixed over all folds.
```
