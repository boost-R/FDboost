# Cross-validation for FDboostLSS

Multidimensional cross-validated estimation of the empirical risk for
hyper-parameter selection, for an object of class `FDboostLSS` setting
the folds per default to resampling curves.

## Usage

``` r
# S3 method for class 'FDboostLSS'
cvrisk(
  object,
  folds = cvLong(id = object[[1]]$id, weights = model.weights(object[[1]])),
  grid = NULL,
  papply = mclapply,
  trace = TRUE,
  fun = NULL,
  ...
)
```

## Arguments

- object:

  an object of class `FDboostLSS`.

- folds:

  a weight matrix a weight matrix with number of rows equal to the
  number of observations. The number of columns corresponds to the
  number of cross-validation runs, defaults to 25 bootstrap samples,
  resampling whole curves

- grid:

  defaults to a grid up to the current number of boosting iterations.
  The default generates the grid according to the defaults of
  [`cvrisk.mboostLSS`](https://rdrr.io/pkg/gamboostLSS/man/cvrisk.html)
  which are different for models with cyclic or noncyclic fitting.

- papply:

  (parallel) apply function, defaults to
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html), see
  [`cvrisk.mboostLSS`](https://rdrr.io/pkg/gamboostLSS/man/cvrisk.html)
  for details.

- trace:

  print status information during cross-validation? Defaults to `TRUE`.

- fun:

  if `fun` is `NULL`, the out-of-sample risk is returned. `fun`, as a
  function of `object`, may extract any other characteristic of the
  cross-validated models. These are returned as is.

- ...:

  additional arguments passed to
  [`mclapply`](https://rdrr.io/r/parallel/mclapply.html).

## Value

An object of class `cvriskLSS` (when `fun` was not specified), basically
a matrix containing estimates of the empirical risk for a varying number
of bootstrap iterations. `plot` and `print` methods are available as
well as an `mstop` method, see
[`cvrisk.mboostLSS`](https://rdrr.io/pkg/gamboostLSS/man/cvrisk.html).

## Details

The function `cvrisk.FDboostLSS` is a wrapper for `cvrisk.mboostLSS` in
package `gamboostLSS`. It overrides the default for the folds, so that
the folds are sampled on the level of curves (not on the level of single
observations, which does not make sense for functional response).

## See also

[`cvrisk.mboostLSS`](https://rdrr.io/pkg/gamboostLSS/man/cvrisk.html) in
package `gamboostLSS`.
