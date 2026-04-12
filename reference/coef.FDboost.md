# Coefficients of boosted functional regression model

Takes a fitted `FDboost`-object produced by
[`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md)
and returns estimated coefficient functions/surfaces \\\beta(t),
\beta(s,t)\\ and estimated smooth effects \\f(z), f(x,z)\\ or \\f(x, z,
t)\\. Not implemented for smooths in more than 3 dimensions.

## Usage

``` r
# S3 method for class 'FDboost'
coef(
  object,
  raw = FALSE,
  which = NULL,
  computeCoef = TRUE,
  returnData = FALSE,
  n1 = 40,
  n2 = 40,
  n3 = 20,
  n4 = 10,
  ...
)
```

## Arguments

- object:

  a fitted `FDboost`-object

- raw:

  logical defaults to `FALSE`. If `raw = FALSE` for each effect the
  estimated function/surface is calculated. If `raw = TRUE` the
  coefficients of the model are returned.

- which:

  a subset of base-learners for which the coefficients should be
  computed (numeric vector), defaults to NULL which is the same as
  `which=seq_along(object$baselearner)`. In the special case of
  `which=0`, only the coefficients of the offset are returned.

- computeCoef:

  defaults to `TRUE`, if `FALSE` only the names of the terms are
  returned

- returnData:

  return the dataset which is used to get the coefficient estimates as
  predictions, see Details.

- n1:

  see below

- n2:

  see below

- n3:

  n1, n2, n3 give the number of grid-points for 1-/2-/3-dimensional
  smooth terms used in the marginal equidistant grids over the range of
  the covariates at which the estimated effects are evaluated.

- n4:

  gives the number of points for the third dimension in a 3-dimensional
  smooth term

- ...:

  other arguments, not used.

## Value

If `raw = FALSE`, a list containing

- `offset` a list with plot information for the offset.

- `smterms` a named list with one entry for each smooth term in the
  model. Each entry contains

  - `x, y, z` the unique grid-points used to evaluate the
    smooth/coefficient function/coefficient surface

  - `xlim, ylim, zlim` the extent of the x/y/z-axes

  - `xlab, ylab, zlab` the names of the covariates for the x/y/z-axes

  - `value` a vector/matrix/list of matrices containing the coefficient
    values

  - `dim` the dimensionality of the effect

  - `main` the label of the smooth term (a short label)

If `raw = TRUE`, a list containing the estimated spline coefficients.

## Details

If `raw = FALSE` the function `coef.FDboost` generates adequate dummy
data and uses the function `predict.FDboost` to compute the estimated
coefficient functions.
