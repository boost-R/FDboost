# Functional MSE

Calculates the functional MSE for a fitted FDboost-object

## Usage

``` r
funMSE(
  object,
  overTime = TRUE,
  breaks = object$yind,
  global = FALSE,
  relative = FALSE,
  root = FALSE,
  ...
)
```

## Arguments

- object:

  fitted FDboost-object

- overTime:

  per default the functional R-squared is calculated over time if
  `overTime=FALSE`, the R-squared is calculated per curve

- breaks:

  an optional vector or number giving the time-points at which the model
  is evaluated. Can be specified as number of equidistant time-points or
  as vector of time-points. Defaults to the index of the response in the
  model.

- global:

  logical. defaults to `FALSE`, if TRUE the global R-squared like in a
  normal linear model is calculated

- relative:

  logical. defaults to `FALSE`. If `TRUE` the MSE is standardized by the
  global variance of the response  
  \\ n^{-1} \int \sum_i (Y_i(t) - \bar{Y})^2 dt \approx G^{-1} n^{-1}
  \sum_g \sum_i (Y_i(t_g) - \bar{Y})^2 \\

- root:

  take the square root of the MSE

- ...:

  currently not used

## Value

Returns a vector with the calculated MSE and some extra information in
attributes.

## Details

Formula to calculate MSE over time, `overTime=TRUE`:  
\\ MSE(t) = n^{-1} \sum_i (Y_i(t) - \hat{Y}\_i(t))^2 \\

Formula to calculate MSE over subjects, `overTime=FALSE`:  
\\ MSE_i = \int (Y_i(t) - \hat{Y}\_i(t))^2 dt \approx G^{-1} \sum_g
(Y_i(t_g) - \hat{Y}\_i(t_g))^2\\

## Note

`breaks` cannot be changed in the case the
[`bsignal()`](https://boost-R.github.io/FDboost/reference/bsignal.md) is
used over the same domain as the response! In that case you would have
to rename the index of the response or that of the covariates.
