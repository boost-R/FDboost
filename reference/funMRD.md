# Functional MRD

Calculates the functional MRD for a fitted FDboost-object

## Usage

``` r
funMRD(object, overTime = TRUE, breaks = object$yind, global = FALSE, ...)
```

## Arguments

- object:

  fitted FDboost-object with regular response

- overTime:

  per default the functional MRD is calculated over time if
  `overTime=FALSE`, the MRD is calculated per curve

- breaks:

  an optional vector or number giving the time-points at which the model
  is evaluated. Can be specified as number of equidistant time-points or
  as vector of time-points. Defaults to the index of the response in the
  model.

- global:

  logical. defaults to `FALSE`, if TRUE the global MRD like in a normal
  linear model is calculated

- ...:

  currently not used

## Value

Returns a vector with the calculated MRD and some extra information in
attributes.

## Details

Formula to calculate MRD over time, `overTime=TRUE`:  
\\ MRD(t) = n^{-1} \sum_i \|Y_i(t) - \hat{Y}\_i(t)\| / \|Y_i(t)\| \\

Formula to calculate MRD over subjects, `overTime=FALSE`:  
\\ MRD\_{i} = \int \|Y_i(t) - \hat{Y}\_i(t)\| / \|Y_i(t)\| dt \approx
G^{-1} \sum_g \|Y_i(t_g) - \hat{Y}\_i(t_g)\| / \|Y_i(t)\|\\

## Note

`breaks` cannot be changed in the case the
[`bsignal()`](https://boost-R.github.io/FDboost/reference/bsignal.md) is
used over the same domain as the response! In that case you would have
to rename the index of the response or that of the covariates.
