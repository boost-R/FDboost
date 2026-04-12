# Functional R-squared

Calculates the functional R-squared for a fitted FDboost-object

## Usage

``` r
funRsquared(object, overTime = TRUE, breaks = object$yind, global = FALSE, ...)
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

- ...:

  currently not used

## Value

Returns a vector with the calculated R-squared and some extra
information in attributes.

## Details

`breaks` should be set to some grid, if there are many missing values or
time-points with very few observations in the dataset. Otherwise at
these points of t the variance will be almost 0 (or even 0 if there is
only one observation at a time-point), and then the prediction by the
local means \\\mu(t)\\ is locally very good. The observations are
interpolated linearly if necessary.

Formula to calculate R-squared over time, `overTime=TRUE`:  
\\R^2(t) = 1 - \sum\_{i}( Y_i(t) - \hat{Y}\_i(t))^2 / \sum\_{i}(
Y_i(t) - \bar{Y}(t) )^2 \\

Formula to calculate R-squared over subjects, `overTime=FALSE`:  
\\R^2_i = 1 - \int (Y_i(t) - \hat{Y}\_i(t))^2 dt / \int (Y_i(t) -
\bar{Y}\_i )^2 dt \\

## Note

`breaks` cannot be changed in the case the
[`bsignal()`](https://boost-R.github.io/FDboost/reference/bsignal.md) is
used over the same domain as the response! In that case you would have
to rename the index of the response or that of the covariates.

## References

Ramsay, J., Silverman, B. (2006). Functional data analysis. Wiley Online
Library. chapter 16.3
