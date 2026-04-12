# Fitted values of a boosted functional regression model

Takes a fitted `FDboost`-object and computes the fitted values.

## Usage

``` r
# S3 method for class 'FDboost'
fitted(object, toFDboost = TRUE, ...)
```

## Arguments

- object:

  a fitted `FDboost`-object

- toFDboost:

  logical, defaults to `TRUE`. In case of regular response in wide
  format (i.e., response is supplied as matrix): should the predictions
  be returned as matrix, or list of matrices instead of vectors

- ...:

  additional arguments passed on to
  [`predict.FDboost`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)

## Value

matrix or vector of fitted values

## See also

[`FDboost`](https://boost-R.github.io/FDboost/reference/FDboost.md) for
the model fit.
