# Residual values of a boosted functional regression model

Takes a fitted `FDboost`-object and computes the residuals, more
precisely the current value of the negative gradient is returned.

## Usage

``` r
# S3 method for class 'FDboost'
residuals(object, ...)
```

## Arguments

- object:

  a fitted `FDboost`-object

- ...:

  not used

## Value

matrix of residual values

## Details

The residual is missing if the corresponding value of the response was
missing.

## See also

[`FDboost`](https://boost-R.github.io/FDboost/reference/FDboost.md) for
the model fit.
