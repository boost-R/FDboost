# Prediction for boosted functional regression model

Takes a fitted `FDboost`-object produced by
[`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md)
and produces predictions given a new set of values for the model
covariates or the original values used for the model fit. This is a
wrapper function for
[`predict.mboost()`](https://rdrr.io/pkg/mboost/man/methods.html)

## Usage

``` r
# S3 method for class 'FDboost'
predict(object, newdata = NULL, which = NULL, toFDboost = TRUE, ...)
```

## Arguments

- object:

  a fitted `FDboost`-object

- newdata:

  a named list or a data frame containing the values of the model
  covariates at which predictions are required. If this is not provided
  then predictions corresponding to the original data are returned. If
  `newdata` is provided then it should contain all the variables needed
  for prediction, in the format supplied to `FDboost`, i.e., functional
  predictors must be supplied as matrices with each row corresponding to
  one observed function.

- which:

  a subset of base-learners to take into account for computing
  predictions or coefficients. If which is given (as an integer vector
  corresponding to base-learners) a list is returned.

- toFDboost:

  logical, defaults to `TRUE`. In case of regular response in wide
  format (i.e. response is supplied as matrix): should the predictions
  be returned as matrix, or list of matrices instead of vectors

- ...:

  additional arguments passed on to
  [`predict.mboost()`](https://rdrr.io/pkg/mboost/man/methods.html).

## Value

a matrix or list of predictions depending on values of unlist and which

## See also

[`FDboost`](https://boost-R.github.io/FDboost/reference/FDboost.md) for
the model fit and
[`plotPredicted`](https://boost-R.github.io/FDboost/reference/plot.FDboost.md)
for a plot of the observed values and their predictions.
