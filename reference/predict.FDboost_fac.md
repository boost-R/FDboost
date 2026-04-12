# Prediction and plotting for factorized FDboost model components

Prediction and plotting for factorized FDboost model components

## Usage

``` r
# S3 method for class 'FDboost_fac'
predict(object, newdata = NULL, which = NULL, ...)

# S3 method for class 'FDboost_fac'
plot(x, which = NULL, main = NULL, ...)
```

## Arguments

- object, x:

  a model-factor given as a `FDboost_fac` object

- newdata:

  optionally, a data frame or list in which to look for variables with
  which to predict. See
  [`predict.mboost`](https://rdrr.io/pkg/mboost/man/methods.html).

- which:

  a subset of base-learner components to take into account for computing
  predictions or coefficients. Different components are never aggregated
  to a joint prediction, but always returned as a matrix or list. Select
  the k-th component by name in the format `bl(x, ...)[k]` or all
  components of a base-learner by dropping the index or all
  base-learners of a variable by using the variable name.

- ...:

  additional arguments passed to underlying methods.

- main:

  the plot title. By default, base-learner names are used with component
  numbers `[k]`.

## Value

A matrix of predictions (for predict method) or no return value (plot
method)

## See also

\[factorize(), factorize.FDboost()\]
