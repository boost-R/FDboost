# Methods for objects of class bootstrapCI

Methods for objects that are fitted to compute bootstrap confidence
intervals.

## Usage

``` r
# S3 method for class 'bootstrapCI'
plot(
  x,
  which = NULL,
  pers = TRUE,
  commonRange = TRUE,
  showNumbers = FALSE,
  showQuantiles = TRUE,
  ask = TRUE,
  probs = c(0.25, 0.5, 0.75),
  ylim = NULL,
  ...
)

# S3 method for class 'bootstrapCI'
print(x, ...)
```

## Arguments

- x:

  an object of class `bootstrapCI`.

- which:

  base-learners that are plotted

- pers:

  plot coefficient surfaces as persp-plots? Defaults to `TRUE`.

- commonRange, :

  plot predicted coefficients on a common range, defaults to `TRUE`.

- showNumbers:

  show number of curve in plot of predicted coefficients, defaults to
  `FALSE`

- showQuantiles:

  plot the 0.05 and the 0.95 Quantile of coefficients in 1-dim effects.

- ask:

  defaults to `TRUE`, ask for next plot using `par(ask = ask)`?

- probs:

  vector of quantiles to be used in the plotting of 2-dimensional
  coefficients surfaces, defaults to `probs = c(0.25, 0.5, 0.75)`

- ylim:

  values for limits of y-axis

- ...:

  additional arguments passed to callies.

## Value

No return value (plot method) or `x` itself (print method)

## Details

`plot.bootstrapCI` plots the bootstrapped coefficients.
