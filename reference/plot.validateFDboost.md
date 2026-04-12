# Methods for objects of class validateFDboost

Methods for objects that are fitted to determine the optimal mstop and
the prediction error of a model fitted by FDboost.

## Usage

``` r
# S3 method for class 'validateFDboost'
mstop(object, riskopt = c("mean", "median"), ...)

# S3 method for class 'validateFDboost'
print(x, ...)

# S3 method for class 'validateFDboost'
plot(
  x,
  riskopt = c("mean", "median"),
  ylab = attr(x, "risk"),
  xlab = "Number of boosting iterations",
  ylim = range(x$oobrisk),
  which = 1,
  modObject = NULL,
  predictNA = FALSE,
  names.arg = NULL,
  ask = TRUE,
  ...
)

plotPredCoef(
  x,
  which = NULL,
  pers = TRUE,
  commonRange = TRUE,
  showNumbers = FALSE,
  showQuantiles = TRUE,
  ask = TRUE,
  terms = TRUE,
  probs = c(0.25, 0.5, 0.75),
  ylim = NULL,
  ...
)
```

## Arguments

- object:

  object of class `validateFDboost`

- riskopt:

  how the risk is minimized to obtain the optimal stopping iteration;
  defaults to the mean, can be changed to the median.

- ...:

  additional arguments passed to callies.

- x:

  an object of class `validateFDboost`.

- ylab:

  label for y-axis

- xlab:

  label for x-axis

- ylim:

  values for limits of y-axis

- which:

  In the case of `plotPredCoef()` the subset of base-learners to take
  into account for plotting. In the case of `plot.validateFDboost()` the
  diagnostic plots that are given (1: empirical risk per fold as a
  funciton of the boosting iterations, 2: empirical risk per fold, 3:
  MRD per fold, 4: observed and predicted values, 5: residuals; 2-5 for
  the model with the optimal number of boosting iterations).

- modObject:

  if the original model object of class `FDboost` is given predicted
  values of the whole model can be compared to the predictions of the
  cross-validated models

- predictNA:

  should missing values in the response be predicted? Defaults to
  `FALSE`.

- names.arg:

  names of the observed curves

- ask:

  defaults to `TRUE`, ask for next plot using `par(ask = ask)` ?

- pers:

  plot coefficient surfaces as persp-plots? Defaults to `TRUE`.

- commonRange, :

  plot predicted coefficients on a common range, defaults to `TRUE`.

- showNumbers:

  show number of curve in plot of predicted coefficients, defaults to
  `FALSE`

- showQuantiles:

  plot the 0.05 and the 0.95 Quantile of coefficients in 1-dim effects.

- terms:

  logical, defaults to `TRUE`; plot the added terms (default) or the
  coefficients?

- probs:

  vector of quantiles to be used in the plotting of 2-dimensional
  coefficients surfaces, defaults to `probs = c(0.25, 0.5, 0.75)`

## Value

No return value (plot method) or the object itself (print method)

## Details

The function `mstop.validateFDboost` extracts the optimal mstop by
minimizing the mean (or the median) risk. `plot.validateFDboost` plots
cross-validated risk, RMSE, MRD, measured and predicted values and
residuals as determined by `validateFDboost`. The function
`plotPredCoef` plots the coefficients that were estimated in the folds -
only possible if the argument getCoefCV is `TRUE` in the call to
`validateFDboost`.
