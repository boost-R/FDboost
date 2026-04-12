# Plot the fit or the coefficients of a boosted functional regression model

Takes a fitted `FDboost`-object produced by
[`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md)
and plots the fitted effects or the coefficient-functions/surfaces.

## Usage

``` r
# S3 method for class 'FDboost'
plot(
  x,
  raw = FALSE,
  rug = TRUE,
  which = NULL,
  includeOffset = TRUE,
  ask = TRUE,
  n1 = 40,
  n2 = 40,
  n3 = 20,
  n4 = 11,
  onlySelected = TRUE,
  pers = FALSE,
  commonRange = FALSE,
  ...
)

plotPredicted(
  x,
  subset = NULL,
  posLegend = "topleft",
  lwdObs = 1,
  lwdPred = 1,
  ...
)

plotResiduals(x, subset = NULL, posLegend = "topleft", ...)
```

## Arguments

- x:

  a fitted `FDboost`-object

- raw:

  logical defaults to `FALSE`. If `raw = FALSE` for each effect the
  estimated function/surface is calculated. If `raw = TRUE` the
  coefficients of the model are returned.

- rug:

  when `TRUE` (default) then the covariate to which the plot applies is
  displayed as a rug plot at the foot of each plot of a 1-d smooth, and
  the locations of the covariates are plotted as points on the contour
  plot representing a 2-d smooth.

- which:

  a subset of base-learners to take into account for plotting.

- includeOffset:

  logical, defaults to `TRUE`. Should the offset be included in the plot
  of the intercept (default) or should it be plotted separately.

- ask:

  logical, defaults to `TRUE`, if several effects are plotted the user
  has to hit Return to see next plot.

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

- onlySelected, :

  logical, defaults to `TRUE`. Only plot effects that where selected in
  at least one boosting iteration.

- pers:

  logical, defaults to `FALSE`, If `TRUE`, perspective plots
  ([`persp`](https://rdrr.io/r/graphics/persp.html)) for 2- and
  3-dimensional effects are drawn. If `FALSE`, image/contour-plots
  ([`image`](https://rdrr.io/r/graphics/image.html),
  [`contour`](https://rdrr.io/r/graphics/contour.html)) are drawn for 2-
  and 3-dimensional effects.

- commonRange:

  logical, defaults to `FALSE`, if `TRUE` the range over all effects is
  the same (does not affect perspecitve or image plots).

- ...:

  other arguments, passed to `funplot` (only used in plotPredicted)

- subset:

  subset of the observed response curves and their predictions that is
  plotted. Per default all observations are plotted.

- posLegend:

  location of the legend, if a legend is drawn automatically (only used
  in plotPredicted). The default is "topleft".

- lwdObs:

  lwd of observed curves (only used in plotPredicted)

- lwdPred:

  lwd of predicted curves (only used in plotPredicted)

## Value

no return value (plot method)

## See also

[`FDboost`](https://boost-R.github.io/FDboost/reference/FDboost.md) for
the model fit and
[`coef.FDboost`](https://boost-R.github.io/FDboost/reference/coef.FDboost.md)
for the calculation of the coefficient functions.
