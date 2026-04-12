# Changelog

## FDboost 1.1-4 (2026-03-24)

CRAN release: 2026-03-24

### Bug fixes

- Suppressed deprecation warnings from
  [`refund::pffrSim()`](https://rdrr.io/pkg/refund/man/pffrSim.html) in
  examples to keep checks clean with recent `refund` versions.
- Stabilized `tests/general_tests.R` by forcing serial evaluation
  (`mc.cores = 1`) to avoid parallel
  [`applyFolds()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  crashes on some check platforms.

## FDboost 1.1.0 (2022-07-12)

### Miscellaneous

- Anisotropic tensor-product operators `b1 %A0% b2` and `b1 %Xa0% b2`
  now also work when `lambda` is specified for `b1` and `df` is
  specified for `b2` (or vice versa).

### New features

- New function
  [`clr()`](https://boost-R.github.io/FDboost/reference/clr.md) to
  compute the centered-log-ratio transform and its inverse for
  density-on-scalar regression in Bayes spaces.
- New dataset `birthDistribution`.
- New vignette illustrating density-on-function regression on the
  `birthDistribution` data.
- Function
  [`factorize()`](https://boost-R.github.io/FDboost/reference/factorize.md)
  added for tensor-product factorization of estimated effects or models.

## FDboost 0.3.4 (2020-08-31)

### Bug fixes

- Fix [`predict()`](https://rdrr.io/r/stats/predict.html) for
  [`bsignal()`](https://boost-R.github.io/FDboost/reference/bsignal.md)
  with `newdata` and the functional covariate given as a numeric matrix,
  raised in [\#17](https://github.com/boost-R/FDboost/issues/17).
- Deprecated argument `LINPACK` in
  [`solve()`](https://rdrr.io/r/base/solve.html) removed.

## FDboost 0.3.3 (2020-06-13)

### New features

- It is now possible to specify several time variables as well as factor
  time variables in the `timeformula`. This feature is needed for the
  manifoldboost package.

### Miscellaneous

- The function
  [`stabsel.FDboost()`](https://boost-R.github.io/FDboost/reference/stabsel.FDboost.md)
  now uses
  [`applyFolds()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  instead of
  [`validateFDboost()`](https://boost-R.github.io/FDboost/reference/validateFDboost.md)
  to do cross-validation with recomputation of the smooth offset. This
  is only relevant for models with a functional response. This will
  change results if the model contains base-learners like
  [`bbsc()`](https://boost-R.github.io/FDboost/reference/bbsc.md) or
  [`bolsc()`](https://boost-R.github.io/FDboost/reference/bbsc.md), as
  [`applyFolds()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  also recomputes the Z-matrix for those base-learners.

### Bug fixes

- Adapted functions
  [`integrationWeights()`](https://boost-R.github.io/FDboost/reference/integrationWeights.md)
  and
  [`integrationWeightsLeft()`](https://boost-R.github.io/FDboost/reference/integrationWeights.md)
  for unsorted time variables.
- Changed code in
  [`predict.FDboost()`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)
  such that interaction effects of two functional covariates like
  `bsignal() %X% bsignal()` can be predicted with new data.
- Adapt FDboost to R 4.0.1 by explicitly using the first entry of
  `dots$aggregate` (i.e., `dots$aggregate[1] != "sum"`) in
  [`predict.FDboost()`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)
  so that it also works with the default, where `aggregate` is a vector
  of length 3 and later only the first argument is used via
  [`match.arg()`](https://rdrr.io/r/base/match.arg.html).

## FDboost 0.3.2 (2018-08-04)

### Bug fixes

- Deprecated argument `corrected` in
  [`cvrisk()`](https://rdrr.io/pkg/mboost/man/cvrisk.html) removed.

## FDboost 0.3.1 (2018-05-10)

### Bug fixes

- [`cvrisk()`](https://rdrr.io/pkg/mboost/man/cvrisk.html) has by
  default adequate folds for a noncyclic fitted FDboostLSS model, see
  [\#14](https://github.com/boost-R/FDboost/issues/14).

### Miscellaneous

- Replaced `cBind()` (which is deprecated) with
  [`cbind()`](https://rdrr.io/r/base/cbind.html).

## FDboost 0.3.0 (2017-05-31)

### User-visible changes

- New function
  [`bootstrapCI()`](https://boost-R.github.io/FDboost/reference/bootstrapCI.md)
  to compute bootstrapped coefficients.
- Added the dataset `emotion` containing EEG and EMG measures under
  different experimental conditions.
- With scalar response,
  [`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md)
  now works with the response as a vector (instead of a 1-row matrix);
  thus, [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html) return a vector.

### Bug fixes

- [`update.FDboost()`](https://boost-R.github.io/FDboost/reference/update.FDboost.md)
  now works with a scalar response.
- [`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md)
  works with family `Binomial(type = "glm")`, see
  [\#1](https://github.com/boost-R/FDboost/issues/1).
- [`applyFolds()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  works for factor response, see
  [\#7](https://github.com/boost-R/FDboost/issues/7).
- [`cvLong()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  and `cvMA()` return a matrix for only one resampling fold with `B = 1`
  (proposed by Almond Stoecker).

### Miscellaneous

- Adapt `FDboost` to `mboost` 2.8-0, which allows for `mstop = 0`.
- Restructure
  [`FDboostLSS()`](https://boost-R.github.io/FDboost/reference/FDboostLSS.md)
  such that it calls
  [`mboostLSS_fit()`](https://rdrr.io/pkg/gamboostLSS/man/mboostLSS.html)
  from `gamboostLSS` 2.0-0.
- In `FDboost`, set `options("mboost_indexmin" = +Inf)` to disable
  internal use of ties in model fitting, as this breaks some methods for
  models with responses in long format and for models containing
  [`bhistx()`](https://boost-R.github.io/FDboost/reference/bhistx.md),
  see [\#10](https://github.com/boost-R/FDboost/issues/10).
- Deprecated
  [`validateFDboost()`](https://boost-R.github.io/FDboost/reference/validateFDboost.md),
  use
  [`applyFolds()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  and
  [`bootstrapCI()`](https://boost-R.github.io/FDboost/reference/bootstrapCI.md)
  instead.

## FDboost 0.2.0 (2016-05-26)

### User-visible changes

- Added function
  [`applyFolds()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  to compute the optimal stopping iteration.

### Bug fixes

- Allows for extrapolation in
  [`predict()`](https://rdrr.io/r/stats/predict.html) with
  [`bbsc()`](https://boost-R.github.io/FDboost/reference/bbsc.md).

## FDboost 0.1.2 (2016-04-22)

### Bug fixes

- Fixed a bug in
  [`bolsc()`](https://boost-R.github.io/FDboost/reference/bbsc.md):
  correctly use the index in
  [`bolsc()`](https://boost-R.github.io/FDboost/reference/bbsc.md)/[`bbsc()`](https://boost-R.github.io/FDboost/reference/bbsc.md).
  Previously, each observation was used only once for computing Z.

### User-visible changes

- Added function `%Xa0%` that computes a row-tensor product of two
  base-learners where the penalty in one direction is zero.
- Added function
  [`reweightData()`](https://boost-R.github.io/FDboost/reference/reweightData.md)
  that computes the data for Bootstrap or cross-validation folds.
- Added function
  [`stabsel.FDboost()`](https://boost-R.github.io/FDboost/reference/stabsel.FDboost.md)
  that refits the smooth offset in each fold.
- Added argument `fun` to
  [`validateFDboost()`](https://boost-R.github.io/FDboost/reference/validateFDboost.md).
- Added
  [`update.FDboost()`](https://boost-R.github.io/FDboost/reference/update.FDboost.md)
  that overwrites `update.mboost()`.

### Miscellaneous

- [`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md)
  works with `family = Binomial()`.

## FDboost 0.1.1 (2016-04-06)

### Bug fixes

- Fixed `oobpred` in
  [`validateFDboost()`](https://boost-R.github.io/FDboost/reference/validateFDboost.md)
  for irregular response and resampling at the curve level so that
  [`plot.validateFDboost()`](https://boost-R.github.io/FDboost/reference/plot.validateFDboost.md)
  works for that case.
- Fixed scope of formula in
  [`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md):
  now the formula given to
  [`mboost()`](https://rdrr.io/pkg/mboost/man/gamboost.html) within
  [`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md)
  uses the variables in the environment of the formula specified in
  [`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md).

### Miscellaneous

- [`plot.FDboost()`](https://boost-R.github.io/FDboost/reference/plot.FDboost.md)
  works for more effects, especially for effects like
  `bolsc() %X% bhistx()`.

## FDboost 0.1.0 (2016-03-10)

### User-visible changes

- New operator `%A0%` for Kronecker product of two base-learners with an
  anisotropic penalty for the special case where `lambda1` or `lambda2`
  is zero.
- The base-learner
  [`bbsc()`](https://boost-R.github.io/FDboost/reference/bbsc.md) can be
  used with `center = TRUE` (derived by Almond Stoecker).
- In
  [`FDboostLSS()`](https://boost-R.github.io/FDboost/reference/FDboostLSS.md),
  a list of one-sided formulas can be specified for `timeformula`.

### Bug fixes

- [`FDboostLSS()`](https://boost-R.github.io/FDboost/reference/FDboostLSS.md)
  works with `families = GammaLSS()`.

### Miscellaneous

- Operator `%A%` uses weights in the model call. This only works
  correctly for weights on the level of `blg1` and `blg2` (same as
  weights on rows and columns of the response matrix).
- Calls to internal functions of `mboost` are done using
  [`mboost_intern()`](https://rdrr.io/pkg/mboost/man/mboost_intern.html).
- `hyper_olsc()` is based on `hyper_ols()` from `mboost`.

## FDboost 0.0.17 (2016-02-25)

### User-visible changes

- Changed the operator `%Xc%` for the row tensor product of two scalar
  covariates. The design matrix of the interaction effects is
  constrained such that the interaction is centered around the intercept
  and around the two main effects of the scalar covariates
  (experimental!). Use, for example, `bols(x1) %Xc% bols(x2)`.

## FDboost 0.0.16 (2016-02-22)

### User-visible changes

- Changed the operator `%Xc%` for row tensor product where the
  sum-to-zero constraint is applied to the design matrix resulting from
  the row-tensor product (experimental!). Specifically, an
  intercept-column is first added, and then the sum-to-zero constraint
  is applied. Use, for example, `bolsc(x1) %Xc% bolsc(x2)`.
- The functional index `s` is now used as `argsvals` in the FPCA
  conducted within
  [`bfpc()`](https://boost-R.github.io/FDboost/reference/bsignal.md).

## FDboost 0.0.15 (2016-02-12)

### User-visible changes

- New operator `%A%` that implies anisotropic penalties for differently
  specified `df` in the two base-learners.

### Bug fixes

- No penalty is applied in the direction of `ONEx` in a smooth intercept
  specified implicitly by `~1`, for example,
  `bols(ONEx, intercept=FALSE, df=1) %A% bbs(time)`.

### Miscellaneous

- Effects containing `%A%` or `%O%` are not expanded with the
  `timeformula`, allowing for different effects over time in the model.

## FDboost 0.0.14 (2016-02-11)

### User-visible changes

- Added the function
  [`FDboostLSS()`](https://boost-R.github.io/FDboost/reference/FDboostLSS.md)
  to fit GAMLSS models with functional data using R-package
  `gamboostLSS`.
- New operator `%Xc%` for row tensor product where the sum-to-zero
  constraint is applied to the design matrix resulting from the
  row-tensor product (experimental!).
- Allowed `newdata` to be a list in
  [`predict.FDboost()`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)
  when used with signal base-learners.
- Expanded
  [`coef.FDboost()`](https://boost-R.github.io/FDboost/reference/coef.FDboost.md)
  so that it works for 3-dimensional tensor products of the form
  `bhistx() %X% bolsc() %X% bolsc()` (with David Ruegamer).
- Added a new possibility for scalar-on-function regression: if
  `timeformula=NULL`, no Kronecker product with `1` is used, which
  changes the penalty (otherwise, the direction of `1` would also be
  penalized).

### Miscellaneous

- New dependency on R-package `gamboostLSS`.
- Removed dependency on R-package `MASS`.
- Used the argument `prediction` in the internal computation of the
  base-learners (work in progress).
- Throw an error if `timeLab` of the `hmatrix`-object in
  [`bhistx()`](https://boost-R.github.io/FDboost/reference/bhistx.md) is
  not equal to the time variable in `timeformula`.

## FDboost 0.0.13 (2015-11-17)

### User-visible changes

- In function
  [`FDboost()`](https://boost-R.github.io/FDboost/reference/FDboost.md),
  the offset is supplied differently. For a scalar offset, use
  `offset = "scalar"`. The default remains `offset = NULL`.
- [`predict.FDboost()`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)
  has a new argument `toFDboost` (logical).
- [`fitted.FDboost()`](https://boost-R.github.io/FDboost/reference/fitted.FDboost.md)
  has argument `toFDboost` explicitly (not only via `...`).
- New base-learner
  [`bhistx()`](https://boost-R.github.io/FDboost/reference/bhistx.md),
  especially suited for effects used with `%X%`, e.g.,
  `bhistx() %X% bolsc()`.
- [`coef.FDboost()`](https://boost-R.github.io/FDboost/reference/coef.FDboost.md)
  and
  [`plot.FDboost()`](https://boost-R.github.io/FDboost/reference/plot.FDboost.md)
  now handle effects like `bhistx() %X% bolsc()`.
- For
  [`predict.FDboost()`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)
  with effects
  [`bhistx()`](https://boost-R.github.io/FDboost/reference/bhistx.md)
  and newdata, the latest `mboostPatch` is necessary.

### Bug fixes

- The check for the necessity of a smooth offset works for missing
  values in a regular response (spotted by Tore Erdmann).

## FDboost 0.0.12 (2015-09-15)

- Internal experimental version.

## FDboost 0.0.11 (2015-06-01)

### User-visible changes

- [`integrationWeights()`](https://boost-R.github.io/FDboost/reference/integrationWeights.md)
  now gives equal weights for regular grids.
- New base-learner
  [`bfpc()`](https://boost-R.github.io/FDboost/reference/bsignal.md) for
  a functional covariate where both the functional covariate and the
  coefficient are expanded using fPCA (experimental feature!). Only
  works for regularly observed functional covariate.

### Bug fixes

- [`coef.FDboost()`](https://boost-R.github.io/FDboost/reference/coef.FDboost.md)
  only works for
  [`bhist()`](https://boost-R.github.io/FDboost/reference/bsignal.md) if
  the time variable is the same in the timeformula and in
  [`bhist()`](https://boost-R.github.io/FDboost/reference/bsignal.md).
- [`predict.FDboost()`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)
  now checks that only `type = "link"` can be predicted for newdata.

## FDboost 0.0.10 (2015-04-16)

### User-visible changes

- Changed the default difference penalties to first-order difference
  (`differences = 1`), improving identifiability.
- New method
  [`cvrisk.FDboost()`](https://boost-R.github.io/FDboost/reference/applyFolds.md)
  that uses (by default) sampling on the levels of curves, which is
  important for functional responses.
- Reorganized documentation of
  [`cvrisk()`](https://rdrr.io/pkg/mboost/man/cvrisk.html) and
  [`validateFDboost()`](https://boost-R.github.io/FDboost/reference/validateFDboost.md).
- In
  [`bhist()`](https://boost-R.github.io/FDboost/reference/bsignal.md),
  an effect can be standardized.

### Miscellaneous

- Added a `CITATION` file.
- Uses `mboost 2.4-2`, which exports all important functions.

### Bug fixes

- `main` argument is always passed in
  [`plot.FDboost()`](https://boost-R.github.io/FDboost/reference/plot.FDboost.md).
- [`bhist()`](https://boost-R.github.io/FDboost/reference/bsignal.md)
  and
  [`bconcurrent()`](https://boost-R.github.io/FDboost/reference/bsignal.md)
  now work for equal `time` and `s`.
- [`predict.FDboost()`](https://boost-R.github.io/FDboost/reference/predict.FDboost.md)
  works with tensor-product base-learners like `bl1 %X% bl2`.
