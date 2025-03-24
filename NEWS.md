# FDboost 1.1.0 (2022-07-12)

## Miscellaneous

- Anisotropic tensor-product operators `b1 %A0% b2` and `b1 %Xa0% b2` now also work when `lambda` is specified for `b1` and `df` is specified for `b2` (or vice versa).

## New features

- New function `clr()` to compute the centered-log-ratio transform and its inverse for density-on-scalar regression in Bayes spaces.
- New dataset `birthDistribution`.
- New vignette illustrating density-on-function regression on the `birthDistribution` data.
- Function `factorize()` added for tensor-product factorization of estimated effects or models.

# FDboost 0.3.4 (2020-08-31)

## Bug fixes

- Fix `predict()` for `bsignal()` with `newdata` and the functional covariate given as a numeric matrix, raised in [#17](https://github.com/boost-R/FDboost/issues/17).
- Deprecated argument `LINPACK` in `solve()` removed.

# FDboost 0.3.3 (2020-06-13)

## New features

- It is now possible to specify several time variables as well as factor time variables in the `timeformula`. This feature is needed for the manifoldboost package.

## Miscellaneous

- The function `stabsel.FDboost()` now uses `applyFolds()` instead of `validateFDboost()` to do cross-validation with recomputation of the smooth offset. This is only relevant for models with a functional response. This will change results if the model contains base-learners like `bbsc()` or `bolsc()`, as `applyFolds()` also recomputes the Z-matrix for those base-learners.

## Bug fixes

- Adapted functions `integrationWeights()` and `integrationWeightsLeft()` for unsorted time variables.
- Changed code in `predict.FDboost()` such that interaction effects of two functional covariates like `bsignal() %X% bsignal()` can be predicted with new data.
- Adapt FDboost to R 4.0.1 by explicitly using the first entry of `dots$aggregate` (i.e., `dots$aggregate[1] != "sum"`) in `predict.FDboost()` so that it also works with the default, where `aggregate` is a vector of length 3 and later only the first argument is used via `match.arg()`.

# FDboost 0.3.2 (2018-08-04)

## Bug fixes

- Deprecated argument `corrected` in `cvrisk()` removed.

# FDboost 0.3.1 (2018-05-10)

## Bug fixes

- `cvrisk()` has by default adequate folds for a noncyclic fitted FDboostLSS model, see [#14](https://github.com/boost-R/FDboost/issues/14).

## Miscellaneous

- Replaced `cBind()` (which is deprecated) with `cbind()`.

# FDboost 0.3.0 (2017-05-31)

## User-visible changes

- New function `bootstrapCI()` to compute bootstrapped coefficients.
- Added the dataset `emotion` containing EEG and EMG measures under different experimental conditions.
- With scalar response, `FDboost()` now works with the response as a vector (instead of a 1-row matrix); thus, `fitted()` and `predict()` return a vector.

## Bug fixes

- `update.FDboost()` now works with a scalar response.
- `FDboost()` works with family `Binomial(type = "glm")`, see [#1](https://github.com/boost-R/FDboost/issues/1).
- `applyFolds()` works for factor response, see [#7](https://github.com/boost-R/FDboost/issues/7).
- `cvLong()` and `cvMA()` return a matrix for only one resampling fold with `B = 1` (proposed by Almond Stoecker).

## Miscellaneous

- Adapt `FDboost` to `mboost` 2.8-0, which allows for `mstop = 0`.
- Restructure `FDboostLSS()` such that it calls `mboostLSS_fit()` from `gamboostLSS` 2.0-0.
- In `FDboost`, set `options("mboost_indexmin" = +Inf)` to disable internal use of ties in model fitting, as this breaks some methods for models with responses in long format and for models containing `bhistx()`, see [#10](https://github.com/boost-R/FDboost/issues/10).
- Deprecated `validateFDboost()`, use `applyFolds()` and `bootstrapCI()` instead.

# FDboost 0.2.0 (2016-05-26)

## User-visible changes

- Added function `applyFolds()` to compute the optimal stopping iteration.

## Bug fixes

- Allows for extrapolation in `predict()` with `bbsc()`.

# FDboost 0.1.2 (2016-04-22)

## Bug fixes

- Fixed a bug in `bolsc()`: correctly use the index in `bolsc()`/`bbsc()`. Previously, each observation was used only once for computing Z.

## User-visible changes

- Added function `%Xa0%` that computes a row-tensor product of two base-learners where the penalty in one direction is zero.
- Added function `reweightData()` that computes the data for Bootstrap or cross-validation folds.
- Added function `stabsel.FDboost()` that refits the smooth offset in each fold.
- Added argument `fun` to `validateFDboost()`.
- Added `update.FDboost()` that overwrites `update.mboost()`.

## Miscellaneous

- `FDboost()` works with `family = Binomial()`.

# FDboost 0.1.1 (2016-04-06)

## Bug fixes

- Fixed `oobpred` in `validateFDboost()` for irregular response and resampling at the curve level so that `plot.validateFDboost()` works for that case.
- Fixed scope of formula in `FDboost()`: now the formula given to `mboost()` within `FDboost()` uses the variables in the environment of the formula specified in `FDboost()`.

## Miscellaneous

- `plot.FDboost()` works for more effects, especially for effects like `bolsc() %X% bhistx()`.

# FDboost 0.1.0 (2016-03-10)

## User-visible changes

- New operator `%A0%` for Kronecker product of two base-learners with an anisotropic penalty for the special case where `lambda1` or `lambda2` is zero.
- The base-learner `bbsc()` can be used with `center = TRUE` (derived by Almond Stoecker).
- In `FDboostLSS()`, a list of one-sided formulas can be specified for `timeformula`.

## Bug fixes

- `FDboostLSS()` works with `families = GammaLSS()`.

## Miscellaneous

- Operator `%A%` uses weights in the model call. This only works correctly for weights on the level of `blg1` and `blg2` (same as weights on rows and columns of the response matrix).
- Calls to internal functions of `mboost` are done using `mboost_intern()`.
- `hyper_olsc()` is based on `hyper_ols()` from `mboost`.

# FDboost 0.0.17 (2016-02-25)

## User-visible changes

- Changed the operator `%Xc%` for the row tensor product of two scalar covariates. The design matrix of the interaction effects is constrained such that the interaction is centered around the intercept and around the two main effects of the scalar covariates (experimental!). Use, for example, `bols(x1) %Xc% bols(x2)`.

# FDboost 0.0.16 (2016-02-22)

## User-visible changes

- Changed the operator `%Xc%` for row tensor product where the sum-to-zero constraint is applied to the design matrix resulting from the row-tensor product (experimental!). Specifically, an intercept-column is first added, and then the sum-to-zero constraint is applied. Use, for example, `bolsc(x1) %Xc% bolsc(x2)`.
- The functional index `s` is now used as `argsvals` in the FPCA conducted within `bfpc()`.

# FDboost 0.0.15 (2016-02-12)

## User-visible changes

- New operator `%A%` that implies anisotropic penalties for differently specified `df` in the two base-learners.

## Bug fixes

- No penalty is applied in the direction of `ONEx` in a smooth intercept specified implicitly by `~1`, for example, `bols(ONEx, intercept=FALSE, df=1) %A% bbs(time)`.

## Miscellaneous

- Effects containing `%A%` or `%O%` are not expanded with the `timeformula`, allowing for different effects over time in the model.

# FDboost 0.0.14 (2016-02-11)

## User-visible changes

- Added the function `FDboostLSS()` to fit GAMLSS models with functional data using R-package `gamboostLSS`.
- New operator `%Xc%` for row tensor product where the sum-to-zero constraint is applied to the design matrix resulting from the row-tensor product (experimental!).
- Allowed `newdata` to be a list in `predict.FDboost()` when used with signal base-learners.
- Expanded `coef.FDboost()` so that it works for 3-dimensional tensor products of the form `bhistx() %X% bolsc() %X% bolsc()` (with David Ruegamer).
- Added a new possibility for scalar-on-function regression: if `timeformula=NULL`, no Kronecker product with `1` is used, which changes the penalty (otherwise, the direction of `1` would also be penalized).

## Miscellaneous

- New dependency on R-package `gamboostLSS`.
- Removed dependency on R-package `MASS`.
- Used the argument `prediction` in the internal computation of the base-learners (work in progress).
- Throw an error if `timeLab` of the `hmatrix`-object in `bhistx()` is not equal to the time variable in `timeformula`.

# FDboost 0.0.13 (2015-11-17)

## User-visible changes

- In function `FDboost()`, the offset is supplied differently. For a scalar offset, use `offset = "scalar"`. The default remains `offset = NULL`.
- `predict.FDboost()` has a new argument `toFDboost` (logical).
- `fitted.FDboost()` has argument `toFDboost` explicitly (not only via `...`).
- New base-learner `bhistx()`, especially suited for effects used with `%X%`, e.g., `bhistx() %X% bolsc()`.
- `coef.FDboost()` and `plot.FDboost()` now handle effects like `bhistx() %X% bolsc()`.
- For `predict.FDboost()` with effects `bhistx()` and newdata, the latest `mboostPatch` is necessary.

## Bug fixes

- The check for the necessity of a smooth offset works for missing values in a regular response (spotted by Tore Erdmann).

# FDboost 0.0.12 (2015-09-15)

- Internal experimental version.

# FDboost 0.0.11 (2015-06-01)

## User-visible changes

- `integrationWeights()` now gives equal weights for regular grids.
- New base-learner `bfpc()` for a functional covariate where both the functional covariate and the coefficient are expanded using fPCA (experimental feature!). Only works for regularly observed functional covariate.

## Bug fixes

- `coef.FDboost()` only works for `bhist()` if the time variable is the same in the timeformula and in `bhist()`.
- `predict.FDboost()` now checks that only `type = "link"` can be predicted for newdata.

# FDboost 0.0.10 (2015-04-16)

## User-visible changes

- Changed the default difference penalties to first-order difference (`differences = 1`), improving identifiability.
- New method `cvrisk.FDboost()` that uses (by default) sampling on the levels of curves, which is important for functional responses.
- Reorganized documentation of `cvrisk()` and `validateFDboost()`.
- In `bhist()`, an effect can be standardized.

## Miscellaneous

- Added a `CITATION` file.
- Uses `mboost 2.4-2`, which exports all important functions.

## Bug fixes

- `main` argument is always passed in `plot.FDboost()`.
- `bhist()` and `bconcurrent()` now work for equal `time` and `s`.
- `predict.FDboost()` works with tensor-product base-learners like `bl1 %X% bl2`.
