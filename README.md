# FDboost

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/FDboost)](https://CRAN.R-project.org/package=FDboost)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/FDboost)](https://www.r-pkg.org/pkg/FDboost)

<!-- badges: end -->

`FDboost` Boosting Functional Regression Models.

The package FDboost fits regression models for functional data, i.e.,
scalar-on-function, function-on-scalar, and function-on-function regression models,
by a component-wise gradient boosting algorithm.
Furthermore, it can be used to fit density-on-scalar regression models.

## Using FDboost

For installation instructions see below.

Instructions on how to use `FDboost` can be found in various places:

- Read the tutorial paper [doi:10.18637/jss.v094.i10](doi:10.18637/jss.v094.i10)
- Have a look at the manual, which also contains example code
- Check the vignettes:
  - [function-on-function regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_canada.pdf)
  - [scalar-on-function regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_fuel.pdf)
  - [function-on-scalar regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_viscosity.pdf)
  - [density-on-scalar regression](https://github.com/Eva2703/FDboost/blob/BayesSpace/vignettes/density-on-scalar_birth.pdf)

## Issues & Feature Requests

For issues, bugs, feature requests etc. please use the [GitHub Issues](https://github.com/boost-R/FDboost/issues).

## Installation

Install the last release from [CRAN](https://cran.r-project.org):

```r
install.packages("FDboost")
```

Install the development version from [GitHub](https://github.com/):

```r
# install.packages("pak")
pak::pak("boost-R/FDboost")
```
