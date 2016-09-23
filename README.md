FDboost
======

<!-- [![Build Status (Linux)](https://travis-ci.org/boost-R/mboost.svg?branch=master)](https://travis-ci.org/boost-R/mboost)
[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/5mkvicgin1j6pfc6/branch/master?svg=true)](https://ci.appveyor.com/project/hofnerb/mboost-h73a1/branch/master) -->
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/FDboost)](http://cran.r-project.org/package=FDboost)
<!--[![Coverage Status](https://coveralls.io/repos/github/boost-R/mboost/badge.svg?branch=master)](https://coveralls.io/github/boost-R/mboost?branch=master) -->
[![](http://cranlogs.r-pkg.org/badges/FDboost)](http://cran.rstudio.com/web/packages/FDboost/index.html)

`FDboost` Boosting Functional Regression Models.

## Using FDboost

For installation instructions see below.

Instructions on how to use `FDboost` can be found in various places:
- Have a look at the vignettes:
  - [function-on-function regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_canada.pdf)
  - [scalar-on-function regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_fuel.pdf)
  - [function-on-scalar regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_viscosity.pdf)

## Issues & Feature Requests

For issues, bugs, feature requests etc. please use the [GitHub Issues](https://github.com/fdboost/FDboost/issues).

## Installation Instructions

- Current version (from CRAN):
  ```r
  install.packages("FDboost")
  ```

- Latest **patch version** (patched version of CRAN package; under development) from GitHub:
  ```r
  library("devtools")
  install_github("fdboost/FDboost", build_vignettes = TRUE)
  library("FDboost")
  ```

<!-- - Latest **development version** (version with new features; under development) from GitHub:
  ```r
  library("devtools")
  install_github("boost-R/mboost", ref = "devel")
  library("mboost")
  ```
-->

  To be able to use the `install_github()` command, one needs to install `devtools` first:
  ```r
  install.packages("devtools")
  ```

