FDboost
======

[![Build Status (Linux)](https://travis-ci.org/boost-R/FDboost.svg?branch=master)](https://travis-ci.org/boost-R/FDboost)
<!--[![Build status (Windows)](https://ci.appveyor.com/api/projects/status/5mkvicgin1j6pfc6/branch/master?svg=true)](https://ci.appveyor.com/project/hofnerb/mboost-h73a1/branch/master)  -->
[![CRAN Status Badge](https://www.r-pkg.org/badges/version/FDboost)](https://cran.r-project.org/package=FDboost)
<!--[![Coverage Status](https://coveralls.io/repos/github/boost-R/mboost/badge.svg?branch=master)](https://coveralls.io/github/boost-R/mboost?branch=master) -->
[![](https://cranlogs.r-pkg.org/badges/FDboost)](https://cran.rstudio.com/web/packages/FDboost/index.html)

`FDboost` Boosting Functional Regression Models.

The package FDboost fits regression models for functional data, i.e., 
scalar-on-function, function-on-scalar and function-on-function regression models, 
by a component-wise gradient boosting algorithm. 

## Using FDboost

For installation instructions see below.

Instructions on how to use `FDboost` can be found in various places:
- Read the tutorial paper [arXiv:1705.10662](https://arxiv.org/abs/1705.10662)
- Have a look at the manual, which also contains example code
- Check the vignettes: 
  - [function-on-function regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_canada.pdf)
  - [scalar-on-function regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_fuel.pdf)
  - [function-on-scalar regression](https://cran.r-project.org/web/packages/FDboost/vignettes/FLAM_viscosity.pdf)

## Issues & Feature Requests

For issues, bugs, feature requests etc. please use the [GitHub Issues](https://github.com/boost-R/FDboost/issues).

## Installation Instructions

- Current version (from CRAN):
  ```r
  install.packages("FDboost")
  ```

- Latest **patch version** (patched version of CRAN package; under development) from GitHub:
  ```r
  library("devtools")
  install_github("boost-R/FDboost") 
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

