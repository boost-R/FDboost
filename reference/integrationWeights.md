# Functions to compute integration weights

Computes trapezoidal integration weights (Riemann sums) for a functional
variable `X1` that has evaluation points `xind`.

## Usage

``` r
integrationWeights(X1, xind, id = NULL)

integrationWeightsLeft(X1, xind, leftWeight = c("first", "mean", "zero"))
```

## Arguments

- X1:

  for functional data that is observed on one common grid, a matrix
  containing the observations of the functional variable. For a
  functional variable that is observed on curve specific grids, a long
  vector.

- xind:

  evaluation points (index) of functional variable

- id:

  defaults to `NULL`. Only necessary for response in long format. In
  this case `id` specifies which curves belong together.

- leftWeight:

  one of `c("mean", "first", "zero")`. With left Riemann sums different
  assumptions for the weight of the first observation are possible. The
  default is to use the mean over all integration weights, `"mean"`.
  Alternatively one can use the first integration weight, `"first"`, or
  use the distance to zero, `"zero"`.

## Value

Matrix with integration

## Details

The function `integrationWeights()` computes trapezoidal integration
weights, that are symmetric. Per default those weights are used in the
[`bsignal`](https://boost-R.github.io/FDboost/reference/bsignal.md)-base-learner.
In the special case of evaluation points (`xind`) with equal distances,
all integration weights are equal.

The function `integrationWeightsLeft()` computes weights, that take into
account only the distance to the prior observation point. Thus one has
to decide what to do with the first observation. The left weights are
adequate for historical effects like in
[`bhist`](https://boost-R.github.io/FDboost/reference/bsignal.md).

## See also

[`bsignal`](https://boost-R.github.io/FDboost/reference/bsignal.md) and
[`bhist`](https://boost-R.github.io/FDboost/reference/bsignal.md) for
the base-learners.

## Examples

``` r
## Example for trapezoidal integration weights
xind0 <- seq(0,1,l = 5)
xind <- c(0, 0.1, 0.3, 0.7, 1)
X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)

# Regualar observation points
integrationWeights(X1, xind0)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.25 0.25 0.25 0.25 0.25
#> [2,] 0.25 0.25 0.25 0.25 0.25
# Irregular observation points
integrationWeights(X1, xind)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.05 0.15  0.3 0.35 0.15
#> [2,] 0.05 0.15  0.3 0.35 0.15

# with missing value
X1[1,2] <- NA
integrationWeights(X1, xind0)
#>      [,1] [,2]  [,3] [,4]  [,5]
#> [1,] 0.25 0.00 0.375 0.25 0.125
#> [2,] 0.25 0.25 0.250 0.25 0.250
integrationWeights(X1, xind)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.15 0.00 0.35 0.35 0.15
#> [2,] 0.05 0.15 0.30 0.35 0.15

## Example for left integration weights
xind0 <- seq(0,1,l = 5)
xind <- c(0, 0.1, 0.3, 0.7, 1)
X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)

# Regular observation points
integrationWeightsLeft(X1, xind0, leftWeight = "mean") 
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.25 0.25 0.25 0.25 0.25
#> [2,] 0.25 0.25 0.25 0.25 0.25
integrationWeightsLeft(X1, xind0, leftWeight = "first") 
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.25 0.25 0.25 0.25 0.25
#> [2,] 0.25 0.25 0.25 0.25 0.25
integrationWeightsLeft(X1, xind0, leftWeight = "zero")
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0 0.25 0.25 0.25 0.25
#> [2,]    0 0.25 0.25 0.25 0.25

# Irregular observation points
integrationWeightsLeft(X1, xind, leftWeight = "mean") 
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,] 0.25  0.1  0.2  0.4  0.3
#> [2,] 0.25  0.1  0.2  0.4  0.3
integrationWeightsLeft(X1, xind, leftWeight = "first") 
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  0.1  0.1  0.2  0.4  0.3
#> [2,]  0.1  0.1  0.2  0.4  0.3
integrationWeightsLeft(X1, xind, leftWeight = "zero")
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    0  0.1  0.2  0.4  0.3
#> [2,]    0  0.1  0.2  0.4  0.3

# obervation points that do not start with 0
xind2 <- xind + 0.5
integrationWeightsLeft(X1, xind2, leftWeight = "zero")
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  0.5  0.1  0.2  0.4  0.3
#> [2,]  0.5  0.1  0.2  0.4  0.3
 
```
