# Function to control estimation of smooth offset

Function to control estimation of smooth offset

## Usage

``` r
o_control(k_min = 20, rule = 2, silent = TRUE, cyclic = FALSE, knots = NULL)
```

## Arguments

- k_min:

  maximal number of k in s()

- rule:

  which rule to use in approx() of the response before calculating the
  global mean, rule=1 means no extrapolation, rule=2 means to
  extrapolate the closest non-missing value, see
  [`approx`](https://rdrr.io/r/stats/approxfun.html)

- silent:

  print error messages of model fit?

- cyclic:

  defaults to FALSE, if TRUE cyclic splines are used

- knots:

  arguments knots passed to
  [`gam`](https://rdrr.io/pkg/mgcv/man/gam.html)

## Value

a list with controls
