# Plot functional data with linear interpolation of missing values

Plot functional data with linear interpolation of missing values

## Usage

``` r
funplot(x, y, id = NULL, rug = TRUE, ...)
```

## Arguments

- x:

  optional, time-vector for plotting

- y:

  matrix of functional data with functions in rows and measured times in
  columns; or vector or functional observations, in this case id has to
  be specified

- id:

  defaults to NULL for y matrix, is id-variables for y in long format

- rug:

  logical. Should rugs be plotted? Defaults to TRUE.

- ...:

  further arguments passed to
  [`matplot`](https://rdrr.io/r/graphics/matplot.html).

## Value

see [`matplot`](https://rdrr.io/r/graphics/matplot.html)

## Details

All observations are marked by a small cross (`pch=3`). Missing values
are imputed by linear interpolation. Parts that are interpolated are
plotted by dotted lines, parts with non-missing values as solid lines.

## Examples

``` r
# \donttest{
### examples for regular data in wide format
data(viscosity)
with(viscosity, funplot(timeAll, visAll, pch=20))

if(require(fda)){
  with(fda::growth, funplot(age, t(hgtm)))
}

# }
```
