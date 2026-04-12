# Extract information of a base-learner

Takes a base-learner and extracts information.

## Usage

``` r
# S3 method for class 'blg'
extract(
  object,
  what = c("design", "penalty", "index"),
  asmatrix = FALSE,
  expand = FALSE,
  ...
)
```

## Arguments

- object:

  a base-learner

- what:

  a character specifying the quantities to extract. This can be a subset
  of "design" (default; design matrix), "penalty" (penalty matrix) and
  "index" (index of ties used to expand the design matrix)

- asmatrix:

  a logical indicating whether the the returned matrix should be coerced
  to a matrix (default) or if the returned object stays as it is (i.e.,
  potentially a sparse matrix). This option is only applicable if
  `extract` returns matrices, i.e., `what = "design"` or
  `what = "penalty"`.

- expand:

  a logical indicating whether the design matrix should be expanded
  (default: `FALSE`). This is useful if ties were taken into account
  either manually (via argument `index` in a base-learner) or
  automatically for data sets with many observations. `expand = TRUE` is
  equivalent to `extract(B)[extract(B, what = "index"),]` for a
  base-learner `B`.

- ...:

  currently not used

## See also

`extract` for the `extract` function of the package `mboost`.
