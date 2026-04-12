# Extract or replace parts of a hmatrix-object

Operator acting on hmatrix preserving the attributes when rows are
extracted.

## Usage

``` r
# S3 method for class 'hmatrix'
x[i, j, ..., drop = FALSE]
```

## Arguments

- x:

  object from which to extract element(s) or in which to replace
  element(s).

- i, j:

  indices specifying elements to extract or replace. Indices are numeric
  vectors or empty (missing) or NULL. Numeric values are coerced to
  integer as by as.integer (and hence truncated towards zero).

- ...:

  not used

- drop:

  If `TRUE` the result is coerced to the lowest possible dimension (or
  just a matrix). This only works for extracting elements, not for the
  replacement, defaults to `FALSE`.

## Value

a `"hmatrix"` object

## Details

If used on columns or rows/columns a matrix is returned. If used on rows
only, i.e. x\[i,\] an object of class hmatrix is returned. The id is
changed so that it runs from 1, ..., nNew, where nNew is the number of
different id values in the new hmatrix-object. From the functional
covariate `x` rows are selected accordingly.

## See also

?"\["
