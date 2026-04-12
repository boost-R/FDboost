# Generic functions to asses attributes of functional data objects

Extract attributes of an object.

## Usage

``` r
getTime(object)

getId(object)

getX(object)

getArgvals(object)

getTimeLab(object)

getIdLab(object)

getXLab(object)

getArgvalsLab(object)
```

## Arguments

- object:

  an R-object, currently implemented for hmatrix and fmatrix

## Value

properties of a hmatrix or fmatrix

## Details

Extract the time variable `getTime`, the id`getId`, the functional
covariate `getX`, its argument values `getArgvals`. Or the names of the
different variables `getTimeLab`, `getIdLab`, `getXLab`,
`getArgvalsLab`.

## See also

[`hmatrix`](https://boost-R.github.io/FDboost/reference/hmatrix.md) for
the h.atrix class.
