# Extract attributes of hmatrix

Extract attributes of an object of class `hmatrix`.

## Usage

``` r
# S3 method for class 'hmatrix'
getTime(object)

# S3 method for class 'hmatrix'
getId(object)

# S3 method for class 'hmatrix'
getX(object)

# S3 method for class 'hmatrix'
getArgvals(object)

# S3 method for class 'hmatrix'
getTimeLab(object)

# S3 method for class 'hmatrix'
getIdLab(object)

# S3 method for class 'hmatrix'
getXLab(object)

# S3 method for class 'hmatrix'
getArgvalsLab(object)
```

## Arguments

- object:

  object of class hmatrix

## Value

properties of a hmatrix

## Details

Extract the time variable `getTime`, the id`getId`, the functional
covariate `getX`, its argument values `getArgvals`. Or the names of the
different variables `getTimeLab`, `getIdLab`, `getXLab`, `getArgvalsLab`
for an object of class `hmatrix`.
