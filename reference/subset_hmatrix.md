# Subsets hmatrix according to an index

Subsets hmatrix according to an index

## Usage

``` r
subset_hmatrix(x, index, compress = TRUE)
```

## Arguments

- x:

  hmatix object that should be subsetted

- index:

  integer vector with (possibly duplicated) indices for each curve to
  select

- compress:

  logical, defaults to `TRUE`. Only used to force a meaningful behaviour
  of `applyFolds` with hmatrix objects when using nested resampling.

## Value

a `hmatrix` object

## Details

This methods is primary useful when subsetting repeatedly.

## Examples

``` r
t1 <- rep((1:5)/2, each = 3)
id1 <- rep(1:3, 5)
x1 <- matrix(1:15, ncol = 5) 
s1 <- (1:5)/2 
hmat <- hmatrix(time = t1, id = id1, x = x1, argvals = s1, timeLab = "t1", 
                argvalsLab = "s1", xLab = "test")

index1 <- c(1, 1, 3)
index2 <- c(2, 3, 3)
resMat <- subset_hmatrix(hmat, index = index1)
try(resMat2 <- subset_hmatrix(resMat, index = index2))
#> Error in attrTemp$x[index, , drop = FALSE] : subscript out of bounds
resMat <- subset_hmatrix(hmat, index = index1, compress = FALSE)
try(resMat2 <- subset_hmatrix(resMat, index = index2))
```
