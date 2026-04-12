# A S3 class for univariate functional data on a common grid

The hmatrix class represents data for a functional historical effect.
The class is basically a matrix containing the time and the id for the
observations of the functional response. The functional covariate is
contained as attribute.

## Usage

``` r
hmatrix(
  time,
  id,
  x,
  argvals = seq_len(ncol(x)),
  timeLab = "t",
  idLab = "wideIndex",
  xLab = "x",
  argvalsLab = "s"
)
```

## Arguments

- time:

  set of argument values of the response in long format, i.e. at which
  `t` the response curve is observed

- id:

  specify to which curve the point belongs to, id from 1, 2, ..., n.

- x:

  matrix of functional covariate, each trajectory is in one row

- argvals:

  set of argument values, i.e., the common gird at which the functional
  covariate is observed, by default `seq_len(ncol(x))`

- timeLab:

  name of the time axis, by default `t`

- idLab:

  name of the id variable, by default `wideIndex`

- xLab:

  name of the functional variable, by default NULL

- argvalsLab:

  name of the argument for the covariate by default `s`

## Value

An matrix object of type `"hmatrix"`

## Details

In the hmatrix class the id has to run from i=1, 2, ..., n including all
integers from 1 to n. The rows of the functional covariate x correspond
to those observations.

## See also

[`getTime.hmatrix`](https://boost-R.github.io/FDboost/reference/getTime.hmatrix.md)
to extract attributes, and ?"\[.hmatrix" for the extract method.

## Examples

``` r
## Example for a hmatrix object
t1 <- rep((1:5)/2, each = 3)
id1 <- rep(1:3, 5)
x1 <- matrix(1:15, ncol = 5) 
s1 <- (1:5)/2 
myhmatrix <- hmatrix(time = t1, id = id1, x = x1, argvals = s1, 
                     timeLab = "t1", argvalsLab = "s1", xLab = "test")

# extract with [ keeps attributes 
# select observations of subjects 2 and 3
myhmatrixSub <- myhmatrix[id1 %in% c(2, 3), ]  
str(myhmatrixSub)
#>  'hmatrix' num [1:10, 1:2] 0.5 0.5 1 1 1.5 1.5 2 2 2.5 2.5 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "time" "id"
#>  - attr(*, "x")= int [1:2, 1:5] 2 3 5 6 8 9 11 12 14 15
#>  - attr(*, "argvals")= num [1:5] 0.5 1 1.5 2 2.5
#>  - attr(*, "timeLab")= chr "t1"
#>  - attr(*, "idLab")= chr "wideIndex"
#>  - attr(*, "xLab")= chr "test"
#>  - attr(*, "argvalsLab")= chr "s1"
getX(myhmatrixSub)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    2    5    8   11   14
#> [2,]    3    6    9   12   15
getX(myhmatrix)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    4    7   10   13
#> [2,]    2    5    8   11   14
#> [3,]    3    6    9   12   15

# get time
myhmatrix[ , 1] # as column matrix as drop = FALSE
#>       time
#>  [1,]  0.5
#>  [2,]  0.5
#>  [3,]  0.5
#>  [4,]  1.0
#>  [5,]  1.0
#>  [6,]  1.0
#>  [7,]  1.5
#>  [8,]  1.5
#>  [9,]  1.5
#> [10,]  2.0
#> [11,]  2.0
#> [12,]  2.0
#> [13,]  2.5
#> [14,]  2.5
#> [15,]  2.5
getTime(myhmatrix) # as vector
#>  [1] 0.5 0.5 0.5 1.0 1.0 1.0 1.5 1.5 1.5 2.0 2.0 2.0 2.5 2.5 2.5

# get id
myhmatrix[ , 2] # as column matrix as drop = FALSE
#>       id
#>  [1,]  1
#>  [2,]  2
#>  [3,]  3
#>  [4,]  1
#>  [5,]  2
#>  [6,]  3
#>  [7,]  1
#>  [8,]  2
#>  [9,]  3
#> [10,]  1
#> [11,]  2
#> [12,]  3
#> [13,]  1
#> [14,]  2
#> [15,]  3
getId(myhmatrix) # as vector
#>  [1] 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3

# subset hmatrix on the basis of an index, which is defined on the curve level
reweightData(data = list(hmat = myhmatrix), vars = "hmat", index = c(1, 1, 2))
#> $hmat
#>       time id
#>  [1,]  0.5  1
#>  [2,]  0.5  2
#>  [3,]  0.5  3
#>  [4,]  1.0  1
#>  [5,]  1.0  2
#>  [6,]  1.0  3
#>  [7,]  1.5  1
#>  [8,]  1.5  2
#>  [9,]  1.5  3
#> [10,]  2.0  1
#> [11,]  2.0  2
#> [12,]  2.0  3
#> [13,]  2.5  1
#> [14,]  2.5  2
#> [15,]  2.5  3
#> attr(,"x")
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    4    7   10   13
#> [2,]    1    4    7   10   13
#> [3,]    2    5    8   11   14
#> attr(,"argvals")
#> [1] 0.5 1.0 1.5 2.0 2.5
#> attr(,"timeLab")
#> [1] "t1"
#> attr(,"idLab")
#> [1] "wideIndex"
#> attr(,"xLab")
#> [1] "test"
#> attr(,"argvalsLab")
#> [1] "s1"
#> attr(,"class")
#> [1] "hmatrix" "matrix"  "array"  
#> 
# this keeps only the unique x values in attr(,'x') but multiplies the corresponding
# ids and times in the time id matrix 
# for bhistx baselearner, there may be an additional id variable for the tensor product
newdat <- reweightData(data = list(hmat = myhmatrix, 
  repIDx = rep(seq_len(nrow(attr(myhmatrix,'x'))), length(attr(myhmatrix,"argvals")))), 
  vars = "hmat", index = c(1,1,2), idvars="repIDx")
length(newdat$repIDx) 
#> [1] 0

## use hmatrix within a data.frame
mydat <- data.frame(I(myhmatrix), z=rnorm(3)[id1])
str(mydat)
#> 'data.frame':    15 obs. of  2 variables:
#>  $ myhmatrix: 'AsIs' num [1:15, 1:2] 0.5 0.5 0.5 1 1 1 1.5 1.5 1.5 2 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "time" "id"
#>   ..- attr(*, "x")= int [1:3, 1:5] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..- attr(*, "argvals")= num [1:5] 0.5 1 1.5 2 2.5
#>   ..- attr(*, "timeLab")= chr "t1"
#>   ..- attr(*, "idLab")= chr "wideIndex"
#>   ..- attr(*, "xLab")= chr "test"
#>   ..- attr(*, "argvalsLab")= chr "s1"
#>  $ z        : num  0.545 -1.47 0.919 0.545 -1.47 ...
str(mydat[id1 %in% c(2, 3), ])
#> 'data.frame':    10 obs. of  2 variables:
#>  $ myhmatrix: 'AsIs' num [1:10, 1:2] 0.5 0.5 1 1 1.5 1.5 2 2 2.5 2.5 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "time" "id"
#>   ..- attr(*, "x")= int [1:2, 1:5] 2 3 5 6 8 9 11 12 14 15
#>   ..- attr(*, "argvals")= num [1:5] 0.5 1 1.5 2 2.5
#>   ..- attr(*, "timeLab")= chr "t1"
#>   ..- attr(*, "idLab")= chr "wideIndex"
#>   ..- attr(*, "xLab")= chr "test"
#>   ..- attr(*, "argvalsLab")= chr "s1"
#>  $ z        : num  -1.47 0.919 -1.47 0.919 -1.47 ...
str(myhmatrix[id1 %in% c(2, 3), ])
#>  'hmatrix' num [1:10, 1:2] 0.5 0.5 1 1 1.5 1.5 2 2 2.5 2.5 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : chr [1:2] "time" "id"
#>  - attr(*, "x")= int [1:2, 1:5] 2 3 5 6 8 9 11 12 14 15
#>  - attr(*, "argvals")= num [1:5] 0.5 1 1.5 2 2.5
#>  - attr(*, "timeLab")= chr "t1"
#>  - attr(*, "idLab")= chr "wideIndex"
#>  - attr(*, "xLab")= chr "test"
#>  - attr(*, "argvalsLab")= chr "s1"
```
