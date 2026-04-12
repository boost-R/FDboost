# Function to Reweight Data

Function to Reweight Data

## Usage

``` r
reweightData(
  data,
  argvals,
  vars,
  longvars = NULL,
  weights,
  index,
  idvars = NULL,
  compress = FALSE
)
```

## Arguments

- data:

  a named list or data.frame.

- argvals:

  character (vector); name(s) for entries in data giving the index for
  observed grid points; must be supplied if `vars` is not supplied.

- vars:

  character (vector); name(s) for entries in data, which are subsetted
  according to weights or index. Must be supplied if `argvals` is not
  supplied.

- longvars:

  variables in long format, e.g., a response that is observed at curve
  specific grids.

- weights:

  vector of weights for observations. Must be supplied if `index` is not
  supplied.

- index:

  vector of indices for observations. Must be supplied if `weights` is
  not supplied.

- idvars:

  character (vector); index, which is needed to expand `vars` to be
  conform with the `hmatrix` structure when using `bhistx`-base-learners
  or to be conform with variables in long format specified in
  `longvars`.

- compress:

  logical; whether `hmatrix` objects are saved in compressed form or
  not. Default is `TRUE`. Should be set to `FALSE` when using
  `reweightData` for nested resampling.

## Value

A list with the reweighted or subsetted data.

## Details

`reweightData` indexes the rows of matrices and / or positions of
vectors by using either the `index` or the `weights`-argument. To
prevent the function from indexing the list entry / entries, which serve
as time index for observed grid points of each trajectory of functional
observations, the `argvals` argument (vector of character names for
these list entries) can be supplied. If `argvals` is not supplied,
`vars` must be supplied and it is assumed that `argvals` is equal to
`names(data)[!names(data) %in% vars]`.

When using `weights`, a weight vector of length N must be supplied,
where N is the number of observations. When using `index`, the vector
must contain the index of each row as many times as it shall be included
in the new data set.

## Author

David Ruegamer, Sarah Brockhaus

## Examples

``` r
## load data
data("viscosity", package = "FDboost")
interval <- "101"
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[ , 1:end])
viscosity$time <- viscosity$timeAll[1:end]

## what does data look like
str(viscosity)
#> List of 9
#>  $ visAll : 'AsIs' num [1:64, 1:132] 41.5 25.2 63.7 35.6 17.8 12.3 38.6 22 18.2 36 ...
#>  $ timeAll: num [1:132] 11 13 15 17 19 21 23 25 27 29 ...
#>  $ T_C    : Factor w/ 2 levels "low","high": 1 1 2 2 2 2 1 1 1 1 ...
#>  $ T_A    : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ T_B    : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 2 2 ...
#>  $ rspeed : Factor w/ 2 levels "low","high": 1 2 1 2 1 2 2 1 2 1 ...
#>  $ mflow  : Factor w/ 2 levels "low","high": 2 1 1 2 1 2 1 2 2 1 ...
#>  $ vis    : 'AsIs' num [1:64, 1:46] 3.73 3.23 4.15 3.57 2.88 ...
#>  $ time   : num [1:46] 11 13 15 17 19 21 23 25 27 29 ...

## do some reweighting
# correct weights
str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = c(0, 32, 32, rep(0, 61))))
#> List of 6
#>  $ vis   : 'AsIs' num [1:64, 1:46] 3.23 3.23 3.23 3.23 3.23 ...
#>  $ T_C   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ T_A   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ rspeed: Factor w/ 2 levels "low","high": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ mflow : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ time  : num [1:46] 11 13 15 17 19 21 23 25 27 29 ...

str(visNew <- reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = c(0, 32, 32, rep(0, 61))))
#> List of 6
#>  $ vis   : 'AsIs' num [1:64, 1:46] 3.23 3.23 3.23 3.23 3.23 ...
#>  $ T_C   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ T_A   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ rspeed: Factor w/ 2 levels "low","high": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ mflow : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ time  : num [1:46] 11 13 15 17 19 21 23 25 27 29 ...
# check the result
# visNew$vis[1:5, 1:5] ## image(visNew$vis)

# incorrect weights
str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = sample(1:64, replace = TRUE)), 1)
#> List of 6
#>  $ vis   : 'AsIs' num [1:2248, 1:46] 3.73 3.73 3.73 3.73 3.73 ...
#>  $ T_C   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ T_A   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ rspeed: Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ mflow : Factor w/ 2 levels "low","high": 2 2 2 2 2 2 2 2 2 2 ...
#>  $ time  : num [1:46] 11 13 15 17 19 21 23 25 27 29 ...

# supply meaningful index
str(visNew <- reweightData(viscosity, vars = c("vis", "T_C", "T_A", "rspeed", "mflow"), 
              argvals = "time", index = rep(1:32, each = 2)))
#> List of 6
#>  $ vis   : 'AsIs' num [1:64, 1:46] 3.73 3.73 3.23 3.23 4.15 ...
#>  $ T_C   : Factor w/ 2 levels "low","high": 1 1 1 1 2 2 2 2 2 2 ...
#>  $ T_A   : Factor w/ 2 levels "low","high": 1 1 1 1 1 1 1 1 1 1 ...
#>  $ rspeed: Factor w/ 2 levels "low","high": 1 1 2 2 1 1 2 2 1 1 ...
#>  $ mflow : Factor w/ 2 levels "low","high": 2 2 1 1 1 1 2 2 1 1 ...
#>  $ time  : num [1:46] 11 13 15 17 19 21 23 25 27 29 ...
# check the result
# visNew$vis[1:5, 1:5]

# errors
if(FALSE){
   reweightData(viscosity, argvals = "")
   reweightData(viscosity, argvals = "covThatDoesntExist", index = rep(1,64))
   }
   
```
