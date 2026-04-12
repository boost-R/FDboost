# Function to update FDboost objects

Function to update FDboost objects

## Usage

``` r
# S3 method for class 'FDboost'
update(
  object,
  weights = NULL,
  oobweights = NULL,
  risk = NULL,
  trace = NULL,
  ...,
  evaluate = TRUE
)
```

## Arguments

- object:

  fitted FDboost-object

- weights, oobweights, risk, trace:

  see
  [`?FDboost`](https://boost-R.github.io/FDboost/reference/FDboost.md)

- ...:

  Additional arguments to the call, or arguments with changed values.

- evaluate:

  If true evaluate the new call else return the call.

## Value

Returns the call of (`evaluate = FALSE`) or the updated
(`evaluate = TRUE`) FDboost model

## Author

David Ruegamer

## Examples

``` r
######## Example from \code{?FDboost}
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df = 2) + bolsc(T_A, df = 2),
               timeformula = ~ bbs(time, df = 4),
               numInt = "equal", family = QuantReg(),
               offset = NULL, offset_control = o_control(k_min = 9),
               data = viscosity, control=boost_control(mstop = 10, nu = 0.4))
#> Use a smooth offset.
               
# update nu
mod2 <- update(mod1, control=boost_control(nu = 1)) # mstop will stay the same
#> Use a smooth offset.
# update mstop
mod3 <- update(mod2, control=boost_control(mstop = 100)) # nu=1 does not get changed
#> Use a smooth offset.
mod4 <- update(mod1, formula = vis ~ 1 + bolsc(T_C, df = 2)) # drop one term
#> Use a smooth offset.
```
