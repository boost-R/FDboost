# Viscosity of resin over time

In an experimental setting the viscosity of resin was measured over time
to asses the curing process depending on 5 binary factors (low-high).

## Usage

``` r
data("viscosity")
```

## Format

A data list with 64 observations on the following 7 variables.

- `visAll`:

  viscosity measures over all available time points

- `timeAll`:

  time points of viscosity measures

- `T_C`:

  temperature of tools

- `T_A`:

  temperature of resin

- `T_B`:

  temperature of curing agent

- `rspeed`:

  rotational speed

- `mflow`:

  mass flow

## Details

The aim is to determine factors that affect the curing process in the
mold. The desired viscosity-curve has low values in the beginning
followed by a sharp increase. Due to technical reasons the measuring
method of the rheometer has to be changed in a certain range of
viscosity. The first observations are measured by rotation of a blade
giving observations every two seconds, the later observations are
measured through oscillation of a blade giving observations every ten
seconds. In the later observations the resin is quite hard so the
measurements should be interpreted as a qualitative measure of
hardening.

## Source

Wolfgang Raffelt, Technical University of Munich, Institute for Carbon
Composites

## Examples

``` r
 data("viscosity", package = "FDboost")
 ## set time-interval that should be modeled
 interval <- "101"

 ## model time until "interval" and take log() of viscosity
 end <- which(viscosity$timeAll==as.numeric(interval))
 viscosity$vis <- log(viscosity$visAll[,1:end])
 viscosity$time <- viscosity$timeAll[1:end]
 
 ## fit median regression model with 100 boosting iterations,
 ## step-length 0.4 and smooth time-specific offset
 ## the factors are in effect coding -1, 1 for the levels
 mod <- FDboost(vis ~ 1 + bols(T_C, contrasts.arg = "contr.sum", intercept=FALSE)
                + bols(T_A, contrasts.arg = "contr.sum", intercept=FALSE),
                timeformula=~bbs(time, lambda=100),
                numInt="equal", family=QuantReg(),
                offset=NULL, offset_control = o_control(k_min = 9),
                data=viscosity, control=boost_control(mstop = 100, nu = 0.4))
#> Use a smooth offset.
 summary(mod)
#> $object
#> 
#>   Model-based Boosting with Functional Response 
#> 
#> Call:
#> FDboost(formula = vis ~ 1 + bols(T_C, contrasts.arg = "contr.sum",     intercept = FALSE) + bols(T_A, contrasts.arg = "contr.sum",     intercept = FALSE), timeformula = ~bbs(time, lambda = 100),     numInt = "equal", data = viscosity, offset = NULL, offset_control = o_control(k_min = 9),     family = QuantReg(), control = boost_control(mstop = 100,         nu = 0.4))
#> 
#> 10
#>   Quantile Regression 
#> 
#> Loss function: tau * (y - f) * ((y - f) >= 0) - (1 - tau) * (y - f) * ((y -  
#>      f) < 0) 
#>  
#> 
#> Number of boosting iterations: mstop = 100 
#> Step size:  0.4 
#> Offset:  3.303 3.264 3.227 ... 5.371 5.502 5.649 
#> Number of baselearners:  3 
#> 
#> 
#> $selprob
#> "bols(T_C, intercept = FALSE, contrasts.arg = \\"contr.sum\\") %O% bbs(time, lambda = 100)" 
#>                                                                                        0.85 
#>            "bols(ONEx, intercept = FALSE, lambda = 0) %A0% bbs(time, df = 4, lambda = 100)" 
#>                                                                                        0.08 
#> "bols(T_A, intercept = FALSE, contrasts.arg = \\"contr.sum\\") %O% bbs(time, lambda = 100)" 
#>                                                                                        0.07 
#> 
#> attr(,"class")
#> [1] "summary.FDboost"
```
