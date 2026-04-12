# Stability Selection

Function for stability selection with functional response. Per default
the sampling is done on the level of curves and if the model contains a
smooth functional intercept, this intercept is refittedn in each
sampling fold.

## Usage

``` r
# S3 method for class 'FDboost'
stabsel(
  x,
  refitSmoothOffset = TRUE,
  cutoff,
  q,
  PFER,
  folds = cvLong(x$id, weights = rep(1, l = length(x$id)), type = "subsampling", B = B),
  B = ifelse(sampling.type == "MB", 100, 50),
  assumption = c("unimodal", "r-concave", "none"),
  sampling.type = c("SS", "MB"),
  papply = mclapply,
  verbose = TRUE,
  eval = TRUE,
  ...
)
```

## Arguments

- x:

  fitted FDboost-object

- refitSmoothOffset:

  logical, should the offset be refitted in each learning sample?
  Defaults to `TRUE`.

- cutoff:

  cutoff between 0.5 and 1. Preferably a value between 0.6 and 0.9
  should be used.

- q:

  number of (unique) selected variables (or groups of variables
  depending on the model) that are selected on each subsample.

- PFER:

  upper bound for the per-family error rate. This specifies the amount
  of falsely selected base-learners, which is tolerated. See details of
  [`stabsel`](https://rdrr.io/pkg/mboost/man/stabsel.html).

- folds:

  a weight matrix with number of rows equal to the number of
  observations, see `{cvLong}`. Usually one should not change the
  default here as subsampling with a fraction of 1/2 is needed for the
  error bounds to hold. One usage scenario where specifying the folds by
  hand might be the case when one has dependent data (e.g. clusters) and
  thus wants to draw clusters (i.e., multiple rows together) not
  individuals.

- B:

  number of subsampling replicates. Per default, we use 50 complementary
  pairs for the error bounds of Shah & Samworth (2013) and 100 for the
  error bound derived in Meinshausen & Buehlmann (2010). As we use `B`
  complementary pairs in the former case this leads to `2B` subsamples.

- assumption:

  Defines the type of assumptions on the distributions of the selection
  probabilities and simultaneous selection probabilities. Only
  applicable for `sampling.type = "SS"`. For `sampling.type = "MB"` we
  always use `"none"`.

- sampling.type:

  use sampling scheme of of Shah & Samworth (2013), i.e., with
  complementary pairs (`sampling.type = "SS"`), or the original sampling
  scheme of Meinshausen & Buehlmann (2010).

- papply:

  (parallel) apply function, defaults to mclapply. Alternatively,
  parLapply can be used. In the latter case, usually more setup is
  needed (see example of cvrisk for some details).

- verbose:

  logical (default: TRUE) that determines wether warnings should be
  issued.

- eval:

  logical. Determines whether stability selection is evaluated
  (`eval = TRUE`; default) or if only the parameter combination is
  returned.

- ...:

  additional arguments to
  [`cvrisk`](https://rdrr.io/pkg/mboost/man/cvrisk.html) or
  [`validateFDboost`](https://boost-R.github.io/FDboost/reference/validateFDboost.md).

## Value

An object of class `stabsel` with a special print method. For the
elements of the object, see
[`stabsel`](https://rdrr.io/pkg/mboost/man/stabsel.html)

## Details

The number of boosting iterations is an important hyper-parameter of the
boosting algorithms and can be chosen using the functions
`cvrisk.FDboost` and `validateFDboost` as they compute honest, i.e.
out-of-bag, estimates of the empirical risk for different numbers of
boosting iterations. The weights (zero weights correspond to test cases)
are defined via the folds matrix, see
[`cvrisk`](https://rdrr.io/pkg/mboost/man/cvrisk.html) in package
mboost. See Hofner et al. (2015) for the combination of stability
selection and component-wise boosting.

## References

B. Hofner, L. Boccuto and M. Goeker (2015), Controlling false
discoveries in high-dimensional situations: boosting with stability
selection. BMC Bioinformatics, 16, 1-17.

N. Meinshausen and P. Buehlmann (2010), Stability selection. Journal of
the Royal Statistical Society, Series B, 72, 417-473.

R.D. Shah and R.J. Samworth (2013), Variable selection with error
control: another look at stability selection. Journal of the Royal
Statistical Society, Series B, 75, 55-80.

## See also

[`stabsel`](https://rdrr.io/pkg/mboost/man/stabsel.html) to perform
stability selection for a mboost-object.

## Examples

``` r
######## Example for function-on-scalar-regression
data("viscosity", package = "FDboost")
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## fit a model cotaining all main effects 
modAll <- FDboost(vis ~ 1 
          + bolsc(T_C, df=1) %A0% bbs(time, df=5) 
          + bolsc(T_A, df=1) %A0% bbs(time, df=5)
          + bolsc(T_B, df=1) %A0% bbs(time, df=5)
          + bolsc(rspeed, df=1) %A0% bbs(time, df=5)
          + bolsc(mflow, df=1) %A0% bbs(time, df=5), 
       timeformula = ~bbs(time, df=5), 
       numInt = "Riemann", family = QuantReg(), 
       offset = NULL, offset_control = o_control(k_min = 10),
       data = viscosity, 
       control = boost_control(mstop = 100, nu = 0.2))
#> Use a smooth offset.


## create folds for stability selection  
## only 5 folds for a fast example, usually use 50 folds 
set.seed(1911)
folds <- cvLong(modAll$id, weights = rep(1, l = length(modAll$id)), 
                type = "subsampling", B = 5) 
    
# \donttest{        
## stability selection with refit of the smooth intercept 
stabsel_parameters(q = 3, PFER = 1, p = 6, sampling.type = "SS")
#> Stability selection with unimodality assumption
#> 
#> Cutoff: 0.84; q: 3; PFER (*):  1 
#>    (*) or expected number of low selection probability variables
#> PFER (specified upper bound):  1 
#> PFER corresponds to signif. level 0.167 (without multiplicity adjustment)
sel1 <- stabsel(modAll, q = 3, PFER = 1, folds = folds, grid = 1:200, sampling.type = "SS")
#> Use applyFolds() to recompute the smooth offset in each fold.
sel1
#>  Stability Selection with unimodality assumption
#> 
#> Selected variables:
#> "bols(ONEx, intercept = FALSE, df = 1) %A0% bbs(time, df = 5)" 
#>                                                              1 
#>                    "bolsc(T_C, df = 1) %A0% bbs(time, df = 5)" 
#>                                                              2 
#> 
#> Selection probabilities:
#>                 "bolsc(rspeed, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.0 
#>                  "bolsc(mflow, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.0 
#>                    "bolsc(T_B, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.1 
#>                    "bolsc(T_A, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.9 
#> "bols(ONEx, intercept = FALSE, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            1.0 
#>                    "bolsc(T_C, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            1.0 
#> 
#> ---
#> Cutoff: 1; q: 3; PFER (*):  0.5 
#>    (*) or expected number of low selection probability variables
#> PFER corresponds to signif. level 0.0833 (without multiplicity adjustment)
#> 

## stability selection without refit of the smooth intercept 
sel2 <- stabsel(modAll, refitSmoothOffset = FALSE, q = 3, PFER = 1, 
                folds = folds, grid = 1:200, sampling.type = "SS")
#> The smooth offset is fixed over all folds.
#> For bolsc, the transformation matrix Z is fixed over all folds.
sel2
#>  Stability Selection with unimodality assumption
#> 
#> Selected variables:
#> "bolsc(T_C, df = 1) %A0% bbs(time, df = 5)" 
#>                                           2 
#> "bolsc(T_A, df = 1) %A0% bbs(time, df = 5)" 
#>                                           3 
#> 
#> Selection probabilities:
#>                    "bolsc(T_B, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.0 
#>                 "bolsc(rspeed, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.0 
#>                  "bolsc(mflow, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.1 
#> "bols(ONEx, intercept = FALSE, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            0.9 
#>                    "bolsc(T_C, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            1.0 
#>                    "bolsc(T_A, df = 1) %A0% bbs(time, df = 5)" 
#>                                                            1.0 
#> 
#> ---
#> Cutoff: 1; q: 3; PFER (*):  0.5 
#>    (*) or expected number of low selection probability variables
#> PFER corresponds to signif. level 0.0833 (without multiplicity adjustment)
#> 
# }
```
