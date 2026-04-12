# Constrained row tensor product

Combining single base-learners to form new, more complex base-learners,
with an identifiability constraint to center the interaction around the
intercept and around the two main effects. Suitable for functional
response.

## Usage

``` r
bl1 %Xc% bl2
```

## Arguments

- bl1:

  base-learner 1, e.g. `bols(x1)`

- bl2:

  base-learner 2, e.g. `bols(x2)`

## Value

An object of class `blg` (base-learner generator) with a `dpp` function
as for other
[`baselearners`](https://rdrr.io/pkg/mboost/man/baselearners.html).

## Details

Similar to `%X%` in package `mboost`, see
[`%X%`](https://rdrr.io/pkg/mboost/man/baselearners.html), a row tensor
product of linear base-learners is returned by `%Xc%`. `%Xc%` applies a
sum-to-zero constraint to the design matrix suitable for functional
response if an interaction of two scalar covariates is specified in the
case that the model contains a global intercept and both main effects,
as the interaction is centered around the intercept and centered around
the two main effects. See Web Appendix A of Brockhaus et al. (2015) for
details on how to enforce the constraint for the functional intercept.
Use, e.g., in a model call to `FDboost`, following the scheme,
`y ~ 1 + bolsc(x1) + bolsc(x2) + bols(x1) %Xc% bols(x2)`, where `1`
induces a global intercept and `x1`, `x2` are factor variables, see
Ruegamer et al. (2018).

## References

Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): The
functional linear array model. Statistical Modelling, 15(3), 279-300.

Ruegamer D., Brockhaus, S., Gentsch K., Scherer, K., Greven, S. (2018).
Boosting factor-specific functional historical models for the detection
of synchronization in bioelectrical signals. Journal of the Royal
Statistical Society: Series C (Applied Statistics), 67, 621-642.

## Author

Sarah Brockhaus, David Ruegamer

## Examples

``` r
 
######## Example for function-on-scalar-regression with interaction effect of two scalar covariates 
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## fit model with interaction that is centered around the intercept 
## and the two main effects 
mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df=1) + bolsc(T_A, df=1) + 
                bols(T_C, df=1) %Xc% bols(T_A, df=1),
                timeformula = ~bbs(time, df=6),
                numInt = "equal", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#> Use a smooth offset.
                
## check centering around intercept
colMeans(predict(mod1, which = 4))
#>  [1]  1.387779e-17  1.040834e-17  1.040834e-17  1.734723e-17  1.040834e-17
#>  [6]  1.040834e-17  5.204170e-18  5.204170e-18  3.469447e-18  2.602085e-18
#> [11]  2.602085e-18  2.602085e-18  2.602085e-18  2.168404e-18  6.505213e-19
#> [16] -4.336809e-19 -4.336809e-19 -1.734723e-18 -1.734723e-18 -1.734723e-18
#> [21] -4.336809e-18 -3.469447e-18 -3.469447e-18 -3.469447e-18 -3.469447e-18
#> [26] -5.204170e-18 -5.204170e-18 -5.204170e-18 -3.469447e-18 -6.938894e-18
#> [31] -6.938894e-18 -8.673617e-18 -6.938894e-18 -6.938894e-18 -3.469447e-18
#> [36] -6.938894e-18 -1.040834e-17 -1.040834e-17 -1.387779e-17 -1.040834e-17
#> [41] -1.387779e-17 -1.734723e-17 -1.040834e-17 -1.387779e-17 -1.040834e-17
#> [46] -1.040834e-17

## check centering around main effects
colMeans(predict(mod1, which = 4)[viscosity$T_A == "low", ])
#>  [1] -2.081668e-17 -2.081668e-17 -3.469447e-17 -2.081668e-17 -2.775558e-17
#>  [6] -2.081668e-17 -1.734723e-17 -1.040834e-17 -5.204170e-18 -5.204170e-18
#> [11] -5.204170e-18 -3.469447e-18 -5.204170e-18 -4.336809e-18 -2.168404e-18
#> [16]  1.084202e-18  2.602085e-18  3.469447e-18  5.204170e-18  6.938894e-18
#> [21]  5.204170e-18  8.673617e-18  1.214306e-17  1.040834e-17  1.040834e-17
#> [26]  1.040834e-17  1.387779e-17  1.040834e-17  1.734723e-17  1.387779e-17
#> [31]  1.387779e-17  1.734723e-17  1.734723e-17  2.428613e-17  2.775558e-17
#> [36]  2.081668e-17  2.775558e-17  2.775558e-17  2.775558e-17  3.469447e-17
#> [41]  4.163336e-17  2.775558e-17  3.469447e-17  2.775558e-17  3.469447e-17
#> [46]  2.081668e-17
colMeans(predict(mod1, which = 4)[viscosity$T_A == "high", ])
#>  [1]  4.857226e-17  4.163336e-17  5.551115e-17  5.551115e-17  4.857226e-17
#>  [6]  4.163336e-17  2.775558e-17  2.081668e-17  1.214306e-17  1.040834e-17
#> [11]  1.040834e-17  8.673617e-18  1.040834e-17  8.673617e-18  3.469447e-18
#> [16] -1.951564e-18 -3.469447e-18 -6.938894e-18 -8.673617e-18 -1.040834e-17
#> [21] -1.387779e-17 -1.561251e-17 -1.908196e-17 -1.734723e-17 -1.734723e-17
#> [26] -2.081668e-17 -2.428613e-17 -2.081668e-17 -2.428613e-17 -2.775558e-17
#> [31] -2.775558e-17 -3.469447e-17 -3.122502e-17 -3.816392e-17 -3.469447e-17
#> [36] -3.469447e-17 -4.857226e-17 -4.857226e-17 -5.551115e-17 -5.551115e-17
#> [41] -6.938894e-17 -6.245005e-17 -5.551115e-17 -5.551115e-17 -5.551115e-17
#> [46] -4.163336e-17
colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0
colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
#>  [1] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
#> [39] 0 0 0 0 0 0 0 0

## find optimal mstop using cvrsik() or validateFDboost()
## ... 

## look at interaction effect in one plot
# funplot(mod1$yind, predict(mod1, which=4))
```
