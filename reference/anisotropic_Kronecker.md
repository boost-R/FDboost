# Kronecker product or row tensor product of two base-learners with anisotropic penalty

Kronecker product or row tensor product of two base-learners allowing
for anisotropic penalties. For the Kronecker product, `%A%` works in the
general case, `%A0%` for the special case where the penalty is zero in
one direction. For the row tensor product, `%Xa0%` works for the special
case where the penalty is zero in one direction.

## Usage

``` r
bl1 %A% bl2

bl1 %A0% bl2

bl1 %Xa0% bl2
```

## Arguments

- bl1:

  base-learner 1, e.g. `bbs(x1)`

- bl2:

  base-learner 2, e.g. `bbs(x2)`

## Value

An object of class `blg` (base-learner generator) with a `dpp` function
as for other
[`baselearners`](https://rdrr.io/pkg/mboost/man/baselearners.html).

## Details

When `%O%` is called with a specification of `df` in both base-learners,
e.g. `bbs(x1, df = df1) %O% bbs(t, df = df2)`, the global `df` for the
Kroneckered base-learner is computed as `df = df1 * df2`. And thus the
penalty has only one smoothness parameter lambda resulting in an
isotropic penalty, \$\$P = lambda \* \[(P1 o I) + (I o P2)\],\$\$ with
overall penalty \\P\\, Kronecker product \\o\\, marginal penalty
matrices \\P1, P2\\ and identity matrices \\I\\. (Currie et al. (2006)
introduced the generalized linear array model, which has a design matrix
that is composed of the Kronecker product of two marginal design
matrices, which was implemented in mboost as `%O%`. See Brockhaus et al.
(2015) for the application of array models to functional data.)

In contrast, a Kronecker product with anisotropic penalty is obtained by
`%A%`, which allows for a different amount of smoothness in the two
directions. For example `bbs(x1, df = df1) %A% bbs(t, df = df2)` results
in computing two different values for lambda for the two marginal design
matrices and a global value of lambda to adjust for the global `df`,
i.e. \$\$P = lambda \* \[(lambda1 \* P1 o I) + (I o lambda2 \*
P2)\],\$\$ with Kronecker product \\o\\, where \\lambda1\\ is computed
individually for \\df1\\ and \\P1\\, \\lambda2\\ is computed
individually for \\df2\\ and \\P2\\, and \\lambda\\ is computed such
that the global \\df\\ hold \\df = df1 \* df2\\. For the computation of
\\lambda1\\ and \\lambda2\\ weights specified in the model call can only
be used when the weights, are such that they are specified on the level
of rows and columns of the response matrix Y, e.g. resampling weights on
the level of rows of Y and integration weights on the columns of Y are
possible. If this the weights cannot be separated to blg1 and blg2 all
weights are set to 1 for the computation of \\lambda1\\ and \\lambda2\\
which implies that \\lambda1\\ and \\lambda2\\ are equal over folds of
`cvrisk`. The computation of the global \\lambda\\ considers the
specified `weights`, such the global \\df\\ are correct.

The operator `%A0%` treats the important special case where \\lambda1 =
0\\ or \\lambda2 = 0\\. In this case it suffices to compute the global
lambda and computation gets faster and arbitrary weights can be
specified. Consider \\lambda1 = 0\\ then the penalty becomes \$\$P =
lambda \* \[(1 \* P1 o I) + (I o lambda2 \* P2)\] = lambda \* lambda2 \*
(I o P2),\$\$ and only one global \\lambda\\ is computed which is then
\\lambda \* lambda2\\.

If the `formula` in `FDboost` contains base-learners connected by `%O%`,
`%A%` or `%A0%`, those effects are not expanded with `timeformula`,
allowing for model specifications with different effects in
time-direction.

`%Xa0%` computes like `%X%` the row tensor product of two base-learners,
with the difference that it sets the penalty for one direction to zero.
Thus, `%Xa0%` behaves to `%X%` analogously like `%A0%` to `%O%`.

## References

Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): The
functional linear array model. Statistical Modelling, 15(3), 279-300.

Currie, I.D., Durban, M. and Eilers P.H.C. (2006): Generalized linear
array models with applications to multidimensional smoothing. Journal of
the Royal Statistical Society, Series B-Statistical Methodology, 68(2),
259-280.

## Examples

``` r
 
######## Example for anisotropic penalty  
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## isotropic penalty, as timeformula is kroneckered to each effect using %O% 
## only for the smooth intercept %A0% is used, as 1-direction should not be penalized 
mod1 <- FDboost(vis ~ 1 + 
                bolsc(T_C, df = 1) + 
                bolsc(T_A, df = 1) + 
                bols(T_C, df = 1) %Xc% bols(T_A, df = 1),
                timeformula = ~ bbs(time, df = 3),
                numInt = "equal", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#> Use a smooth offset.
## cf. the formula that is passed to mboost
mod1$formulaMboost
#> [1] "dresponse ~ bols(ONEx, intercept = FALSE, df = 1) %A0% bbs(time,     df = 3) + bolsc(T_C, df = 1) %O% bbs(time, df = 3) + bolsc(T_A,     df = 1) %O% bbs(time, df = 3) + bols(T_C, df = 1) %Xc% bols(T_A,     df = 1) %O% bbs(time, df = 3)"

## anisotropic effects using %A0%, as lambda1 = 0 for all base-learners
## in this case using %A% gives the same model, but three lambdas are computed explicitly 
mod1a <- FDboost(vis ~ 1 + 
                bolsc(T_C, df = 1) %A0% bbs(time, df = 3) + 
                bolsc(T_A, df = 1) %A0% bbs(time, df = 3) + 
                bols(T_C, df = 1) %Xc% bols(T_A, df = 1) %A0% bbs(time, df = 3),
                timeformula = ~ bbs(time, df = 3),
                numInt = "equal", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = viscosity, control=boost_control(mstop = 100, nu = 0.4)) 
#> Use a smooth offset.
## cf. the formula that is passed to mboost
mod1a$formulaMboost
#> [1] "dresponse ~ bols(ONEx, intercept = FALSE, df = 1) %A0% bbs(time,     df = 3) + bolsc(T_C, df = 1) %A0% bbs(time, df = 3) + bolsc(T_A,     df = 1) %A0% bbs(time, df = 3) + bols(T_C, df = 1) %Xc% bols(T_A,     df = 1) %A0% bbs(time, df = 3)"

## alternative model specification by using a 0-matrix as penalty 
## only works for bolsc() as in bols() one cannot specify K 
## -> model without interaction term 
K0 <- matrix(0, ncol = 2, nrow = 2)
mod1k0 <- FDboost(vis ~ 1 + 
                 bolsc(T_C, df = 1, K = K0) +
                 bolsc(T_A, df = 1, K = K0), 
                 timeformula = ~ bbs(time, df = 3), 
                 numInt = "equal", family = QuantReg(), 
                 offset = NULL, offset_control = o_control(k_min = 9), 
                 data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#> Use a smooth offset.
## cf. the formula that is passed to mboost
mod1k0$formulaMboost
#> [1] "dresponse ~ bols(ONEx, intercept = FALSE, df = 1) %A0% bbs(time,     df = 3) + bolsc(T_C, df = 1, K = K0) %O% bbs(time, df = 3) +     bolsc(T_A, df = 1, K = K0) %O% bbs(time, df = 3)"
                
## optimize mstop for mod1, mod1a and mod1k0
## ...
                
## compare estimated coefficients
# \donttest{
if (interactive()) {
  oldpar <- par(mfrow=c(4, 2))
  plot(mod1, which = 1)
  plot(mod1a, which = 1)
  plot(mod1, which = 2)
  plot(mod1a, which = 2)
  plot(mod1, which = 3)
  plot(mod1a, which = 3)
  funplot(mod1$yind, predict(mod1, which=4))
  funplot(mod1$yind, predict(mod1a, which=4))
  par(oldpar)
}
# }
```
