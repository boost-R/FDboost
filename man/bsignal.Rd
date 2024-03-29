% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baselearners.R
\name{bsignal}
\alias{bsignal}
\alias{bconcurrent}
\alias{bhist}
\alias{bfpc}
\title{Base-learners for Functional Covariates}
\usage{
bsignal(
  x,
  s,
  index = NULL,
  inS = c("smooth", "linear", "constant"),
  knots = 10,
  boundary.knots = NULL,
  degree = 3,
  differences = 1,
  df = 4,
  lambda = NULL,
  center = FALSE,
  cyclic = FALSE,
  Z = NULL,
  penalty = c("ps", "pss"),
  check.ident = FALSE
)

bconcurrent(
  x,
  s,
  time,
  index = NULL,
  knots = 10,
  boundary.knots = NULL,
  degree = 3,
  differences = 1,
  df = 4,
  lambda = NULL,
  cyclic = FALSE
)

bhist(
  x,
  s,
  time,
  index = NULL,
  limits = "s<=t",
  standard = c("no", "time", "length"),
  intFun = integrationWeightsLeft,
  inS = c("smooth", "linear", "constant"),
  inTime = c("smooth", "linear", "constant"),
  knots = 10,
  boundary.knots = NULL,
  degree = 3,
  differences = 1,
  df = 4,
  lambda = NULL,
  penalty = c("ps", "pss"),
  check.ident = FALSE
)

bfpc(
  x,
  s,
  index = NULL,
  df = 4,
  lambda = NULL,
  penalty = c("identity", "inverse", "no"),
  pve = 0.99,
  npc = NULL,
  npc.max = 15,
  getEigen = TRUE
)
}
\arguments{
\item{x}{matrix of functional variable x(s). The functional covariate has to be 
supplied as n by <no. of evaluations> matrix, i.e., each row is one functional observation.}

\item{s}{vector for the index of the functional variable x(s) giving the 
measurement points of the functional covariate.}

\item{index}{a vector of integers for expanding the covariate in \code{x} 
For example, \code{bsignal(X, s, index = index)} is equal to \code{bsignal(X[index,], s)}, 
where index is an integer of length greater or equal to \code{NROW(x)}.}

\item{inS}{the functional effect can be smooth, linear or constant in s, 
which is the index of the functional covariates x(s).}

\item{knots}{either the number of knots or a vector of the positions 
of the interior knots (for more details see \code{\link[mboost:baselearners]{bbs}}).}

\item{boundary.knots}{boundary points at which to anchor the B-spline basis 
(default the range of the data). A vector (of length 2) 
for the lower and the upper boundary knot can be specified.}

\item{degree}{degree of the regression spline.}

\item{differences}{a non-negative integer, typically 1, 2 or 3. Defaults to 1.
If \code{differences} = \emph{k}, \emph{k}-th-order differences are used as 
a penalty (\emph{0}-th order differences specify a ridge penalty).}

\item{df}{trace of the hat matrix for the base-learner defining the 
base-learner complexity. Low values of \code{df} correspond to a 
large amount of smoothing and thus to "weaker" base-learners.}

\item{lambda}{smoothing parameter of the penalty, computed from \code{df} when \code{df} is specified.}

\item{center}{See \code{\link[mboost:baselearners]{bbs}}. 
The effect is re-parameterized such that the unpenalized part of the fit is subtracted and only 
the penalized effect is fitted, using a spectral decomposition of the penalty matrix.  
The unpenalized, parametric part has then to be included in separate 
base-learners using \code{bsignal(..., inS = 'constant')} or \code{bsignal(..., inS = 'linear')} 
for first (\code{difference = 1}) and second (\code{difference = 2}) order difference penalty respectively. 
See the help on the argument \code{center} of \code{\link[mboost:baselearners]{bbs}}.}

\item{cyclic}{if \code{cyclic = TRUE} the fitted coefficient function coincides at the boundaries 
(useful for cyclic covariates such as day time etc.).}

\item{Z}{a transformation matrix for the design-matrix over the index of the covariate.
\code{Z} can be calculated as the transformation matrix for a sum-to-zero constraint in the case
that all trajectories have the same mean 
(then a shift in the coefficient function is not identifiable).}

\item{penalty}{for \code{bsignal}, by default, \code{penalty = "ps"}, the difference penalty for P-splines is used, 
for \code{penalty = "pss"} the penalty matrix is transformed to have full rank, 
so called shrinkage approach by Marra and Wood (2011). 
For \code{bfpc} the penalty can be either \code{"identity"} for a ridge penalty 
(the default) or \code{"inverse"} to use the matrix with the inverse eigenvalues 
on the diagonal as penalty matrix or \code{"no"} for no penalty.}

\item{check.ident}{use checks for identifiability of the effect, based on Scheipl and Greven (2016) 
for linear functional effect using \code{bsignal} and 
based on Brockhaus et al. (2017) for historical effects using \code{bhist}}

\item{time}{vector for the index of the functional response y(time) 
giving the measurement points of the functional response.}

\item{limits}{defaults to \code{"s<=t"} for an historical effect with s<=t;
either one of \code{"s<t"} or \code{"s<=t"} for [l(t), u(t)] = [T1, t]; 
otherwise specify limits as a function for integration limits [l(t), u(t)]: 
function that takes \eqn{s} as the first and \code{t} as the second argument and returns 
\code{TRUE} for combinations of values (s,t) if \eqn{s} falls into the integration range for 
the given \eqn{t}.}

\item{standard}{the historical effect can be standardized with a factor. 
"no" means no standardization, "time" standardizes with the current value of time and 
"length" standardizes with the length of the integral}

\item{intFun}{specify the function that is used to compute integration weights in \code{s} 
over the functional covariate \eqn{x(s)}}

\item{inTime}{the historical effect can be smooth, linear or constant in time, 
which is the index of the functional response y(time).}

\item{pve}{proportion of variance explained by the first K functional principal components (FPCs): 
used to choose the number of functional principal components (FPCs).}

\item{npc}{prespecified value for the number K of FPCs (if given, this overrides \code{pve}).}

\item{npc.max}{maximal number K of FPCs to use; defaults to 15.}

\item{getEigen}{save the eigenvalues and eigenvectors, defaults to \code{TRUE}.}
}
\value{
Equally to the base-learners of package \code{mboost}: 

An object of class \code{blg} (base-learner generator) with a 
\code{dpp()} function (dpp, data pre-processing). 

The call of \code{dpp()} returns an object of class 
\code{bl} (base-learner) with a \code{fit()} function. The call to 
\code{fit()} finally returns an object of class \code{bm} (base-model).
}
\description{
Base-learners that fit effects of functional covariates.
}
\details{
\code{bsignal()} implements a base-learner for functional covariates to  
estimate an effect of the form \eqn{\int x_i(s)\beta(s)ds}. Defaults to a cubic  
B-spline basis with first difference penalties for \eqn{\beta(s)} and numerical 
integration over the entire range by using trapezoidal Riemann weights. 
If \code{bsignal()} is used within \code{FDboost()}, the base-learner of 
\code{timeformula} is attached, resulting in an effect varying over the index
of the response \eqn{\int x_i(s)\beta(s, t)ds} if \code{timeformula = bbs(t)}. 
The functional variable must be observed on one common grid \code{s}.  

\code{bconcurrent()} implements a concurrent effect for a functional covariate
on a functional response, i.e., an effect of the form \eqn{x_i(t)\beta(t)} for
a functional response \eqn{Y_i(t)} and concurrently observed covariate \eqn{x_i(t)}. 
\code{bconcurrent()} can only be used if \eqn{Y(t)} and \eqn{x(s)} are observed over
the same domain \eqn{s,t \in [T1, T2]}.  

\code{bhist()} implements a base-learner for functional covariates with 
flexible integration limits \code{l(t)}, \code{r(t)} and the possibility to
standardize the effect by \code{1/t} or the length of the integration interval. 
The effect is \eqn{stand * \int_{l(t)}^{r_{t}} x(s)\beta(t,s)ds}, where \eqn{stand} is 
the chosen standardization which defaults to 1. 
The base-learner defaults to a historical effect of the form 
\eqn{\int_{T1}^{t} x_i(s)\beta(t,s)ds}, 
where \eqn{T1} is the minimal index of \eqn{t} of the response \eqn{Y(t)}. 
The functional covariate must be observed on one common grid \code{s}.  
See Brockhaus et al. (2017) for details on historical effects.   

\code{bfpc()} is a base-learner for a linear effect of functional covariates based on 
functional principal component analysis (FPCA). 
For the functional linear effect \eqn{\int x_i(s)\beta(s)ds} the functional covariate 
and the coefficient function are both represented by a FPC basis. 
The functional covariate
\eqn{x(s)} is decomposed into \eqn{x(s) \approx \sum_{k=1}^K \xi_{ik} \Phi_k(s)} using 
\code{\link[refund]{fpca.sc}} for the truncated Karhunen-Loeve decomposition. 
Then \eqn{\beta(s)} is represented in the function
space spanned by \eqn{\Phi_k(s)}, k=1,...,K, see Scheipl et al. (2015) for details. 
As penalty matrix, the identity matrix is used. 
The implementation is similar to \code{\link[refund]{ffpc}}.  

It is recommended to use centered functional covariates with 
\eqn{\sum_i x_i(s) = 0} for all \eqn{s} in \code{bsignal()}-, 
\code{bhist()}- and \code{bconcurrent()}-terms. 
For centered covariates, the effects are centered per time-point of the response. 
If all effects are centered, the functional intercept 
can be interpreted as the global mean function. 

The base-learners for functional covariates cannot deal with any missing 
values in the covariates.
}
\examples{
######## Example for scalar-on-function-regression with bsignal()  
data("fuelSubset", package = "FDboost")

## center the functional covariates per observed wavelength
fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)

## to make mboost:::df2lambda() happy (all design matrix entries < 10)
## reduce range of argvals to [0,1] to get smaller integration weights
fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
                                  (max(uvvis.lambda) - min(uvvis.lambda) ))
fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
                                (max(nir.lambda) - min(nir.lambda) ))

## model fit with scalar response and two functional linear effects 
## include no intercept 
## as all base-learners are centered around 0 
mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df=4, check.ident = FALSE), 
               timeformula = NULL, data = fuelSubset) 
summary(mod2) 

 
###############################################
### data simulation like in manual of pffr::ff

if(require(refund)){

#########
# model with linear functional effect, use bsignal()
# Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
set.seed(2121)
data1 <- pffrSim(scenario = "ff", n = 40)
data1$X1 <- scale(data1$X1, scale = FALSE)
dat_list <- as.list(data1)
dat_list$t <- attr(data1, "yindex")
dat_list$s <- attr(data1, "xindex")

## model fit by FDboost 
m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 5), 
              timeformula = ~ bbs(t, knots = 5), data = dat_list, 
              control = boost_control(mstop = 21))

## search optimal mSTOP
\donttest{
  set.seed(123)
  cv <- validateFDboost(m1, grid = 1:100) # 21 iterations
}

## model fit by pffr
t <- attr(data1, "yindex")
s <- attr(data1, "xindex")
m1_pffr <- pffr(Y ~ ff(X1, xind = s), yind = t, data = data1)

\donttest{
  oldpar <- par(mfrow = c(2, 2))
  plot(m1, which = 1); plot(m1, which = 2) 
  plot(m1_pffr, select = 1, shift = m1_pffr$coefficients["(Intercept)"]) 
  plot(m1_pffr, select = 2)
  par(oldpar)
}


############################################
# model with functional historical effect, use bhist() 
# Y(t) = f(t)  + \int_0^t X1(s)\beta(s,t)ds + eps
set.seed(2121)
mylimits <- function(s, t){
  (s < t) | (s == t)
}
data2 <- pffrSim(scenario = "ff", n = 40, limits = mylimits)
data2$X1 <- scale(data2$X1, scale = FALSE)
dat2_list <- as.list(data2)
dat2_list$t <- attr(data2, "yindex")
dat2_list$s <- attr(data2, "xindex")

## model fit by FDboost 
m2 <- FDboost(Y ~ 1 + bhist(x = X1, s = s, time = t, knots = 5), 
              timeformula = ~ bbs(t, knots = 5), data = dat2_list, 
              control = boost_control(mstop = 40))
              
## search optimal mSTOP
\donttest{
  set.seed(123)
  cv2 <- validateFDboost(m2, grid = 1:100) # 40 iterations
}               

## model fit by pffr
t <- attr(data2, "yindex")
s <- attr(data2, "xindex")
m2_pffr <- pffr(Y ~ ff(X1, xind = s, limits = "s<=t"), yind = t, data = data2)

\donttest{
oldpar <- par(mfrow = c(2, 2))
plot(m2, which = 1); plot(m2, which = 2)
## plot of smooth intercept does not contain m1_pffr$coefficients["(Intercept)"]
plot(m2_pffr, select = 1, shift = m2_pffr$coefficients["(Intercept)"]) 
plot(m2_pffr, select = 2) 
par(oldpar)
}


}


}
\references{
Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
The functional linear array model. Statistical Modelling, 15(3), 279-300.

Brockhaus, S., Melcher, M., Leisch, F. and Greven, S. (2017): 
Boosting flexible functional regression models with a high number of functional historical effects,  
Statistics and Computing, 27(4), 913-926.   

Marra, G. and Wood, S.N. (2011): Practical variable selection for generalized additive models. 
Computational Statistics & Data Analysis, 55, 2372-2387.

Scheipl, F., Staicu, A.-M. and Greven, S. (2015): 
Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 24(2), 477-501. 

Scheipl, F. and Greven, S. (2016): Identifiability in penalized function-on-function regression models. 
Electronic Journal of Statistics, 10(1), 495-526.
}
\seealso{
\code{\link{FDboost}} for the model fit.
}
\keyword{models}
