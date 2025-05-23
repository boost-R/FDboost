% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossvalidation.R
\name{validateFDboost}
\alias{validateFDboost}
\title{Cross-Validation and Bootstrapping over Curves}
\usage{
validateFDboost(
  object,
  response = NULL,
  folds = cv(rep(1, length(unique(object$id))), type = "bootstrap"),
  grid = 1:mstop(object),
  fun = NULL,
  getCoefCV = TRUE,
  riskopt = c("mean", "median"),
  mrdDelete = 0,
  refitSmoothOffset = TRUE,
  showProgress = TRUE,
  ...
)
}
\arguments{
\item{object}{fitted FDboost-object}

\item{response}{optional, specify a response vector for the computation of the prediction errors.
Defaults to \code{NULL} which means that the response of the fitted model is used.}

\item{folds}{a weight matrix with number of rows equal to the number of observed trajectories.}

\item{grid}{the grid over which the optimal number of boosting iterations (mstop) is searched.}

\item{fun}{if \code{fun} is \code{NULL}, the out-of-bag risk is returned. 
\code{fun}, as a function of \code{object}, 
may extract any other characteristic of the cross-validated models. These are returned as is.}

\item{getCoefCV}{logical, defaults to \code{TRUE}. Should the coefficients and predictions
be computed for all the models on the sampled data?}

\item{riskopt}{how is the optimal stopping iteration determined. Defaults to the mean, 
but median is possible as well.}

\item{mrdDelete}{Delete values that are \code{mrdDelete} percent smaller than the mean
of the response. Defaults to 0 which means that only response values being 0 
are not used in the calculation of the MRD (= mean relative deviation).}

\item{refitSmoothOffset}{logical, should the offset be refitted in each learning sample? 
Defaults to \code{TRUE}. In \code{\link[mboost]{cvrisk}} the offset of the original model fit in  
\code{object} is used in all folds.}

\item{showProgress}{logical, defaults to \code{TRUE}.}

\item{...}{further arguments passed to \code{\link{mclapply}}}
}
\value{
The function \code{validateFDboost} returns a \code{validateFDboost}-object, 
which is a named list containing: 
\item{response}{the response}
\item{yind}{the observation points of the response}
\item{id}{the id variable of the response}
\item{folds}{folds that were used}
\item{grid}{grid of possible numbers of boosting iterations}
\item{coefCV}{if \code{getCoefCV} is \code{TRUE} the estimated coefficient functions in the folds}
\item{predCV}{if \code{getCoefCV} is \code{TRUE} the out-of-bag predicted values of the response}
\item{oobpreds}{if the type of folds is curves the out-of-bag predictions for each trajectory}
\item{oobrisk}{the out-of-bag risk}
\item{oobriskMean}{the out-of-bag risk at the minimal mean risk}
\item{oobmse}{the out-of-bag mean squared error (MSE)}
\item{oobrelMSE}{the out-of-bag relative mean squared error (relMSE)}
\item{oobmrd}{the out-of-bag mean relative deviation (MRD)}
\item{oobrisk0}{the out-of-bag risk without consideration of integration weights}
\item{oobmse0}{the out-of-bag mean squared error (MSE) without consideration of integration weights}
\item{oobmrd0}{the out-of-bag mean relative deviation (MRD) without consideration of integration weights}
\item{format}{one of "FDboostLong" or "FDboost" depending on the class of the object}
\item{fun_ret}{list of what fun returns if fun was specified}
}
\description{
DEPRECATED! 
The function \code{validateFDboost()} is deprecated,  
use \code{\link{applyFolds}} and \code{\link{bootstrapCI}} instead.
}
\details{
The number of boosting iterations is an important hyper-parameter of boosting  
and can be chosen using the function \code{validateFDboost} as they compute
honest, i.e., out-of-bag, estimates of the empirical risk for different numbers of boosting iterations. 

The function \code{validateFDboost} is especially suited to models with functional response. 
Using the option \code{refitSmoothOffset} the offset is refitted on each fold. 
Note, that the function \code{validateFDboost} expects folds that give weights
per curve without considering integration weights. The integration weights of 
\code{object} are used to compute the empirical risk as integral. The argument \code{response} 
can be useful in simulation studies where the true value of the response is known but for 
the model fit the response is used with noise.
}
\examples{
\donttest{
if(require(fda)){
 ## load the data
 data("CanadianWeather", package = "fda")
 
 ## use data on a daily basis 
 canada <- with(CanadianWeather, 
                list(temp = t(dailyAv[ , , "Temperature.C"]),
                     l10precip = t(dailyAv[ , , "log10precip"]),
                     l10precip_mean = log(colMeans(dailyAv[ , , "Precipitation.mm"]), base = 10),
                     lat = coordinates[ , "N.latitude"],
                     lon = coordinates[ , "W.longitude"],
                     region = factor(region),
                     place = factor(place),
                     day = 1:365,  ## corresponds to t: evaluation points of the fun. response 
                     day_s = 1:365))  ## corresponds to s: evaluation points of the fun. covariate
 
## center temperature curves per day 
canada$tempRaw <- canada$temp
canada$temp <- scale(canada$temp, scale = FALSE) 
rownames(canada$temp) <- NULL ## delete row-names 
  
## fit the model  
mod <- FDboost(l10precip ~ 1 + bolsc(region, df = 4) + 
                 bsignal(temp, s = day_s, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
               timeformula = ~ bbs(day, cyclic = TRUE, boundary.knots = c(0.5, 365.5)), 
               data = canada)
mod <- mod[75]

  #### create folds for 3-fold bootstrap: one weight for each curve
  set.seed(124)
  folds_bs <- cv(weights = rep(1, mod$ydim[1]), type = "bootstrap", B = 3)

  ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
  cvr <- applyFolds(mod, folds = folds_bs, grid = 1:75)

  ## compute out-of-bag risk and coefficient estimates on folds  
  cvr2 <- validateFDboost(mod, folds = folds_bs, grid = 1:75)

  ## weights per observation point  
  folds_bs_long <- folds_bs[rep(seq_len(nrow(folds_bs)), times = mod$ydim[2]), ]
  attr(folds_bs_long, "type") <- "3-fold bootstrap"
  ## compute out-of-bag risk on the 3 folds for 1 to 75 boosting iterations  
  cvr3 <- cvrisk(mod, folds = folds_bs_long, grid = 1:75)

  ## plot the out-of-bag risk
  oldpar <- par(mfrow = c(1,3))
  plot(cvr); legend("topright", lty=2, paste(mstop(cvr)))
  plot(cvr2)
  plot(cvr3); legend("topright", lty=2, paste(mstop(cvr3)))

  ## plot the estimated coefficients per fold
  ## more meaningful for higher number of folds, e.g., B = 100 
  par(mfrow = c(2,2))
  plotPredCoef(cvr2, terms = FALSE, which = 1)
  plotPredCoef(cvr2, terms = FALSE, which = 3)
  
  ## compute out-of-bag risk and predictions for leaving-one-curve-out cross-validation
  cvr_jackknife <- validateFDboost(mod, folds = cvLong(unique(mod$id), 
                                   type = "curves"), grid = 1:75)
  plot(cvr_jackknife)
  ## plot oob predictions per fold for 3rd effect 
  plotPredCoef(cvr_jackknife, which = 3) 
  ## plot coefficients per fold for 2nd effect
  plotPredCoef(cvr_jackknife, which = 2, terms = FALSE)
  
  par(oldpar)

}
}

}
