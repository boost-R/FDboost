pkgname <- "FDboost"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('FDboost')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("FDboost")
### * FDboost

flush(stderr()); flush(stdout())

### Name: FDboost
### Title: Model-based Gradient Boosting for Functional Response
### Aliases: FDboost
### Keywords: models nonlinear regression smooth

### ** Examples

######## Example for function-on-scalar-regression 
data("viscosity", package = "FDboost") 
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## fit median regression model with 100 boosting iterations,
## step-length 0.4 and smooth time-specific offset
## the factors are coded such that the effects are zero for each timepoint t
## no integration weights are used!
mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df = 2) + bolsc(T_A, df = 2),
               timeformula = ~ bbs(time, df = 4),
               numInt = "equal", family = QuantReg(),
               offset = NULL, offset_control = o_control(k_min = 9),
               data = viscosity, control=boost_control(mstop = 100, nu = 0.4))


######## Example for scalar-on-function-regression 
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

## model fit with scalar response 
## include no intercept as all base-learners are centered around 0
mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE), 
               timeformula = NULL, data = fuelSubset, control = boost_control(mstop = 200)) 
               
## additionally include a non-linear effect of the scalar variable h2o 
mod2s <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bbs(h2o, df = 4), 
               timeformula = NULL, data = fuelSubset, control = boost_control(mstop = 200)) 
               
## alternative model fit as FLAM model with scalar response; as timeformula = ~ bols(1)  
## adds a penalty over the index of the response, i.e., here a ridge penalty
## thus, mod2f and mod2 have different penalties 
mod2f <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df = 4, check.ident = FALSE), 
               timeformula = ~ bols(1), data = fuelSubset, control = boost_control(mstop = 200))
               

## Example for function-on-function-regression 
if(require(fda)){

  data("CanadianWeather", package = "fda")
  CanadianWeather$l10precip <- t(log(CanadianWeather$monthlyPrecip))
  CanadianWeather$temp <- t(CanadianWeather$monthlyTemp)
  CanadianWeather$region <- factor(CanadianWeather$region)
  CanadianWeather$month.s <- CanadianWeather$month.t <- 1:12
  
  ## center the temperature curves per time-point
  CanadianWeather$temp <- scale(CanadianWeather$temp, scale = FALSE)
  rownames(CanadianWeather$temp) <- NULL ## delete row-names
  
  ## fit model with cyclic splines over the year
  mod3 <- FDboost(l10precip ~ bols(region, df = 2.5, contrasts.arg = "contr.dummy") 
                   + bsignal(temp, month.s, knots = 11, cyclic = TRUE, 
                             df = 2.5, boundary.knots = c(0.5,12.5), check.ident = FALSE), 
                  timeformula = ~ bbs(month.t, knots = 11, cyclic = TRUE, 
                                      df = 3, boundary.knots = c(0.5, 12.5)), 
                  offset = "scalar", offset_control = o_control(k_min = 5), 
                  control = boost_control(mstop = 60), 
                  data = CanadianWeather) 
 
}

######## Example for functional response observed on irregular grid
######## Delete part of observations in viscosity data-set
data("viscosity", package = "FDboost")
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))

## only keep one eighth of the observation points
set.seed(123)
selectObs <- sort(sample(x = 1:(64*46), size = 64*46/4, replace = FALSE))
dataIrregular <- with(viscosity, list(vis = c(vis)[selectObs], 
                                      T_A = T_A, T_C = T_C,  
                                      time = rep(time, each = 64)[selectObs], 
                                      id = rep(1:64, 46)[selectObs]))

## fit median regression model with 50 boosting iterations,
## step-length 0.4 and smooth time-specific offset
## the factors are in effect coding -1, 1 for the levels
## no integration weights are used!
mod4 <- FDboost(vis ~ 1 + bols(T_C, contrasts.arg = "contr.sum", intercept = FALSE)
                + bols(T_A, contrasts.arg = "contr.sum", intercept=FALSE),
                timeformula = ~ bbs(time, lambda = 100), id = ~id, 
                numInt = "Riemann", family = QuantReg(),
                offset = NULL, offset_control = o_control(k_min = 9),
                data = dataIrregular, control = boost_control(mstop = 50, nu = 0.4))
## summary(mod4)
## plot(mod4)
## plotPredicted(mod4, lwdPred = 2)


## Be careful if you want to predict newdata with irregular response,  
## as the argument index is not considered in the prediction of newdata. 
## Thus, all covariates have to be repeated according to the number of observations 
## in each response trajectroy. 
## Predict four response curves with full time-observations 
## for the four combinations of T_A and T_C. 
newd <- list(T_A = factor(c(1,1,2,2), levels = 1:2, 
                        labels = c("low", "high"))[rep(1:4, length(viscosity$time))], 
             T_C = factor(c(1,2,1,2), levels = 1:2, 
                        labels = c("low", "high"))[rep(1:4, length(viscosity$time))], 
             time = rep(viscosity$time, 4))
             
pred <- predict(mod4, newdata = newd)
## funplot(x = rep(viscosity$time, 4), y = pred, id = rep(1:4, length(viscosity$time)))
                  
                



cleanEx()
nameEx("FDboostLSS")
### * FDboostLSS

flush(stderr()); flush(stdout())

### Name: FDboostLSS
### Title: Model-based Gradient Boosting for Functional GAMLSS
### Aliases: FDboostLSS
### Keywords: models nonlinear regression smooth

### ** Examples

########### simulate Gaussian scalar-on-function data
n <- 500 ## number of observations
G <- 120 ## number of observations per functional covariate
set.seed(123) ## ensure reproducibility
z <- runif(n) ## scalar covariate
z <- z - mean(z)
s <- seq(0, 1, l=G) ## index of functional covariate
## generate functional covariate
if(require(splines)){
   x <- t(replicate(n, drop(bs(s, df = 5, int = TRUE) %*% runif(5, min = -1, max = 1))))
}else{
  x <- matrix(rnorm(n*G), ncol = G, nrow = n)
}
x <- scale(x, center = TRUE, scale = FALSE) ## center x per observation point

mu <- 2 + 0.5*z + (1/G*x) %*% sin(s*pi)*5 ## true functions for expectation
sigma <- exp(0.5*z - (1/G*x) %*% cos(s*pi)*2) ## for standard deviation

y <- rnorm(mean = mu, sd = sigma, n = n) ## draw respone y_i ~ N(mu_i, sigma_i)

## save data as list containing s as well 
dat_list <- list(y = y, z = z, x = I(x), s = s)

## model fit with noncyclic algorithm assuming Gaussian location scale model 
m_boost <- FDboostLSS(list(mu = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16), 
                           sigma = y ~ bols(z, df = 2) + bsignal(x, s, df = 2, knots = 16)), 
                           timeformula = NULL, data = dat_list, method = "noncyclic")
summary(m_boost)




cleanEx()
nameEx("anisotropic_Kronecker")
### * anisotropic_Kronecker

flush(stderr()); flush(stdout())

### Name: anisotropic_Kronecker
### Title: Kronecker product or row tensor product of two base-learners
###   with anisotropic penalty
### Aliases: anisotropic_Kronecker %A% %A0% %Xa0%

### ** Examples

 
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
## cf. the formula that is passed to mboost
mod1$formulaMboost

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
## cf. the formula that is passed to mboost
mod1a$formulaMboost

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
## cf. the formula that is passed to mboost
mod1k0$formulaMboost
                
## optimize mstop for mod1, mod1a and mod1k0
## ...
                
## compare estimated coefficients




cleanEx()
nameEx("applyFolds")
### * applyFolds

flush(stderr()); flush(stdout())

### Name: applyFolds
### Title: Cross-Validation and Bootstrapping over Curves
### Aliases: applyFolds cvMa cvLong cvrisk.FDboost

### ** Examples

Ytest <- matrix(rnorm(15), ncol = 3) # 5 trajectories, each with 3 observations 
Ylong <- as.vector(Ytest)
## 4-folds for bootstrap for the response in long format without integration weights
cvMa(ydim = c(5,3), type = "bootstrap", B = 4)  
cvLong(id = rep(1:5, times = 3), type = "bootstrap", B = 4)

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



}




cleanEx()
nameEx("bbsc")
### * bbsc

flush(stderr()); flush(stdout())

### Name: bbsc
### Title: Constrained Base-learners for Scalar Covariates
### Aliases: bbsc brandomc bolsc
### Keywords: models

### ** Examples

#### simulate data with functional response and scalar covariate (functional ANOVA)
n <- 60   ## number of cases
Gy <- 27  ## number of observation poionts per response curve 
dat <- list()
dat$t <- (1:Gy-1)^2/(Gy-1)^2
set.seed(123)
dat$z1 <- rep(c(-1, 1), length = n)
dat$z1_fac <- factor(dat$z1, levels = c(-1, 1), labels = c("1", "2"))
# dat$z1 <- runif(n)
# dat$z1 <- dat$z1 - mean(dat$z1)

# mean and standard deviation for the functional response 
mut <- matrix(2*sin(pi*dat$t), ncol = Gy, nrow = n, byrow = TRUE) + 
        outer(dat$z1, dat$t, function(z1, t) z1*cos(pi*t) ) # true linear predictor
sigma <- 0.1

# draw respone y_i(t) ~ N(mu_i(t), sigma)
dat$y <- apply(mut, 2, function(x) rnorm(mean = x, sd = sigma, n = n)) 

## fit function-on-scalar model with a linear effect of z1
m1 <- FDboost(y ~ 1 + bolsc(z1_fac, df = 1), timeformula = ~ bbs(t, df = 6), data = dat)

# look for optimal mSTOP using cvrisk() or validateFDboost()
m1[200] # use 200 boosting iterations 

# plot true and estimated coefficients 
plot(dat$t, 2*sin(pi*dat$t), col = 2, type = "l", main = "intercept")
plot(m1, which = 1, lty = 2, add = TRUE)

plot(dat$t, 1*cos(pi*dat$t), col = 2, type = "l", main = "effect of z1")
lines(dat$t, -1*cos(pi*dat$t), col = 2, type = "l")
plot(m1, which = 2, lty = 2, col = 1, add = TRUE)





cleanEx()
nameEx("bhistx")
### * bhistx

flush(stderr()); flush(stdout())

### Name: bhistx
### Title: Base-learners for Functional Covariates
### Aliases: bhistx
### Keywords: models

### ** Examples

if(require(refund)){
## simulate some data from a historical model
## the interaction effect is in this case not necessary
n <- 100
nygrid <- 35
data1 <- suppressWarnings(pffrSim(scenario = c("int", "ff"), limits = function(s,t){ s <= t }, 
                n = n, nygrid = nygrid))
data1$X1 <- scale(data1$X1, scale = FALSE) ## center functional covariate                  
dataList <- as.list(data1)
dataList$tvals <- attr(data1, "yindex")

## create the hmatrix-object
X1h <- with(dataList, hmatrix(time = rep(tvals, each = n), id = rep(1:n, nygrid), 
                             x = X1, argvals = attr(data1, "xindex"), 
                             timeLab = "tvals", idLab = "wideIndex", 
                             xLab = "myX", argvalsLab = "svals"))
dataList$X1h <- I(X1h)   
dataList$svals <- attr(data1, "xindex")
## add a factor variable 
dataList$zlong <- factor(gl(n = 2, k = n/2, length = n*nygrid), levels = 1:2)  
dataList$z <- factor(gl(n = 2, k = n/2, length = n), levels = 1:2)

## do the model fit with main effect of bhistx() and interaction of bhistx() and bolsc()
mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) + 
               bhistx(x = X1h, df = 5, knots = 5) %X% bolsc(zlong), 
              timeformula = ~ bbs(tvals, knots = 10), data = dataList)
              
## alternative parameterization: interaction of bhistx() and bols()
mod <- FDboost(Y ~ 1 + bhistx(x = X1h, df = 5, knots = 5) %X% bols(zlong), 
              timeformula = ~ bbs(tvals, knots = 10), data = dataList)

}




cleanEx()
nameEx("birthDistribution")
### * birthDistribution

flush(stderr()); flush(stdout())

### Name: birthDistribution
### Title: Densities of live births in Germany
### Aliases: birthDistribution
### Keywords: datasets

### ** Examples

data("birthDistribution", package = "FDboost")

# Plot densities
year_col <- rainbow(70, start = 0.5, end = 1)
year_lty <- c(1, 2, 4, 5)
oldpar <- par(mfrow = c(1, 2))
funplot(1:12, birthDistribution$birth_densities[1:70, ], ylab = "densities", xlab = "month", 
        xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Male")
funplot(1:12, birthDistribution$birth_densities[71:140, ], ylab = "densities", xlab = "month", 
        xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Female")
par(mfrow = c(1, 1))

# fit density-on-scalar model with effects for sex and year
model <- FDboost(birth_densities_clr ~ 1 + bolsc(sex, df = 1) + 
                   bbsc(year, df = 1, differences = 1),
                 # use bbsc() in timeformula to ensure integrate-to-zero constraint
                 timeformula = ~bbsc(month, df = 4, 
                                     # December is followed by January of subsequent year
                                     cyclic = TRUE, 
                                     # knots = {1, ..., 12} with additional boundary knot
                                     # 0 (coinciding with 12) due to cyclic = TRUE
                                     knots = 1:11, boundary.knots = c(0, 12), 
                                     # degree = 1 with these knots yields identity matrix 
                                     # as design matrix
                                     degree = 1),
                 data = birthDistribution, offset = 0, 
                 control = boost_control(mstop = 1000))

# Plotting 'model' yields the clr-transformed effects
par(mfrow = c(1, 3))
plot(model, n1 = 12, n2 = 12)

# Use inverse clr transformation to get effects in Bayes Hilbert space, e.g. for intercept
intercept_clr <- predict(model, which = 1)[1, ]
intercept <- clr(intercept_clr, w = 1, inverse = TRUE)
funplot(1:12, intercept, xlab = "month", xaxp = c(1, 12, 11), pch = 20,
        main = "Intercept", ylab = expression(hat(beta)[0]), id = rep(1, 12))

# Same with predictions
predictions_clr <- predict(model)
predictions <- t(apply(predictions_clr, 1, clr, inverse = TRUE))
pred_ylim <- range(birthDistribution$birth_densities)
par(mfrow = c(1, 2))
funplot(1:12, predictions[1:70, ], ylab = "predictions", xlab = "month", ylim = pred_ylim,
        xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Male")
funplot(1:12, predictions[71:140, ], ylab = "predictions", xlab = "month", ylim = pred_ylim,
        xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Female")
par(oldpar)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("bootstrapCI")
### * bootstrapCI

flush(stderr()); flush(stdout())

### Name: bootstrapCI
### Title: Function to compute bootstrap confidence intervals
### Aliases: bootstrapCI

### ** Examples

if(require(refund)){
#########
# model with linear functional effect, use bsignal()
# Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
set.seed(2121)
data1 <- suppressWarnings(pffrSim(scenario = "ff", n = 40))
data1$X1 <- scale(data1$X1, scale = FALSE)
dat_list <- as.list(data1)
dat_list$t <- attr(data1, "yindex")
dat_list$s <- attr(data1, "xindex")

## model fit by FDboost 
m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 8, df = 3), 
              timeformula = ~ bbs(t, knots = 8), data = dat_list)

}
              

my_inner_fun <- function(object){ 
cvrisk(object, folds = cvLong(id = object$id, weights = 
model.weights(object), B = 2) # 10-fold for inner resampling
) 
}


## We can also use the ... argument to parallelize the applyFolds
## function in the outer resampling 


## Now let's parallelize the outer resampling and use 
## crossvalidation instead of bootstrap for the inner resampling

my_inner_fun <- function(object){ 
cvrisk(object, folds = cvLong(id = object$id, weights = 
model.weights(object), type = "kfold", # use CV
B = 5, # 5-fold for inner resampling
)) # use five cores
}

# use applyFolds for outer function to avoid messing up weights
my_outer_fun <- function(object, fun){
applyFolds(object = object,
folds = cv(rep(1, length(unique(object$id))), 
type = "bootstrap", B = 10), fun = fun) # parallelize on 10 cores
}


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
## include no intercept as all base-learners are centered around 0    

mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
               + bsignal(NIR, nir.lambda, knots = 40, df=4, check.ident = FALSE), 
               timeformula = NULL, data = fuelSubset) 







cleanEx()
nameEx("bsignal")
### * bsignal

flush(stderr()); flush(stdout())

### Name: bsignal
### Title: Base-learners for Functional Covariates
### Aliases: bsignal bconcurrent bhist bfpc
### Keywords: models

### ** Examples

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
data1 <- suppressWarnings(pffrSim(scenario = "ff", n = 40))
data1$X1 <- scale(data1$X1, scale = FALSE)
dat_list <- as.list(data1)
dat_list$t <- attr(data1, "yindex")
dat_list$s <- attr(data1, "xindex")

## model fit by FDboost 
m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 5), 
              timeformula = ~ bbs(t, knots = 5), data = dat_list, 
              control = boost_control(mstop = 21))

## search optimal mSTOP

## model fit by pffr
t <- attr(data1, "yindex")
s <- attr(data1, "xindex")
m1_pffr <- pffr(Y ~ ff(X1, xind = s), yind = t, data = data1)



############################################
# model with functional historical effect, use bhist() 
# Y(t) = f(t)  + \int_0^t X1(s)\beta(s,t)ds + eps
set.seed(2121)
mylimits <- function(s, t){
  (s < t) | (s == t)
}
data2 <- suppressWarnings(pffrSim(scenario = "ff", n = 40, limits = mylimits))
data2$X1 <- scale(data2$X1, scale = FALSE)
dat2_list <- as.list(data2)
dat2_list$t <- attr(data2, "yindex")
dat2_list$s <- attr(data2, "xindex")

## model fit by FDboost 
m2 <- FDboost(Y ~ 1 + bhist(x = X1, s = s, time = t, knots = 5), 
              timeformula = ~ bbs(t, knots = 5), data = dat2_list, 
              control = boost_control(mstop = 40))
              
## search optimal mSTOP
               

## model fit by pffr
t <- attr(data2, "yindex")
s <- attr(data2, "xindex")
m2_pffr <- pffr(Y ~ ff(X1, xind = s, limits = "s<=t"), yind = t, data = data2)



}





cleanEx()
nameEx("clr")
### * clr

flush(stderr()); flush(stdout())

### Name: clr
### Title: Clr and inverse clr transformation
### Aliases: clr

### ** Examples

### Continuous case (T = [0, 1] with Lebesgue measure):
# evaluate density of a Beta distribution on an equidistant grid
g <- seq(from = 0.005, to = 0.995, by = 0.01)
f <- dbeta(g, 2, 5)
# compute clr transformation with distance of two grid points as integration weight
f_clr <- clr(f, w = 0.01)
# visualize result
plot(g, f_clr , type = "l")
abline(h = 0, col = "grey")
# compute inverse clr transformation (w as above)
f_clr_inv <- clr(f_clr, w = 0.01, inverse = TRUE)
# visualize result
plot(g, f, type = "l")
lines(g, f_clr_inv, lty = 2, col = "red")

### Discrete case (T = {1, ..., 12} with sum of dirac measures at t in T):
data("birthDistribution", package = "FDboost")
# fit density-on-scalar model with effects for sex and year
model <- FDboost(birth_densities_clr ~ 1 + bolsc(sex, df = 1) + 
                   bbsc(year, df = 1, differences = 1),
                 # use bbsc() in timeformula to ensure integrate-to-zero constraint
                 timeformula = ~bbsc(month, df = 4, 
                                     # December is followed by January of subsequent year
                                     cyclic = TRUE, 
                                     # knots = {1, ..., 12} with additional boundary knot
                                     # 0 (coinciding with 12) due to cyclic = TRUE
                                     knots = 1:11, boundary.knots = c(0, 12), 
                                     # degree = 1 with these knots yields identity matrix 
                                     # as design matrix
                                     degree = 1),
                 data = birthDistribution, offset = 0, 
                 control = boost_control(mstop = 1000))
# Extract predictions (clr-transformed!) and transform them to Bayes Hilbert space
predictions_clr <- predict(model)
predictions <- t(apply(predictions_clr, 1, clr, inverse = TRUE))




cleanEx()
nameEx("emotion")
### * emotion

flush(stderr()); flush(stdout())

### Name: emotion
### Title: EEG and EMG recordings in a computerised gambling study
### Aliases: emotion
### Keywords: datasets

### ** Examples

data("emotion", package = "FDboost")

# fit function-on-scalar model with random effect and power effect
fos_random_power <- FDboost(EMG ~ 1 + brandomc(subject, df = 2)
                            + bolsc(power, df = 2),
                            timeformula = ~ bbs(t, df = 3),
                            data = emotion)
## Not run: 
##D                             
##D # fit function-on-function model with intercept and historical EEG effect
##D # where limits specifies the used lag between EMG and EEG signal
##D fof_historical <- FDboost(EMG ~ 1 + bhist(EEG, s = s, time = t,
##D                           limits = function(s,t) s < t - 3),
##D                           timeformula = ~ bbs(t, df = 3), data = emotion,
##D                           control = boost_control(mstop = 200))                            
## End(Not run)



cleanEx()
nameEx("factorize")
### * factorize

flush(stderr()); flush(stdout())

### Name: factorize
### Title: Factorize tensor product model
### Aliases: factorize factorise factorize.FDboost

### ** Examples

library(FDboost)

# generate irregular toy data -------------------------------------------------------

n <- 100
m <- 40
# covariates
x <- seq(0,2,len = n)
# time & id
set.seed(90384)
t <- runif(n = n*m, -pi,pi)
id <- sample(1:n, size = n*m, replace = TRUE)

# generate components
fx <- ft <- list()
fx[[1]] <- exp(x)
d <- numeric(2)
d[1] <- sqrt(c(crossprod(fx[[1]])))
fx[[1]] <- fx[[1]] / d[1]
fx[[2]] <- -5*x^2
fx[[2]] <- fx[[2]] - fx[[1]] * c(crossprod(fx[[1]], fx[[2]])) # orthogonalize fx[[2]]
d[2] <- sqrt(c(crossprod(fx[[2]])))
fx[[2]] <- fx[[2]] / d[2]
ft[[1]] <- sin(t)
ft[[2]] <- cos(t)
ft[[1]] <- ft[[1]] / sqrt(sum(ft[[1]]^2))
ft[[2]] <- ft[[2]] / sqrt(sum(ft[[2]]^2))

mu1 <- d[1] * fx[[1]][id] * ft[[1]]
mu2 <- d[2] * fx[[2]][id] * ft[[2]]
# add linear covariate
ft[[3]] <- t^2 * sin(4*t) 
ft[[3]] <- ft[[3]] - ft[[1]] * c(crossprod(ft[[1]], ft[[3]]))
ft[[3]] <- ft[[3]] - ft[[2]] * c(crossprod(ft[[2]], ft[[3]]))
ft[[3]] <- ft[[3]] / sqrt(sum(ft[[3]]^2))
set.seed(9234)
fx[[3]] <- runif(0,3, n = length(x))
fx[[3]] <- fx[[3]] - fx[[1]] * c(crossprod(fx[[1]], fx[[3]]))
fx[[3]] <- fx[[3]] - fx[[2]] * c(crossprod(fx[[2]], fx[[3]]))
d[3] <- sqrt(sum(fx[[3]]^2))
fx[[3]] <- fx[[3]] / d[3]

mu3 <- d[3] * fx[[3]][id] * ft[[3]]

mu <- mu1 + mu2 + mu3
# add some noise
y <- mu + rnorm(length(mu), 0, .01)
# and noise covariate
z <- rnorm(n)

# fit FDboost model -------------------------------------------------------

dat <- list(y = y, x = x, t = t, x_lin = fx[[3]], id = id)
m <- FDboost(y ~ bbs(x, knots = 5, df = 2, differences = 0) + 
               # bbs(z, knots = 2, df = 2, differences = 0) + 
               bols(x_lin, intercept = FALSE, df = 2)
               , ~ bbs(t),
             id = ~ id, 
             offset = 0, #numInt = "Riemann", 
             control = boost_control(nu = 1), 
             data = dat)
MU <- split(mu, id)
PRED <- split(predict(m), id)
Ti <- split(t, id)
t0 <- seq(-pi, pi, length.out = 40)
MU <- do.call(cbind, Map(function(mu, t) approx(t, mu, t0)$y,
          MU, Ti)) 
PRED <- do.call(cbind, Map(function(mu, t) approx(t, mu, t0)$y,
            PRED, Ti))

opar <- par(mfrow = c(2,2))
image(t0, x, MU)
contour(t0, x, MU, add = TRUE)
image(t0, x, PRED)
contour(t0, x, PRED, add = TRUE)
persp(t0, x, MU, zlim = range(c(MU, PRED), na.rm = TRUE))
persp(t0, x, PRED, zlim = range(c(MU, PRED), na.rm = TRUE))
par(opar)

# factorize model ---------------------------------------------------------

fac <- factorize(m)

vi <- as.data.frame(varimp(fac$cov))
# if(require(lattice))
#   barchart(variable ~ reduction, group = blearner, vi, stack = TRUE)

cbind(d^2, sort(vi$reduction, decreasing = TRUE)[1:3])


x_plot <- list(x, x, fx[[3]])

cols <- c("cornflowerblue", "darkseagreen", "darkred")
opar <- par(mfrow = c(3,2))
wch <- c(1,2,10)
for(w in seq_along(wch)) {
  plot.mboost(fac$resp, which = wch[w], col = "darkgrey", ask = FALSE,
       main = names(fac$resp$baselearner[wch[w]]))
  lines(sort(t), ft[[w]][order(t)]*max(d), col = cols[w], lty = 2)
  plot(fac$cov, which = wch[w], 
       main = names(fac$cov$baselearner[wch[w]]))
  points(x_plot[[w]], d[w] * fx[[w]] / max(d), col = cols[w], pch = 3)
}
par(opar)

# re-compose predictions
preds <- lapply(fac, predict)
predf <- rowSums(preds$resp * preds$cov[id, ])
PREDf <- split(predf, id)
PREDf <- do.call(cbind, Map(function(mu, t) approx(t, mu, t0)$y,
                            PREDf, Ti))
opar <- par(mfrow = c(1,2))
image(t0,x, PRED, main = "original prediction")
contour(t0,x, PRED, add = TRUE)
image(t0,x,PREDf, main = "recomposed")
contour(t0,x, PREDf, add = TRUE)
par(opar)

stopifnot(all.equal(PRED, PREDf))

# check out other methods
set.seed(8399)
newdata_resp <- list(t = sort(runif(60, min(t), max(t))))
a <- predict(fac$resp, newdata = newdata_resp, which = 1:5)
plot(newdata_resp$t, a[, 1])
# coef method
cf <- coef(fac$resp, which = 1)


# check factorization on a new dataset ------------------------------------

t_grid <- seq(-pi,pi,len = 30)
x_grid <- seq(0,2,len = 30)
x_lin_grid <- seq(min(dat$x_lin), max(dat$x_lin), len = 30)

# use grid data for factorization
griddata <- expand.grid(
  # time 
  t = t_grid,
  # covariates
  x = x_grid,
  x_lin = 0
)

griddata_lin <- expand.grid(
  t = seq(-pi, pi, len = 30),
  x = 0,
  x_lin = x_lin_grid
)

griddata <- rbind(griddata, griddata_lin)

griddata$id <- as.numeric(factor(paste(griddata$x, griddata$x_lin, sep = ":")))

fac2 <- factorize(m, newdata = griddata)

ratio <- -max(abs(predict(fac$resp, which = 1))) / max(abs(predict(fac2$resp, which = 1)))

opar <- par(mfrow = c(3,2))
wch <- c(1,2,10)
for(w in seq_along(wch)) {
  plot.mboost(fac$resp, which = wch[w], col = "darkgrey", ask = FALSE,
                                main = names(fac$resp$baselearner[wch[w]]))
  
  lines(sort(griddata$t), 
        ratio*predict(fac2$resp, which = wch[w])[order(griddata$t)], 
        col = cols[w], lty = 2)
  plot(fac$cov, which = wch[w], 
       main = names(fac$cov$baselearner[wch[w]]))
  this_x <- fac2$cov$model.frame(which = wch[w])[[1]][[1]]
    lines(sort(this_x), 1/ratio*predict(fac2$cov, which = wch[w])[order(this_x)], 
          col = cols[w], lty = 1)
}
par(opar)

# check predictions
p <- predict(fac2$resp, which = 1)
library(FDboost)

# generate regular toy data --------------------------------------------------

n <- 100
m <- 40
# covariates
x <- seq(0,2,len = n)
# time 
t <- seq(-pi,pi,len = m)
# generate components
fx <- ft <- list()
fx[[1]] <- exp(x)
d <- numeric(2)
d[1] <- sqrt(c(crossprod(fx[[1]])))
fx[[1]] <- fx[[1]] / d[1]
fx[[2]] <- -5*x^2
fx[[2]] <- fx[[2]] - fx[[1]] * c(crossprod(fx[[1]], fx[[2]])) # orthogonalize fx[[2]]
d[2] <- sqrt(c(crossprod(fx[[2]])))
fx[[2]] <- fx[[2]] / d[2]
ft[[1]] <- sin(t)
ft[[2]] <- cos(t)
ft[[1]] <- ft[[1]] / sqrt(sum(ft[[1]]^2))
ft[[2]] <- ft[[2]] / sqrt(sum(ft[[2]]^2))
mu1 <- d[1] * fx[[1]] %*% t(ft[[1]])
mu2 <- d[2] * fx[[2]] %*% t(ft[[2]])
# add linear covariate
ft[[3]] <- t^2 * sin(4*t) 
ft[[3]] <- ft[[3]] - ft[[1]] * c(crossprod(ft[[1]], ft[[3]]))
ft[[3]] <- ft[[3]] - ft[[2]] * c(crossprod(ft[[2]], ft[[3]]))
ft[[3]] <- ft[[3]] / sqrt(sum(ft[[3]]^2))
set.seed(9234)
fx[[3]] <- runif(0,3, n = length(x))
fx[[3]] <- fx[[3]] - fx[[1]] * c(crossprod(fx[[1]], fx[[3]]))
fx[[3]] <- fx[[3]] - fx[[2]] * c(crossprod(fx[[2]], fx[[3]]))
d[3] <- sqrt(sum(fx[[3]]^2))
fx[[3]] <- fx[[3]] / d[3]
mu3 <- d[3] * fx[[3]] %*% t(ft[[3]])

mu <- mu1 + mu2 + mu3
# add some noise
y <- mu + rnorm(length(mu), 0, .01)
# and noise covariate
z <- rnorm(n)

# fit FDboost model -------------------------------------------------------

dat <- list(y = y, x = x, t = t, x_lin = fx[[3]])
m <- FDboost(y ~ bbs(x, knots = 5, df = 2, differences = 0) + 
               # bbs(z, knots = 2, df = 2, differences = 0) + 
               bols(x_lin, intercept = FALSE, df = 2)
               , ~ bbs(t), offset = 0, 
             control = boost_control(nu = 1), 
             data = dat)

opar <- par(mfrow = c(1,2))
image(t, x, t(mu))
contour(t, x, t(mu), add = TRUE)
image(t, x, t(predict(m)))
contour(t, x, t(predict(m)), add = TRUE)
par(opar)

# factorize model ---------------------------------------------------------

fac <- factorize(m)

vi <- as.data.frame(varimp(fac$cov))
# if(require(lattice))
#   barchart(variable ~ reduction, group = blearner, vi, stack = TRUE)

cbind(d^2, vi$reduction[c(1:2, 10)])


x_plot <- list(x, x, fx[[3]])

cols <- c("cornflowerblue", "darkseagreen", "darkred")
opar <- par(mfrow = c(3,2))
wch <- c(1,2,10)
for(w in seq_along(wch)) {
  plot.mboost(fac$resp, which = wch[w], col = "darkgrey", ask = FALSE,
       main = names(fac$resp$baselearner[wch[w]]))
  lines(t, ft[[w]]*max(d), col = cols[w], lty = 2)
  plot(fac$cov, which = wch[w], 
       main = names(fac$cov$baselearner[wch[w]]))
  points(x_plot[[w]], d[w] * fx[[w]] / max(d), col = cols[w], pch = 3)
}
par(opar)

# re-compose prediction
preds <- lapply(fac, predict)
PREDSf <- array(0, dim = c(nrow(preds$resp),nrow(preds$cov)))
for(i in seq_len(ncol(preds$resp)))
  PREDSf <- PREDSf + preds$resp[,i] %*% t(preds$cov[,i])

opar <- par(mfrow = c(1,2))
image(t,x, t(predict(m)), main = "original prediction")
contour(t,x, t(predict(m)), add = TRUE)
image(t,x,PREDSf, main = "recomposed")
contour(t,x, PREDSf, add = TRUE)
par(opar)
# => matches
stopifnot(all.equal(as.numeric(t(predict(m))), as.numeric(PREDSf)))

# check out other methods
set.seed(8399)
newdata_resp <- list(t = sort(runif(60, min(t), max(t))))
a <- predict(fac$resp, newdata = newdata_resp, which = 1:5)
plot(newdata_resp$t, a[, 1])
# coef method
cf <- coef(fac$resp, which = 1)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("fuelSubset")
### * fuelSubset

flush(stderr()); flush(stdout())

### Name: fuelSubset
### Title: Spectral data of fossil fuels
### Aliases: fuelSubset
### Keywords: datasets

### ** Examples


    data("fuelSubset", package = "FDboost")
    
    ## center the functional covariates per observed wavelength
    fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
    fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)

    ## to make mboost::df2lambda() happy (all design matrix entries < 10)
    ## reduce range of argvals to [0,1] to get smaller integration weights
    fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
                                          (max(uvvis.lambda) - min(uvvis.lambda) ))
    fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
                                          (max(nir.lambda) - min(nir.lambda) ))


    ### fit mean regression model with 100 boosting iterations,
    ### step-length 0.1 and
    mod <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4, check.ident=FALSE)
                   + bsignal(NIR, nir.lambda, knots=40, df=4, check.ident=FALSE),
                   timeformula = NULL, data = fuelSubset)
    summary(mod)
    ## plot(mod)



cleanEx()
nameEx("funplot")
### * funplot

flush(stderr()); flush(stdout())

### Name: funplot
### Title: Plot functional data with linear interpolation of missing values
### Aliases: funplot

### ** Examples




cleanEx()
nameEx("grapes-Xc-grapes")
### * grapes-Xc-grapes

flush(stderr()); flush(stdout())

### Name: %Xc%
### Title: Constrained row tensor product
### Aliases: %Xc%

### ** Examples

 
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
                
## check centering around intercept
colMeans(predict(mod1, which = 4))

## check centering around main effects
colMeans(predict(mod1, which = 4)[viscosity$T_A == "low", ])
colMeans(predict(mod1, which = 4)[viscosity$T_A == "high", ])
colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])

## find optimal mstop using cvrsik() or validateFDboost()
## ... 

## look at interaction effect in one plot
# funplot(mod1$yind, predict(mod1, which=4))




cleanEx()
nameEx("hmatrix")
### * hmatrix

flush(stderr()); flush(stdout())

### Name: hmatrix
### Title: A S3 class for univariate functional data on a common grid
### Aliases: hmatrix

### ** Examples

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
getX(myhmatrixSub)
getX(myhmatrix)

# get time
myhmatrix[ , 1] # as column matrix as drop = FALSE
getTime(myhmatrix) # as vector

# get id
myhmatrix[ , 2] # as column matrix as drop = FALSE
getId(myhmatrix) # as vector

# subset hmatrix on the basis of an index, which is defined on the curve level
reweightData(data = list(hmat = myhmatrix), vars = "hmat", index = c(1, 1, 2))
# this keeps only the unique x values in attr(,'x') but multiplies the corresponding
# ids and times in the time id matrix 
# for bhistx baselearner, there may be an additional id variable for the tensor product
newdat <- reweightData(data = list(hmat = myhmatrix, 
  repIDx = rep(seq_len(nrow(attr(myhmatrix,'x'))), length(attr(myhmatrix,"argvals")))), 
  vars = "hmat", index = c(1,1,2), idvars="repIDx")
length(newdat$repIDx) 

## use hmatrix within a data.frame
mydat <- data.frame(I(myhmatrix), z=rnorm(3)[id1])
str(mydat)
str(mydat[id1 %in% c(2, 3), ])
str(myhmatrix[id1 %in% c(2, 3), ])




cleanEx()
nameEx("integrationWeights")
### * integrationWeights

flush(stderr()); flush(stdout())

### Name: integrationWeights
### Title: Functions to compute integration weights
### Aliases: integrationWeights integrationWeightsLeft

### ** Examples

## Example for trapezoidal integration weights
xind0 <- seq(0,1,l = 5)
xind <- c(0, 0.1, 0.3, 0.7, 1)
X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)

# Regualar observation points
integrationWeights(X1, xind0)
# Irregular observation points
integrationWeights(X1, xind)

# with missing value
X1[1,2] <- NA
integrationWeights(X1, xind0)
integrationWeights(X1, xind)

## Example for left integration weights
xind0 <- seq(0,1,l = 5)
xind <- c(0, 0.1, 0.3, 0.7, 1)
X1 <- matrix(xind^2, ncol = length(xind0), nrow = 2)

# Regular observation points
integrationWeightsLeft(X1, xind0, leftWeight = "mean") 
integrationWeightsLeft(X1, xind0, leftWeight = "first") 
integrationWeightsLeft(X1, xind0, leftWeight = "zero")

# Irregular observation points
integrationWeightsLeft(X1, xind, leftWeight = "mean") 
integrationWeightsLeft(X1, xind, leftWeight = "first") 
integrationWeightsLeft(X1, xind, leftWeight = "zero")

# obervation points that do not start with 0
xind2 <- xind + 0.5
integrationWeightsLeft(X1, xind2, leftWeight = "zero")
 



cleanEx()
nameEx("reweightData")
### * reweightData

flush(stderr()); flush(stdout())

### Name: reweightData
### Title: Function to Reweight Data
### Aliases: reweightData

### ** Examples

## load data
data("viscosity", package = "FDboost")
interval <- "101"
end <- which(viscosity$timeAll == as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[ , 1:end])
viscosity$time <- viscosity$timeAll[1:end]

## what does data look like
str(viscosity)

## do some reweighting
# correct weights
str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = c(0, 32, 32, rep(0, 61))))

str(visNew <- reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = c(0, 32, 32, rep(0, 61))))
# check the result
# visNew$vis[1:5, 1:5] ## image(visNew$vis)

# incorrect weights
str(reweightData(viscosity, vars=c("vis", "T_C", "T_A", "rspeed", "mflow"), 
    argvals = "time", weights = sample(1:64, replace = TRUE)), 1)

# supply meaningful index
str(visNew <- reweightData(viscosity, vars = c("vis", "T_C", "T_A", "rspeed", "mflow"), 
              argvals = "time", index = rep(1:32, each = 2)))
# check the result
# visNew$vis[1:5, 1:5]

# errors
if(FALSE){
   reweightData(viscosity, argvals = "")
   reweightData(viscosity, argvals = "covThatDoesntExist", index = rep(1,64))
   }
   



cleanEx()
nameEx("stabsel.FDboost")
### * stabsel.FDboost

flush(stderr()); flush(stdout())

### Name: stabsel.FDboost
### Title: Stability Selection
### Aliases: stabsel.FDboost

### ** Examples

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


## create folds for stability selection  
## only 5 folds for a fast example, usually use 50 folds 
set.seed(1911)
folds <- cvLong(modAll$id, weights = rep(1, l = length(modAll$id)), 
                type = "subsampling", B = 5) 
    




cleanEx()
nameEx("subset_hmatrix")
### * subset_hmatrix

flush(stderr()); flush(stdout())

### Name: subset_hmatrix
### Title: Subsets hmatrix according to an index
### Aliases: subset_hmatrix

### ** Examples

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
resMat <- subset_hmatrix(hmat, index = index1, compress = FALSE)
try(resMat2 <- subset_hmatrix(resMat, index = index2))




cleanEx()
nameEx("truncateTime")
### * truncateTime

flush(stderr()); flush(stdout())

### Name: truncateTime
### Title: Function to truncate time in functional data
### Aliases: truncateTime

### ** Examples

if(require(fda)){
  dat <- fda::growth
  dat$hgtm <- t(dat$hgtm[,1:10])
  dat$hgtf <- t(dat$hgtf[,1:10])
  
  ## only use time-points 1:16 of variable age
  datTr <- truncateTime(funVar=c("hgtm","hgtf"), time="age", newtime=1:16, data=dat)
  
}



cleanEx()
nameEx("update.FDboost")
### * update.FDboost

flush(stderr()); flush(stdout())

### Name: update.FDboost
### Title: Function to update FDboost objects
### Aliases: update.FDboost

### ** Examples

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
               
# update nu
mod2 <- update(mod1, control=boost_control(nu = 1)) # mstop will stay the same
# update mstop
mod3 <- update(mod2, control=boost_control(mstop = 100)) # nu=1 does not get changed
mod4 <- update(mod1, formula = vis ~ 1 + bolsc(T_C, df = 2)) # drop one term



cleanEx()
nameEx("validateFDboost")
### * validateFDboost

flush(stderr()); flush(stdout())

### Name: validateFDboost
### Title: Cross-Validation and Bootstrapping over Curves
### Aliases: validateFDboost

### ** Examples





cleanEx()
nameEx("viscosity")
### * viscosity

flush(stderr()); flush(stdout())

### Name: viscosity
### Title: Viscosity of resin over time
### Aliases: viscosity
### Keywords: datasets

### ** Examples


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
 summary(mod)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
