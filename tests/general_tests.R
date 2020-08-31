

library(FDboost)
library(gamboostLSS)

# print(sessionInfo())

# simulated data ----------------------------------------------------------


if(require(refund)){
  
  ## simulate a small data set 
  print("simulate data")
  set.seed(230)
  pffr_data <- pffrSim(n = 25, nxgrid = 21, nygrid = 19)
  pffr_data$X1 <- scale(pffr_data$X1, scale = FALSE)
  
  dat <- as.list(pffr_data)
  dat$tvals <- attr(pffr_data, "yindex")
  dat$svals <- attr(pffr_data, "xindex")
  
  dat$Y_scalar <- dat$Y[ , 10]
  
  dat$Y_long <- c(dat$Y)
  dat$tvals_long <- rep(dat$tvals, each = nrow(dat$Y))
  dat$id_long <- rep(1:nrow(dat$Y), ncol(dat$Y))
  
  # second functional covariate
  dat$s2 <- seq(0, 1, l = 15)
  dat$X2 <- I(matrix(rnorm(25 * 15), nrow = 25))
  dat$X2 <- scale(dat$X2, scale = FALSE)
  

  # model fit ---------------------------------------------------------------

  print("model fit")
  
  ## response matrix for response observed on one common grid 
  m <- FDboost(Y ~ 1 + bhist(X1, svals, tvals, knots = 10, df = 6) 
               + bsignal(X1, svals, knots = 6, df = 3)
               + bbsc(xsmoo, knots = 6, df = 3) 
               + bolsc(xte1, df = 3)
               + brandomc(xte2, df = 3), 
               timeformula = ~ bbs(tvals, knots = 9, df = 2, differences = 1), 
               control = boost_control(mstop = 10), data = dat)
  
  ## response in long format
  ml <- FDboost(Y_long ~ 1 + bhist(X1, svals, tvals_long, knots = 6, df = 12) 
                + bsignal(X1, svals, knots = 6, df = 4)
                + bbsc(xsmoo, knots = 6, df = 4) 
                + bolsc(xte1, df = 4)
                + brandomc(xte2, df = 4), 
                timeformula = ~ bbs(tvals_long, knots = 8, df = 3, differences = 1), 
                id = ~ id_long, 
                offset_control = o_control(k_min = 10), 
                control = boost_control(mstop = 10), data = dat)
  
  ## scalar response 
  ms <- FDboost(Y_scalar ~ 1 + bsignal(X1, svals, knots = 6, df = 2)
                + bbs(xsmoo, knots = 6, df = 2, differences = 1) 
                + bols(xte1, df = 2) 
                + bols(xte2, df = 2)
                + bols(xfactor, df = 2), 
                timeformula = NULL, 
                control = boost_control(mstop = 50), data = dat)
  
  ## scalar response and interaction effect between two functional variables
  ms_funint <- FDboost(Y_scalar ~ 1 + 
                         bsignal(X1, svals, knots = 9, df = 3) %X% bsignal(X2, s2, knots = 9, df = 3), 
                timeformula = NULL, 
                control = boost_control(mstop = 50), data = dat)
  
  ## GAMLSS with functional response 
  mlss <- FDboostLSS(Y ~ 1 + bsignal(X1, svals, knots = 6, df = 3)               
                     + bbsc(xsmoo, knots = 6, df = 3) 
                     + bolsc(xte1, df = 3), 
                     timeformula = ~ bbs(tvals, knots = 9, df = 3, differences = 1), 
                     control = boost_control(mstop = 20), data = dat, 
                     method = "noncyclic")
  
  ## GAMLSS with scalar response 
  mslss <- FDboostLSS(Y_scalar ~ 1 + bsignal(X1, svals, knots = 6, df = 3)
                      + bbs(xsmoo, knots = 6, df = 3, differences = 1), 
                      timeformula = NULL, 
                      control = boost_control(mstop = 50), data = dat, 
                      method = "noncyclic")
  
  ## response matrix with factor + continuous time variable
  
  # linear array model implemented only for matrices
  # => tvals and factor for dimensions have to be flattened
  dat2D <- with(dat, expand.grid(
    tvals = tvals, 
    xfactor = factor(2:3)
  ))
  dat2D <- as.list(dat2D[order(dat2D$xfactor), ])
  dat2D$Y <- cbind(
    dat$Y[dat$xfactor == "2", ],
    dat$Y[dat$xfactor == "3", ]
    )
  dat2D$xsmoo <- dat$xsmoo[dat$xfactor == "2"] 
    
  m2D <- FDboost(Y ~ bbsc(xsmoo), 
                 timeformula = ~ bols(xfactor) %X% bbs(tvals),
                 control = boost_control(mstop = 20), 
                 data = dat2D)
  
  
  # test some methods and utility functions  --------------------------------
  
  ## test plot()
  print("plot effects")
  par(mfrow = c(1,1))
  plot(m, ask = FALSE)
  plot(ml, ask = FALSE)
  plot(ms, ask = FALSE)
  plot(mlss$mu, ask = FALSE) 
  plot(mlss$sigma, ask = FALSE)
  
  ## test applyFolds()
  print("run applyFolds")
  set.seed(123)
  applyFolds(m, folds = cv(rep(1, length(unique(m$id))), B = 2), grid = 0:5)
  #applyFolds(ml, folds = cv(rep(1, length(unique(ml$id))), B = 2), grid = 0:5)
  #applyFolds(ms, folds = cv(rep(1, length(unique(ms$id))), B = 2), grid = 0:5)

  ## test cvrisk()
  print("run cvrisk")
  set.seed(123)
  cvrisk(m, folds = cvLong(id = m$id, weights = model.weights(m), B = 2), grid = 0:5)
  cvrisk(ml, folds = cvLong(id = ml$id, weights = model.weights(ml), B = 2), grid = 0:5)
  cvrisk(ms, folds = cvLong(id = ms$id, weights = model.weights(ms), B = 2), grid = 0:5)
  cvrisk(ms_funint, folds = cvLong(id = ms$id, weights = model.weights(ms), B = 2), grid = 0:5)
  
  cvrisk(mlss, folds = cv(model.weights(mlss[[1]]), B = 2),
         grid = 1:5, trace = FALSE)
  cvrisk(mslss, folds = cv(model.weights(mslss[[1]]), B = 2),
         grid = 1:5, trace = FALSE)
  
  cvrisk(m2D, folds = cv(model.weights(m2D), B = 2),
         grid = 1:5)
  
  ## test stabsel (use very small number of folds, B = 10, to speed up testing)
  print("run stabsel")
  stabsel(m, cutoff=0.8, PFER = 0.1*length(m$baselearner), sampling.type = "SS", eval = TRUE, B = 3)
  stabsel(m, cutoff=0.8, PFER = 0.1*length(m$baselearner), sampling.type = "SS", eval = TRUE, B = 3, 
          refitSmoothOffset = FALSE)
  
  ## FIXME: this stabsel() should also work with refitSmoothOffset = TRUE 
  stabsel(ml, cutoff=0.8, PFER = 0.1*length(ml$baselearner), sampling.type = "SS", eval = TRUE, B = 3, 
          refitSmoothOffset = FALSE)
  stabsel(ms, cutoff=0.8, PFER = 0.1*length(ms$baselearner), sampling.type = "SS", eval = TRUE, B = 3)
  ## FIXME: this is broken again 
  ## fixed in gamboostLSS package on github with commit 4989474 
  ##stabsel(mlss, cutoff=0.8, PFER = 0.1*length(mlss$mu$baselearner), sampling.type = "SS", eval = TRUE, B = 3)
  ##stabsel(mslss, cutoff=0.8, PFER = 0.1*length(mslss$mu$baselearner), sampling.type = "SS", eval = TRUE, B = 3)
  
  
  ## test predict with newdata
  print("predict with new data")
  pred <- predict(m, newdata = dat)
  ## for this predict() with newdata, you need a data.frame where the irregular time of y fits with X
  ## pred <- predict(ml, newdata = dat)
  pred <- predict(ms, newdata = dat)
  pred <- predict(ms_funint, newdata = dat)
  pred <- predict(mlss, newdata = dat)
  pred <- predict(mslss, newdata = dat)
  
}



# sof: fuel data ----------------------------------------------------------

print("run checks with fuel data")

## prediction with functional variable as numeric matrix, see Issue #17
data(fuelSubset)
fuel <- fuelSubset[c('heatan', 'h2o', 'UVVIS', 'uvvis.lambda')]
str(fuel$UVVIS) # numeric matrix

sof <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 20, df = 4) + 
                 bbs(h2o, df = 4),
               timeformula = ~bols(1), data = fuel)

# Predict with newdata
pred <- predict(sof, newdata = fuel)

# with interaction term
sof_int <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 9, df = 9) + 
                 bsignal(NIR, nir.lambda, knots = 9, df = 9) +
                 bsignal(UVVIS, uvvis.lambda, knots = 9, df = 3) %X% bsignal(NIR, nir.lambda, knots = 9, df = 3),
               timeformula = ~bols(1), data = fuelSubset)

# Predict with newdata
pred <- predict(sof_int, newdata = fuelSubset)



# fof: fuel data -----------------------------------------------------------

## model does not make sense, but is good for checking

# function-on-function with bsignal
fof <- FDboost(UVVIS ~ bsignal(NIR, nir.lambda, knots = 9, df = 9),
               timeformula = ~ bbs(uvvis.lambda), data = fuelSubset)

pred <- predict(fof, newdata = fuelSubset)




