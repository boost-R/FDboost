

library(FDboost)

################################################################
######### simulate some data 

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
  dat$X2 <- matrix(rnorm(25 * 15), nrow = 25)
  dat$X2 <- scale(dat$X2, scale = FALSE)
  
  
  ################################################################
  ######### model fit 
  print("model fit")
  
  ## response matrix for response observed on one common grid 
  m <- FDboost(Y ~ 1 + bhist(X1, svals, tvals, knots = 6, df = 12) 
               + bsignal(X1, svals, knots = 6, df = 4)
               + bbsc(xsmoo, knots = 6, df = 4) 
               + bolsc(xte1, df = 4)
               + brandomc(xte2, df = 4), 
               timeformula = ~ bbs(tvals, knots = 9, df = 3, differences = 1), 
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
                + bols(xte2, df = 2), 
                timeformula = NULL, 
                control = boost_control(mstop = 50), data = dat)
  
  ## scalar response and interaction effect between two functional variables
  ms_funint <- FDboost(Y_scalar ~ 1 + 
                         bsignal(X1, svals, knots = 6, df = 3) %X% bsignal(X2, s2, knots = 6, df = 3), 
                timeformula = NULL, 
                control = boost_control(mstop = 50), data = dat)
  
  ## GAMLSS with functional response 
  mlss <- FDboostLSS(Y ~ 1 + bsignal(X1, svals, knots = 6, df = 4)               
                     + bbsc(xsmoo, knots = 6, df = 4) 
                     + bolsc(xte1, df = 4), 
                     timeformula = ~ bbs(tvals, knots = 9, df = 3, differences = 1), 
                     control = boost_control(mstop = 10), data = dat, 
                     method = "noncyclic")
  
  ## GAMLSS with scalar response 
  mslss <- FDboostLSS(Y_scalar ~ 1 + bsignal(X1, svals, knots = 6, df = 4)
                      + bbs(xsmoo, knots = 6, df = 4, differences = 1), 
                      timeformula = NULL, 
                      control = boost_control(mstop = 50), data = dat, 
                      method = "noncyclic")
  
  
  ################################################################
  ######### test some methods and utility functions 
  
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
  applyFolds(ml, folds = cv(rep(1, length(unique(ml$id))), B = 2), grid = 0:5)
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
  
  ## test stabsel (use very small number of folds B = 10 to speed up testing)
  print("run stabsel")
  stabsel(m, cutoff=0.8, PFER = 0.1*length(m$baselearner), sampling.type = "SS", eval = TRUE, B = 3)
  stabsel(ml, cutoff=0.8, PFER = 0.1*length(ml$baselearner), sampling.type = "SS", eval = TRUE, B = 3)
  stabsel(ms, cutoff=0.8, PFER = 0.1*length(ms$baselearner), sampling.type = "SS", eval = TRUE, B = 3)
  ## fixed in gamboostLSS package on github with commit 4989474 
  #try(stabsel(mlss, cutoff=0.8, PFER = 0.1*length(mlss$mu$baselearner), sampling.type = "SS", eval = TRUE, B = 3))
  #try(stabsel(mslss, cutoff=0.8, PFER = 0.1*length(mslss$mu$baselearner), sampling.type = "SS", eval = TRUE, B = 3))
  
  
  ## test predict with newdata
  print("predict with new data")
  pred <- predict(m, newdata = dat)
  ## pred <- predict(ml, newdata = dat) ## you need a data.frame where the irregular time of y fits with X
  pred <- predict(ms, newdata = dat)
  ## FIXME
  ## pred <- predict(ms_funint, newdata = dat)
  ## pred <- predict(mlss, newdata = dat)
  ## pred <- predict(mslss, newdata = dat)
  
}

