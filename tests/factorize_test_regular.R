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
for(w in 1:length(wch)) {
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
for(i in 1:ncol(preds$resp))
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

