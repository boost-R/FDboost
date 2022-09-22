
#' Factorize tensor product model
#' 
#' Factorize an FDboost tensor product model into the response and covariate parts 
#' \deqn{h_j(x, t) = \sum_{k} v_j^{(k)}(t) h_j^{(k)}(x), j = 1, ..., J,}
#' for effect visualization as proposed in Stoecker, Steyer and Greven (2022).
#'
#' @param x a model object of class FDboost.
#' @param ... other arguments passed to methods.
#'
#' @details The mboost infrastructure is used for handling the orthogonal response 
#' directions \eqn{v_j^{(k)}(t)} in one \code{mboost}-object 
#' (with \eqn{k} running over iteration indices) and the effects into the respective 
#' directions \eqn{h_j^{(k)}(t)} in another \code{mboost}-object, 
#' both of subclass \code{FDboost_fac}. 
#' The number of boosting iterations of \code{FDboost_fac}-objects cannot be 
#' further increased as in regular \code{mboost}-objects.
#'
#' @return a list of two mboost models of class \code{FDboost_fac} containing basis functions
#' for response and covariates, respectively, as base-learners.
#' @export
#' 
#' @name factorize
#' @aliases factorise factorize.FDboost
#' @importFrom MASS ginv
#' @importFrom Matrix rankMatrix
#' @seealso [FDboost_fac-class]
#'
#' @references 
#' Stoecker, A., Steyer L. and Greven, S. (2022):
#' Functional additive models on manifolds of planar shapes and forms
#' <arXiv:2109.02624>
#' 
#'
#' @example tests/factorize_test_irregular.R 
#' @example tests/factorize_test_regular.R
#'
factorize <- factorise <- function(x, ...) {
  UseMethod("factorize")
}

#' @param newdata new data the factorization is based on. 
#' By default (\code{NULL}), the factorization is carried out on the data used for fitting.
#' @param newweights vector of the length of the data or length one, 
#' containing new weights used for factorization.
#' @param blwise logical, should the factorization be carried out base-learner-wise (\code{TRUE}, default)
#' or for the whole model simultaneously.
#'
#' @method factorize FDboost
#' @return A factorized model
#' @export 
#' @rdname factorize
factorize.FDboost <- function(x, newdata = NULL, newweights = 1, blwise = TRUE, ...) {

  FDboost_regular <- !inherits(x, c("FDboostScalar", "FDboostLong"))

  nd <- !is.null(newdata)

  # built subdata
  dat <- list()
  dat$cov <- if(!nd) x$data
                else newdata[names(x$data)]
  dat$cov[[x$yname]] <- rep(1, min(lengths(dat$cov)))
  if(is.list(x$yind)) {
    dat$resp <- if(!nd) x$yind else newdata[names(x$yind)]
  } else {
    dat$resp <- if(!nd) setNames(list(x$yind),
                                 attr(x$yind, "nameyind")) else
                                   newdata[attr(x$yind, "nameyind")]
  }
  dat$resp <- as.data.frame(dat$resp)
  dat$resp[[x$yname]] <- 1

  # extract formulae
  formulae <- list()
  formulae$cov <- as.formula(x$formulaFDboost)
  formulae$resp <- as.formula(paste(x$yname, x$timeformula))

  # set up component models
  mod <- list()
  # standard mboost model for covariates
  mod$cov <- mboost(formulae$cov,
                    data = dat$cov,
                    offset = 0,
                    control = boost_control(mstop = 0, nu = 1))
  # artificial FDboost intercept model for response
  mod$resp <- mboost(formulae$resp,
                    data = dat$resp,
                    offset = if(FDboost_regular)
                      matrix(x$offset, nrow = x$ydim[1])[1,] else
                        x$offset,
                    control = boost_control(mstop = 0, nu = 1))
  # copy essential parts from base model to response
  which_vars <- c("yname", "ydim", "predictOffset", "withIntercept",
                  "callEval", "timeformula", "formulaFDboost",
                  "formulaMboost", "family", "(weights)", "id")
  cls <- class(mod$resp)
  mod$resp[which_vars] <- unclass(x)[which_vars]
  if(FDboost_regular) mod$resp$ydim <- c(1, x$ydim[2])
  mod$resp$yind <- range(x$yind)
  attr(mod$resp$yind, "nameyind") <- attr(x$yind, "nameyind")
  if(FDboost_regular)
    class(mod$resp) <- c("FDboostLong", class(x)) else
      class(mod$resp) <- class(x)

  if(nd) {
    if(length(newweights)==1)
      mod$resp[["(weights)"]] <- rep(newweights, length(dat$resp[[x$yname]])) else {
      stopifnot(length(newweights) == length(dat$resp[[x$yname]]))
      mod$resp[["(weights)"]] <- newweights
      }
    mod$resp$id <- newdata[[attr(mod$resp$id, "nameid")]]
  }

  # set to FDboost_fac class
  for(i in names(mod))
    class(mod[[i]]) <- c("FDboost_fac", class(mod[[i]]))

  # get coefficients (only of selected learners)
  bl_selected <- x$which(usedonly = TRUE)
  cf <- coef(x, raw = TRUE, which = bl_selected)

  # extract design matrices

  X <- list(
    cov = extract(mod$cov, what = "design", which = bl_selected),
    resp = extract(mod$resp, what = "design", which = 1)
  )
  index <- list(
    cov = extract(mod$cov, what = "index", which = bl_selected),
    resp = extract(mod$resp, what = "index", which = 1)
  )

  wghts <-  mod$resp$`(weights)`

  if(is.null(wghts)) {
    wghts <- list(cov = 1, resp = 1)
  } else {
    if(FDboost_regular) {

        dim(wghts) <- x$ydim
        wghts <- list(
          cov = rowMeans(wghts),
          resp = wghts[1, ]
        )
    } else {
      wghts <- list(cov = as.vector(tapply(wghts, mod$resp$id, mean)))
      wghts$resp <- mod$resp[["(weights)"]] / wghts$cov[mod$resp$id]
    }
  }

  wghts <- Map(function(w, idx) {
    lapply(idx, function(i) {
      if(is.null(i)) w else
        c(tapply(w, i, sum))
    })
  }, wghts, index)

  # multiply sqrt(weights) to X to take them into account
  X <- Map(function(x,w) {
    Map(function(.x, .w) sqrt(.w) * .x, x,w)
  }, X, wghts)
  # NOTE: X is now sqrt(w) * X !

  # do QR decomposition to achieve orthonormal basis representation
  QR <- lapply(X, lapply, qr)
  ## extract Q as orthonormal version of X
  # Q <- lapply(QR, lapply, qr.Q) # not necessary

  # transform cf accordingly
  R <- lapply(QR, lapply, function(x) {
    if(inherits(x, "qr"))
       qr.R(x)[, order(x$pivot)] else
         qrR(x, backPermute = TRUE) })

  cf <- Map(matrix, cf, nrow = lapply(X$cov, ncol), byrow = !FDboost_regular)
  cf <- Map(function(r1, o) r1 %*% tcrossprod(o, R$resp[[1]]), R$cov, cf)

  # perform SVD on cf
  if(blwise) {
    SVD <- lapply(cf, svd)
    Ud <- lapply(SVD, function(x) sweep(x$u, 2, x$d, "*"))
    d2 <- list(cov = lapply(SVD, function(x) (x$d)^2))
    V <- lapply(SVD, `[[`, "v")
    rm(SVD)
  } else {
      cf <- do.call(rbind, cf)
      SVD <- svd(cf)
      cfidx <- relist(seq_len(nrow(cf)),
                      lapply(X$cov, function(x) numeric(ncol(x))))
      Ud <- lapply(cfidx, function(idx)
        sweep(SVD$u[idx, , drop = FALSE], 2, SVD$d, "*"))
      d2 <- list(
        cov = lapply(Ud, function(ud) colSums(ud^2)),
        resp = SVD$d^2
        )
      V <- list(model = SVD$v)
      rm(SVD)
    }

  # compute new coefs
  d_max <- sqrt(max(unlist(d2)))
  if(d_max == 0) d_max <- 1
  cf <- list()
  my_solve <- function(a, b) {
    ret <- try(solve(a, b), silent = TRUE)
    if(inherits(ret, "try-error")) {
      ret <- ginv(a) %*% b
    }
    ret
  }
  cf$cov <- Map(function(R, du) {
    as.matrix(my_solve(R, du)) / d_max
    }, R$cov, Ud)
  cf$resp <- setNames(
    lapply(V, my_solve, a = R$resp[[1]] / d_max),
    nm = if(blwise)
      paste0(names(X$resp)[1], " [", names(X$cov), "]") else
        names(X$resp)[1]
      )
  .no_mat <- which(!sapply(cf$resp, is.matrix))
  cf$resp[.no_mat] <-
    lapply(cf$resp[.no_mat], as.matrix)
  # drop dimension discrepancies
  if(length(cf$cov) == length(cf$resp)) {
    for(bl in seq_along(cf$cov)) {
      nc <- min(NCOL(cf$cov[[bl]]), NCOL(cf$resp[[bl]]))
      cf$cov[[bl]] <- cf$cov[[bl]][, 1:nc, drop = FALSE]
      cf$resp[[bl]] <- cf$resp[[bl]][, 1:nc, drop = FALSE]
    }
  }
  for(bl in seq_along(cf$cov)) {
    d2$cov[[bl]] <- head(d2$cov[[bl]], NCOL(cf$cov[[bl]]))
  }

  # decomposition complete - now prepare output ---------------

  # get model environments
  e <- lapply(mod, function(m) environment(m$predict))

  # clone and equip baselearners
  bl_dims <- lapply(cf, sapply, NCOL)
  # vector for cloning bls
  bl_mltpl <- list(
   cov =  rep(seq_along(bl_dims$cov), bl_dims$cov),
   resp = rep(1, sum(bl_dims$resp))
  )
  bl_names <- Map(function(.cf, .bl_dims)
    unlist(Map(function(name, len) paste0(name, " [", seq_len(len), "]"),
                         names(.cf), .bl_dims), use.names = FALSE),
    cf, bl_dims)
  # order of newly generated bls with respect to their variance
  d2l <- lapply(d2, unlist)
  bl_order <- Map(function(bmlt, d2) order(bmlt)[order(d2, decreasing = TRUE)],
                  bl_mltpl, d2l)

  for(i in names(mod)) {
    this_select <- if(i=="cov") bl_selected else 1
    mod[[i]]$baselearner <- e[[i]]$blg <- setNames(
      e[[i]]$blg[this_select][bl_mltpl[[i]]], bl_names[[i]])
    mod[[i]]$basemodel <- e[[i]]$bl <- setNames(
      e[[i]]$bl[this_select][bl_mltpl[[i]]], bl_names[[i]])
    e[[i]]$bnames <- bl_names[[i]]
    # fill in coefs with bl order decreasing with explained variance
    e[[i]]$xselect <- bl_order[[i]]
    e[[i]]$ens <- unlist(lapply(cf[[i]], asplit, 2), recursive = FALSE)
    e[[i]]$ens <- Map( function(x, cls) {
      bm <- list(model = x)
      class(bm) <- gsub("bl", "bm", cls)
      bm
    },
    x = e[[i]]$ens[bl_order[[i]]],
    cls = lapply(mod[[i]]$basemodel, class)[bl_order[[i]]])
    # add risk
    this_d2l <- d2l[[i]]
    if(is.null(this_d2l))
      this_d2l <- d2l[[1]]
    e[[i]]$mrisk <- sum(this_d2l) -
      cumsum(c(0,sort(this_d2l, decreasing = TRUE)))
    # engage full number of components
    mod[[i]]$subset(sum(this_d2l>0))
  }

  # return factor models
  mod
}


# define class and methods ----------------------------------------------------

#' @importFrom methods setOldClass
#' @exportClass FDboost_fac

setOldClass("FDboost_fac")

#' `FDboost_fac` S3 class for factorized FDboost model components
#'
#' @description Model factorization with `factorize()` decomposes an
#' `FDboost` model into two objects of class `FDboost_fac` - one for the
#' response and one for the covariate predictor. The first is essentially
#' an `FDboost` object and the second an `mboost` object, however,
#' in a 'read-only' mode and slightly adjusted methods (method defaults).
#'
#' @name FDboost_fac-class
#' @seealso [factorize(), factorize.FDboost()]
NULL



#' Prediction and plotting for factorized FDboost model components
#'
#' @param object,x a model-factor given as a \code{FDboost_fac} object
#' @param newdata optionally, a data frame or list
#' in which to look for variables with which to predict.
#' See \code{\link{predict.mboost}}.
#' @param which a subset of base-learner components to take into
#' account for computing predictions or coefficients. Different
#' components are never aggregated to a joint prediction, but always
#' returned as a matrix or list. Select the k-th component
#' by name in the format \code{bl(x, ...)[k]} or all components of a base-learner
#' by dropping the index or all base-learners of a variable by using
#' the variable name.
#' @param main the plot title. By default, base-learner names are used with 
#' component numbers \code{[k]}. 
#' @param ... additional arguments passed to underlying methods.
#'
#' @method predict FDboost_fac
#'
#' @export
#' @name predict.FDboost_fac
#' @aliases plot.FDboost_fac
#' @return A matrix of predictions (for predict method) or no 
#' return value (plot method)
#' @seealso [factorize(), factorize.FDboost()]
#'
predict.FDboost_fac <- function(object, newdata = NULL, which = NULL, ...) {
  w <- object$which(which)
  if(any(is.na(w)))
    stop("Don't know 'which' base-learner is meant.")
  names(w) <- names(object$baselearner)[w]
  drop(sapply(w,
         function(x) predict.mboost(which = x,
                                    object = object,
                                    newdata = newdata,
                                    aggregate = "sum", ...)))
}

#' @method plot FDboost_fac
#' @rdname predict.FDboost_fac
plot.FDboost_fac <- function(x, which = NULL, main = NULL, ...) {
  w <- x$which(which, usedonly = TRUE)
  if(any(is.na(w)))
    stop(paste("Don't know which base-learner is meant by:",
         which[which.min(is.na(w))]))
  if(is.null(main))
    main <- names(x$baselearner)[w]
  for(i in seq_along(w))
    plot.mboost(x, which = w[i], main = main[i], ...)
}