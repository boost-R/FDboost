#' Function to compute bootstrap confidence intervals
#'
#' @param object a fitted model object of class \code{FDboost}, 
#' for which the confidence intervals should be computed.
#' @param which a subset of base-learners to take into account for 
#' computing confidence intervals. 
#' @param resampling_fun_outer function for the outer resampling procedure.
#' \code{resampling_fun_outer} must be a function with arguments \code{object}
#' and \code{fun}, where \code{object} corresponds to the fitted 
#' \code{FDboost} object and \code{fun} is passed to the \code{fun}
#' argument of the resampling function (see examples).
#' If \code{NULL}, \code{\link{applyFolds}} is used with 100-fold boostrap.
#' Although the function can be defined very flexible, it is recommended 
#' to use \code{applyFolds} and, in particular, not \code{cvrisk}, 
#' as in this case, weights of the inner and outer 
#' fold will interact, probably causing the inner 
#' resampling to crash. For bootstrapped confidence intervals
#' the outer function should usually be a bootstrap type of resampling.
#' @param resampling_fun_inner function for the inner resampling procudure,
#' which determines the optimal stopping iteration in each fold of the
#' outer resampling procedure. Should be a function with one argument
#' \code{object} for the fitted \code{FDboost} object. 
#' If \code{NULL}, \code{cvrisk} is used with 25-fold bootstrap.
#' @param B_outer Number of resampling folds in the outer loop.
#' Argument is overwritten, when a custom \code{resampling_fun_outer}
#' is supplied.
#' @param B_inner Number of resampling folds in the inner loop.
#' Argument is overwritten, when a custom \code{resampling_fun_inner}
#' is supplied.
#' @param levels the confidence levels required. If NULL, the 
#' raw results are returned.
#' @param plot logical; whether or not to plot the results.
#'
#' @author David Ruegamer, Sarah Brockhaus
#' 
#' @note Note that parallelization can be achieved by defining
#' the \code{resampling_fun_outer} or \code{_inner} accordingly.
#' See e.g. \code{\link{cvrisk}} on how to parallelize resampling
#' functions or the examples below. Also note that by defining
#' a custum inner or outer resampling function the respective
#' argument \code{B_inner} or \code{B_outer} is ignored.
#' For models with complex baselearners, e.g. created by combining
#' several baselearners with the Kronecker or row-wise tensor product,
#' it is also recommended to use \code{levels = NULL} in order to
#' let the function return the raw results and then manually compute
#' confidence intervals.
#' If a baselearner is not selected in any fold, the function
#' treats its effect as constant zero.
#' 
#'  
#' @return a list containing the elements \code{raw_results}, the 
#' \code{quantiles} and \code{mstops}. 
#' In both list elements, each baselearner
#' selected with \code{which} in turn corresponds to a list
#' element. The quantiles are given as vector, matrix or list of
#' matrices depending on the nature of the effect. In case of functional
#' effects the list element in\code{quantiles} is a \code{length(levels)} times
#' \code{length(effect)} matrix, i.e. the rows correspond to the quantiles.
#' In case of coefficient surfaces, \code{quantiles} comprises a list of matrices,
#' where each list element corresponds to a quantile.
#' 
#' @examples 
#' if(require(refund)){
#' #########
#' # model with linear functional effect, use bsignal()
#' # Y(t) = f(t) + \int X1(s)\beta(s,t)ds + eps
#' set.seed(2121)
#' data1 <- pffrSim(scenario = "ff", n = 40)
#' data1$X1 <- scale(data1$X1, scale = FALSE)
#' dat_list <- as.list(data1)
#' dat_list$t <- attr(data1, "yindex")
#' dat_list$s <- attr(data1, "xindex")
#' 
#' ## model fit by FDboost 
#' m1 <- FDboost(Y ~ 1 + bsignal(x = X1, s = s, knots = 5), 
#'               timeformula = ~ bbs(t, knots = 5), data = dat_list, 
#'               control = boost_control(mstop = 21))
#'
#'}
#'               
#' \dontrun{             
#' # a short example with not so meaningful number of folds
#' bootCIs <- bootstrapCI(m1, B_inner = 3, B_outer = 5)               
#' }
#' 
#' ## now speed things up by defining the inner resampling
#' ## function with parallelization
#' 
#' my_inner_fun <- function(object)
#' { 
#' cvrisk(object, folds = cvLong(id = object$id, weights = 
#' model.weights(object), 
#' B = 10 # 10-fold for inner resampling
#' ), mc.cores = 10) # use ten cores
#' }
#' 
#' \dontrun{
#' bootCIs <- bootstrapCI(m1, resampling_fun_inner = my_inner_fun)
#' }
#' 
#' ## Now let's parallelize the outer resampling and use 
#' ## crossvalidation instead of bootstrap for the inner resampling
#' 
#' my_inner_fun <- function(object)
#' { 
#' cvrisk(object, folds = cvLong(id = object$id, weights = 
#' model.weights(object), type = "kfold", # use CV
#' B = 10 # 10-fold for inner resampling
#' )) # use ten cores
#' }
#' 
#' # use applyFolds for outer function to avoid
#' # mess weights mess up
#' my_outer_fun <- function(object, fun)
#' {
#' applyFolds(object = object,
#' folds = cv(rep(1, length(unique(object$id))), 
#' type = "bootstrap", B = 100), fun = fun,
#' mc.cores = 10) # parallelize on 10 cores
#' }
#' 
#' ######## Example for scalar-on-function-regression with bsignal() 
#' data("fuelSubset", package = "FDboost")
#' 
#' ## center the functional covariates per observed wavelength
#' fuelSubset$UVVIS <- scale(fuelSubset$UVVIS, scale = FALSE)
#' fuelSubset$NIR <- scale(fuelSubset$NIR, scale = FALSE)
#' 
#' ## to make mboost:::df2lambda() happy (all design matrix entries < 10)
#' ## reduce range of argvals to [0,1] to get smaller integration weights
#' fuelSubset$uvvis.lambda <- with(fuelSubset, (uvvis.lambda - min(uvvis.lambda)) /
#' (max(uvvis.lambda) - min(uvvis.lambda) ))
#' fuelSubset$nir.lambda <- with(fuelSubset, (nir.lambda - min(nir.lambda)) /
#' (max(nir.lambda) - min(nir.lambda) ))
#' 
#' ## model fit with scalar response and two functional linear effects 
#' ## include no intercept as all base-learners are centered around 0    
#' 
#' mod2 <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots = 40, df = 4, check.ident = FALSE) 
#'                + bsignal(NIR, nir.lambda, knots = 40, df=4, check.ident = FALSE), 
#'                timeformula = NULL, data = fuelSubset) 
#' 
#' 
#' \dontrun{
#' bootCIs <- bootstrapCI(mod2)
#' }
#' 
#' ## run with a larger number of outer bootstrap samples
#' ## and only 10-fold for validation of each outer fold
#' ## WARNING: This may take very long!
#' \dontrun{
#' bootCIs <- bootstrapCI(mod2, B_outer = 1000, B_inner = 10)
#' }
#' 
#' @export
bootstrapCI <- function(object, which = NULL, 
                        resampling_fun_outer = NULL,
                        resampling_fun_inner = NULL,
                        B_outer = 100,
                        B_inner = 25,
                        levels = c(0.05, 0.95), 
                        plot = TRUE)
{
  
  ########## check for scalar response #########
  scalarResp <- "FDboostScalar" %in% class(object)
  # isLongFormat <- "FDboostLong" %in% class(object)

  
  
  ########## define outer resampling function if NULL #########
  if(is.null(resampling_fun_outer)){
    
    resampling_fun_outer <- function(object, fun) applyFolds(object = object,
                                                             folds = cv(rep(1, length(unique(object$id))), type =
                                                                          "bootstrap", B = B_outer), fun = fun)

  }else{
    
    message("B_outer is ignored as resampling_fun_outer is not NULL.")
    
  }
      
  ########## define inner resampling function if NULL #########
  if(is.null(resampling_fun_inner)){
    
    if(scalarResp){
      
      resampling_fun_inner <- function(object) cvrisk(object, folds = cvLong(id = object$id, weights =
                                                                               model.weights(object), B = B_inner))
      
    }else{ 
      # 
      # if(isLongFormat){
      #   
      #   resampling_fun_inner <- function(object) applyFolds()
      #   
      # }else{ # not long format
        
      resampling_fun_inner <- function(object) applyFolds(object, folds = cv(rep(1, length(unique(object$id))), type =
                                                                               "bootstrap", B = B_inner))
        
      }
      
    
  }else{
    
    message("B_inner is ignored as resampling_fun_outer is not NULL.")
    
  }
    
  # 'catch' error caused by using the cvrisk function for inner and outer resampling
  if(identical(resampling_fun_outer, cvrisk) & 
     identical(resampling_fun_inner, cvrisk)) 
    stop("Please specify a different outer resampling function.")
  
  ########## get coefficients ##########
  cat("Start computing bootstrap confidence intervals... \n")
  
  results <- resampling_fun_outer(object, 
                                  fun = function(mod)
                                  {
                                    
                                    ms <- mstop( resampling_fun_inner(mod) )
                                    coefs <- coef( mod[ms], which = which )
                                    
                                    return(list(coefs = coefs, ms = ms))
                                    
                                  }
                                  
  )
  
  cat("\n")
  
  coefs <- lapply(results, "[[", "coefs")
  mstops <- sapply(results, "[[", "ms")
  offsets <- t(sapply(coefs, function(x) x$offset$value))
  
  
  ########## format coefficients #########
  # number of baselearners
  nrEffects <- length(coefs[[1]]$smterms)
  
  # extract values
  listOfCoefs <- lapply(1:nrEffects, function(i) lapply(1:length(coefs), 
                                                        function(j) coefs[[j]]$smterms[[i]]$value))
  
  names(listOfCoefs) <- names(object$baselearner)
  
  # check for effect surfaces
  isSurface <- sapply(1:nrEffects, function(i) !is.null(coefs[[1]]$smterms[[i]]$y) )
  
  # reduce lists for non surface effects
  listOfCoefs[!isSurface] <- lapply(listOfCoefs[!isSurface], function(x) do.call("rbind", x))
  
  listOfCoefs <- c(offsets = list(offsets), listOfCoefs)
  isSurface <- c(FALSE, isSurface)

  # add information about the values of the covariate
  # and change format
  for(i in 1:length(listOfCoefs)){
    
    if(i!=1){
      
      atx <- coefs[[1]]$smterms[[i-1]]$x # i-1 because of the offset
      
    }else{
      
      atx <- coefs[[1]]$offset$x
      
    }
    
    aty <- NA
    if(isSurface[i]) aty <- coefs[[1]]$smterms[[i-1]]$y # i-1 because of the offset

    # format functional factors
    if(is.list(listOfCoefs[[i]]) & is.factor(atx)){
      
      # combine each factor level
      listOfCoefs[[i]] <- lapply(1:length(levels(atx)),
                                 function(faclevnr) sapply(listOfCoefs[[i]], function(x) x[faclevnr,]))
      isSurface[i] <- FALSE
      
    }else if(is.list(listOfCoefs[[i]])){ # effect surfaces
      
      listOfCoefs[[i]] <- do.call("rbind", lapply(listOfCoefs[[i]],c))
      
    }
    
    attr(listOfCoefs[[i]], "x") <- atx
    if(!is.na(sum(aty))) attr(listOfCoefs[[i]], "y") <- aty
    
  }

  # return raw results
  if(is.null(levels)) return(listOfCoefs)

  
  
  ########## calculate quantiles #########
  # create list for quantiles
  listOfQuantiles <- vector("list", length(listOfCoefs))
  
  # calculate quantiles
  for(i in 1:length(listOfCoefs)){
    
    # for matrix object
    if(is.matrix(listOfCoefs[[i]])){
      
      listOfQuantiles[[i]] <- apply(listOfCoefs[[i]], 2, quantile, probs = levels)
      attr(listOfQuantiles[[i]], "x") <- attr(listOfCoefs[[i]], "x")
      if(!is.null(attr(listOfCoefs[[i]], "y")))
        attr(listOfQuantiles[[i]], "y") <- attr(listOfCoefs[[i]], "y")

    }else if(is.list(listOfCoefs[[i]])){ # functional factor variables
      
      listOfQuantiles[[i]] <- lapply(listOfCoefs[[i]], function(x) apply(x, 1, quantile, probs = levels))
      
    }else{# scalar case
      
      listOfQuantiles[[i]] <- quantile(listOfCoefs[[i]], probs = levels)
      
    }
         
  }
  
  # since coefficient surfaces are saved as vectors, reconstruct quantiles 
  # as coefficient surfaces and return a list of matrices, where each 
  # matrix corresponds to a quantile in levels
  if(sum(isSurface)!=0) listOfQuantiles[which(isSurface)] <- 
    lapply(listOfQuantiles[isSurface], 
           function(x){ 
             
             retL <- lapply(1:nrow(x), function(i) 
               matrix(x[i,], nrow = length(attr(x, "y"))))
             names(retL) <- levels
             return(retL)
             
           })
  
  # name rows of the matrices for non-surface effects
  if(sum(!isSurface)!=0) listOfQuantiles[which(!isSurface)] <- 
    lapply(listOfQuantiles[!isSurface],
           function(x){
             
             if(is.list(x)){
               
               for(j in 1:length(x)){
                 
                 if(!is.null(dim(x[[j]]))){
                   rownames(x[[j]]) <- levels
                 }else{
                   names(x[[j]]) <- levels
                 }
                 
               }
               
             }else{
             
               if(!is.null(dim(x))){
                 rownames(x) <- levels
               }else{
                 names(x) <- levels
               }
             
             }
             
             return(x)
             
           })
  
  # save names of baselearners
  names(listOfQuantiles) <- names(listOfCoefs)


  ########## return results #########
  ret <- list(raw_results = listOfCoefs,
              quantiles = listOfQuantiles,
              mstops = mstops,
              resampling_fun_outer = resampling_fun_outer,
              resampling_fun_inner = resampling_fun_inner,
              B_outer = B_outer,
              B_inner = B_inner,
              which = which,
              levels = levels)
  
  class(ret) <- "bootstrapCI"
  
  return(ret)
  
}


#' Plot a bootstrapCI object
#' 
#' Takes a \code{bootstrapCI}-object and produces ggplot objects
#' 
#' @param x a fitted \code{bootstrapCI}-object
#' @param ... ignored
#' @seealso \code{\link{bootstrapCI}}
#' @return a list of \code{ggplot} objects
#' @method plot bootstrapCI
#' 
#' @note Quantiles are plotted using \code{ggplot}.
#' @description For scalar non-functional covariates, a boxplot is 
#' used to visualize the distribution. Quantiles are drawn by 
#' vertical red lines. For scalar response, functional covariates
#' are visualized by line plots with quantiles in dashed red.
#' For function-on-function effect surfaces, the quantiles are plotted
#' as surface plots.
#' 
#' @export
plot.bootstrapCI <- function(x, ...)
{

  stopifnot(class(x)=="bootstrapCI")
  if(!requireNamespace("ggplot2")) stop("Please install ggplot2.")
  
  ### prepare objects
  
  plotObjList <- vector("list", length = length(x$raw_results))
  qsclass <- sapply(x$quantiles, class)
  levels <- x$levels
  type <- rep(NA, length(x$raw_results))
  
  for(i in 1:length(x$raw_results)){
 
    if(qsclass[i]=="matrix"){
      
      rr <- x$raw_results[[i]]
      qs <- x$quantiles[[i]]
      
      plotObjList[[i]] <- 
        rbind(data.frame(value = c(t(rr)),
                         group = rep(1:ncol(rr), each = nrow(rr)),
                         x = attr(rr, "x"),
                         what = "raw"),
              data.frame(value = c(t(qs)),
                         group = rep(1:ncol(qs), each = nrow(qs)),
                         x = attr(rr, "x"),
                         what = "quantile")
        )
      
      type[i] <- "line"
      
    }else if(qsclass[i]=="list"){
      
      rr <- x$raw_results[[i]]
      qs <- x$quantiles[[i]]
      
      if(is.factor(attr(rr, "x"))){
        # second dimension is not numeric but a factor
        
        # get all dimensions
        lenf <- length(levels(attr(rr, "x")))
        lennum <- length(attr(rr, "y"))
        lenboot <- ncol(rr[[1]])
        lenlev <- length(levels)

        plotObjList[[i]] <- rbind(
          data.frame(value = c(sapply(rr, c)),
                     group = rep(rep(1:lenboot, each = lennum), lenf),
                     fac = rep(attr(rr, "x"), each = lennum * lenboot),
                     x = rep(attr(rr, "y"), lenf * lenboot),
                     what = "raw"),
          data.frame(value = c(sapply(qs, function(x) c(t(x)))),
                     group = rep(rep(1:lenlev, each = lennum), lenf),
                     fac = rep(attr(rr, "x"), each = lennum * lenlev),
                     x = rep(attr(rr, "y"), lenf * lenlev),
                     what = "quantile")
        )
        
        type[i] <- "multLine"

      }else{ # effect with two numeric dimensions

        plotObjList[[i]] <- data.frame(coefficient = c(sapply(qs, function(x) c(t(x)))),
                                       # group = rep(1:ncol(qs[[1]]), each = nrow(qs[[1]])),
                                       x = rep(attr(rr, "x"), length(attr(rr, "y"))),
                                       y = rep(attr(rr, "y"), each = length(attr(rr, "x"))),
                                       levels = rep(levels, each = prod(dim(qs[[1]]))))
        
        type[i] <- "surface"
        
      }
              
    }else{
      
      rr <- x$raw_results[[i]]
      
      plotObjList[[i]] <- data.frame(value = rr)
      
      type[i] <- "boxplot"
      
    } 
      
  }
  
  # define ggplot function for every kind of type
  plotFun <- function(type) 
    switch(type,
           
           line = function(obj) 
             ggplot(obj) + 
             geom_line(data = obj[obj$what=="raw",], 
                       aes(x = x, y = value, group = group),
                       alpha = 0.67, colour = "grey50") + 
             geom_line(data = obj[obj$what=="quantile",],
                       aes(x = x, y = value, group = group),
                       colour = "red", linetype = "dashed", size = 1.1) + 
             xlab("value") + ylab("coefficient"),
           
           multLine = function(obj)
             ggplot(obj) + 
             geom_line(data = obj[obj$what=="raw",], 
                       aes(x = x, y = value, group = group),
                       alpha = 0.67, colour = "grey50") + 
             geom_line(data = obj[obj$what=="quantile",],
                       aes(x = x, y = value, group = group),
                       colour = "red", linetype = "dashed", size = 1.1) + 
             facet_wrap(~ fac) + 
             xlab("value") + ylab("coefficient"),
           
           surface = function(obj) 
             ggplot(obj, aes(x = x, y = y, z = coefficient, fill = coefficient)) + 
             geom_tile() + 
             scale_fill_distiller(palette = "Spectral") + 
             stat_contour(binwidth = abs((max(obj$coefficient) - min(obj$coefficient))/10),
                          col = "black") +
             facet_wrap(~levels) + theme_bw() + xlab("t") + ylab("s"),
           
           boxplot = function(obj) 
             ggplot(obj) + 
             geom_boxplot() + 
             geom_abline(yintercept = quantiles(rr, probs = levels), col = "red") + 
             ylab("coefficient")
    )

  # call custom ggplot functions
  ggList <- lapply(1:length(type), function(i) plotFun(type[i])(obj = plotObjList[[i]]) + 
                     ggtitle(names(x$raw_results)[i]))
  
  # print ggplot objects
  lapply(ggList, print)
  
  # return
  invisible(ggList)
  
}






#' Print a bootstrapCI object
#' 
#' Takes a \code{bootstrapCI}-object and produces a print on the console.
#' 
#' @param x a fitted \code{bootstrapCI}-object
#' @param ... currently not used
#' @seealso \code{\link{bootstrapCI}}
#' @return a list with information on the model 
#' @method print bootstrapCI
#' @export
print.bootstrapCI <- function(x, ...)
{

  stopifnot(class(x)=="bootstrapCI")

  cat("\n")
     
  cat("\t Bootstrapped confidence interval object of FDboost fit\n")
  
  cat("\n")
  cat("Coefficients:\n\t", names(x$quantiles), sep="\t", fill = TRUE)
  cat("\n")
  cat("\n")
  cat("Summary for stopping iterations of inner validation:\n\n")
  print(summary(x$mstops))
  cat("\n")
  invisible(x)
    
}