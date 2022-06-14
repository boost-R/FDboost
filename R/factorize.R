
#' Factorize tensor product effects
#' 
#' @param x input 
#' @param ... additional parameters
#' 
#' @export
#' @name factorize

factorize <- factorise <- function(x, ...) {
  UseMethod("factorize")
}

factorize.default <- function(x, ...) stop("No default factorization.")
