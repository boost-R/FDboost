
#' Factorize tensor product effects
#' 
#' @export
#' @name factorize

factorize <- factorise <- function(x, ...) {
  UseMethod("factorize")
}

factorize.default <- function(x, ...) stop("No default factorization.")