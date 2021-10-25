#' Clr and inverse clr transformation
#' 
#' \code{clr} computes the clr or inverse clr transformation of a vector \code{f}
#' with respect to integration weights \code{w}, corresponding to a Bayes Hilbert space
#' \eqn{B^2(\mu) = B^2(\mathcal{T}, \mathcal{A}, \mu)}{B^2(\mu) = B^2(T, A, \mu)}.
#' 
#' @param f a vector containing the function values (evaluated on a grid) of the 
#' function \eqn{f} to transform. If \code{inverse = TRUE}, \code{f} must be a density,
#' i.e., all entries must be positive and usually \code{f} integrates to one. 
#' If \code{inverse = FALSE}, \code{f} should integrate to zero, see Details.
#' @param w a vector of length one or of the same length as \code{f} containing 
#' positive integration weights. If \code{w} has length one, this
#' weight is used for all function values. The integral of \eqn{f} is approximated
#' via \eqn{\int_{\mathcal{T}} f \, \mathrm{d}\mu \approx 
#' \sum_{j=1}^m}{\sum_{j=1}^m} \code{w}\eqn{_j} \code{f}\eqn{_j},
#' where \eqn{m} equals the length of \code{f}.
#' @param inverse if \code{TRUE}, the inverse clr transformation is computed.
#' 
#' @details The clr transformation maps a density \eqn{f} from \eqn{B^2(\mu)} to
#' \eqn{L^2_0(\mu) := \{ f \in L^2(\mu) ~|~ \int_{\mathcal{T}} f \, \mathrm{d}\mu = 0\}}{L^2_0(\mu) := {f \in L^2(\mu) | \int_T f d\mu = 0}}
#' via
#' \deqn{\mathrm{clr}(f) := \log f - \frac{1}{\mu (\mathcal{T})} \int_{\mathcal{T}} \log f \, \mathrm{d}\mu.}{clr(f) := log f - 1/\mu(T) * \int_T log f d\mu.}
#' The inverse clr transformation maps a function \eqn{f} from
#' \eqn{L^2_0(\mu)} to \eqn{B^2(\mu)} via
#' \deqn{\mathrm{clr}^{-1}(f) := \frac{\exp f}{\int_{\mathcal{T}} \exp f \, \mathrm{d}\mu}.}{clr^{-1}(f) := (exp f) / (\int_T \exp f d\mu).}
#' Note that in contrast to Maier et al. (2021), this definition of the inverse
#' clr transformation includes normalization, yielding the respective probability 
#' density function (representative of the equivalence class of proportional
#' functions in \eqn{B^2(\mu)}). 
#' 
#' The (inverse) clr transformation depends not only on \eqn{f}, but also on the
#' underlying measure space \eqn{\left( \mathcal{T}, \mathcal{A}, \mu\right)}{(T, A, \mu)}, 
#' which determines the integral. In \code{clr} this is specified via the 
#' integration weights \code{w}. E.g., for a discrete set \eqn{\mathcal{T}}{T}
#' with \eqn{\mathcal{A} = \mathcal{P}(\mathcal{T})}{A = P(T)} the power set of 
#' \eqn{\mathcal{T}}{T} and \eqn{\mu = \sum_{t \in T} \delta_t} the sum of dirac
#' measures at \eqn{t \in \mathcal{T}}{t \in T}, the default \code{w = 1} is
#' the correct choice. In this case, integrals are indeed computed exactly, not
#' only approximately. 
#' For an interval \eqn{\mathcal{T} = [a, b]}{T = [a, b]}
#' with \eqn{\mathcal{A} = \mathcal{B}}{A = B} the Borel \eqn{\sigma}-algebra 
#' restricted to \eqn{\mathcal{T}}{T} and \eqn{\mu = \lambda} the Lebesgue measure,
#' the choice of \code{w} depends on the grid on which the function was evaluated:
#' \code{w}\eqn{_j} must correspond to the length of the subinterval of \eqn{[a, b]}, which 
#' \code{f}\eqn{_j} represents.
#' E.g., for a grid with equidistant distance \eqn{d}, where the boundary grid 
#' values are \eqn{a + \frac{d}{2}}{a + d/2} and \eqn{b - \frac{d}{2}}{b - d/2}
#' (i.e., the grid points are centers of intervals of size \eqn{d}),
#' equal weights \eqn{d} should be chosen for \code{w}. 
#' 
#' The clr transformation is crucial for density-on-scalar regression 
#' since estimating the clr transformed model in \eqn{L^2_0(\mu)} is equivalent
#' to estimating the original model in \eqn{B^2(\mu)} (as the clr transformation
#' is an isometric isomorphism), see also the vignette "FDboost_density-on-scalar_births"
#' and Maier et al. (2021).
#' 
#' @return A vector of the same length as \code{f} containing the (inverse) clr 
#' transformation of \code{f}.
#' 
#' @author Eva-Maria Maier
#' 
#' @references 
#' Maier, E.-M., Stoecker, A., Fitzenberger, B., Greven, S. (2021):
#' Additive Density-on-Scalar Regression in Bayes Hilbert Spaces with an Application to Gender Economics.
#' arXiv preprint arXiv:2110.11771.
#' 
#' @examples
#' ### Continuous case (T = [0, 1] with Lebesgue measure):
#' # evaluate density of a Beta distribution on an equidistant grid
#' g <- seq(from = 0.005, to = 0.995, by = 0.01)
#' f <- dbeta(g, 2, 5)
#' # compute clr transformation with distance of two grid points as integration weight
#' f_clr <- clr(f, w = 0.01)
#' # visualize result
#' plot(g, f_clr , type = "l")
#' abline(h = 0, col = "grey")
#' # compute inverse clr transformation (w as above)
#' f_clr_inv <- clr(f_clr, w = 0.01, inverse = TRUE)
#' # visualize result
#' plot(g, f, type = "l")
#' lines(g, f_clr_inv, lty = 2, col = "red")
#' 
#' ### Discrete case (T = {1, ..., 12} with sum of dirac measures at t in T):
#' data("birthDistribution", package = "FDboost")
#' # fit density-on-scalar model with effects for sex and year
#' model <- FDboost(birth_densities_clr ~ 1 + bolsc(sex, df = 1) + 
#'                    bbsc(year, df = 1, differences = 1),
#'                  # use bbsc() in timeformula to ensure integrate-to-zero constraint
#'                  timeformula = ~bbsc(month, df = 4, 
#'                                      # December is followed by January of subsequent year
#'                                      cyclic = TRUE, 
#'                                      # knots = {1, ..., 12} with additional boundary knot
#'                                      # 0 (coinciding with 12) due to cyclic = TRUE
#'                                      knots = 1:11, boundary.knots = c(0, 12), 
#'                                      # degree = 1 with these knots yields identity matrix 
#'                                      # as design matrix
#'                                      degree = 1),
#'                  data = birthDistribution, offset = 0, 
#'                  control = boost_control(mstop = 1000))
#' # Extract predictions (clr-transformed!) and transform them to Bayes Hilbert space
#' predictions_clr <- predict(model)
#' predictions <- t(apply(predictions_clr, 1, clr, inverse = TRUE))
#' 
#' @export
clr <- function(f, w = 1, inverse = FALSE) {
  stopifnot("inverse must be TRUE or FALSE." = inverse %in% c(TRUE, FALSE))
  stopifnot("f must be numeric." = is.numeric(f))
  stopifnot("f contains missing values." = !(anyNA(f)))
  n <- length(f)
  if (length(w) == 1) {
    w <- rep(w, n)
  }
  stopifnot("w must be contain positive weights of length 1 or the same length as f." = 
              length(w) == n & is.numeric(w) & all(w > 0))
  int_f <- sum(f * w)
  if (!inverse) {
    stopifnot("As a density, f must be positive." = all(f > 0))
    if (!isTRUE(all.equal(int_f, 1, tolerance = 0.01))) {
      warning(paste0("f is not a probability density with respect to w. Its integral is ",
                     round(int_f, 2), "."))
    }
    return(log(f) - 1 / sum(w) * sum(log(f) * w))
  } else {
    if (!isTRUE(all.equal(int_f, 0, tolerance = 0.01))) {
      warning(paste0("f does not integrate to zero with respect to w. Its integral is ",
                     round(int_f, 2), "."))
    }
    return(exp(f) / sum(exp(f) * w))
  }
}


#' Densities of live births in Germany 
#'
#' \code{birthDistribution} contains densities of live births in Germany over the
#' months per year (1950 to 2019) and sex (male and female), resulting in 140
#' densities.
#'
#' @docType data
#'
#' @usage data(birthDistribution, package = "FDboost")
#'
#' @format A list in the correct format to be passed to \code{\link{FDboost}} for 
#' density-on-scalar regression:
#' \describe{
#' \item{\code{birth_densities}}{A 140 x 12 matrix containing the birth densities
#' in its rows. The first 70 rows correspond to male newborns, the second 70 rows 
#' to female ones. Within both of these, the years are ordered increasingly 
#' (1950-2019), see also \code{sex} and \code{year}.}
#' \item{\code{birth_densities_clr}}{A 140 x 12 matrix containing the clr
#' transformed densities in its rows. Same structure as \code{birth_densities}.}
#' \item{\code{sex}}{A factor vector of length 140 with levels \code{"m"} (male)
#' and \code{"f"} (female), corresponding to the sex of the newborns for the rows of
#' \code{birth_densities} and \code{birth_densities_clr}. The first 70 elements 
#' are \code{"m"}, the second 70 \code{"f"}.}
#' \item{\code{year}}{A vector of length 140 containing the integers from 1950 
#' to 2019 two times (\code{c(1950:2019, 1950:2019)}), corresponding to the years
#' for the rows of \code{birth_densities} and \code{birth_densities_clr}.}
#' \item{\code{month}}{A vector containing the integers from 1 to 12, corresponding
#' to the months for the columns of \code{birth_densities} and \code{birth_densities_clr}
#' (domain \eqn{\mathcal{T}}{T} of the (clr-)densities).}
#' }
#' Note that for estimating a density-on-scalar model with \code{FDboost}, the
#' clr transformed densities (\code{birth_densities_clr}) serve as response, see
#' also the vignette "FDboost_density-on-scalar_births".
#' The original densities (\code{birth_densities}) are not needed for estimation,
#' but still included for the sake of completeness.
#' 
#' @details To compensate for the different lengths of the months, the average 
#' number of births per day for each month (by sex and year) was used to compute
#' the birth shares from the absolute birth counts. The 12 shares corresponding
#' to one year and sex form one density in the Bayes Hilbert space 
#' \eqn{B^2(\delta) = B^2\left( \mathcal{T}, \mathcal{A}, \delta\right)}{B^2(\delta) = B^2(T, A, \delta)},
#' where \eqn{\mathcal{T} = \{1, \ldots, 12\}}{T = {1, \ldots, 12}} corresponds
#' to the set of the 12 months, \eqn{\mathcal{A} := \mathcal{P}(\mathcal{T})}{A := P(T)}
#' corresponds to the power set of \eqn{\mathcal{T}}{T}, and the reference measure
#' \eqn{\delta := \sum_{t = 1}^{12} \delta_t} corresponds to the sum of dirac 
#' measures at \eqn{t \in \mathcal{T}}{t \in T}.
#' 
#' @seealso \code{\link{clr}} for the (inverse) clr transformation.
#'
#' @source Statistisches Bundesamt (Destatis), Genesis-Online, data set 
#' \href{https://www-genesis.destatis.de/genesis//online?operation=table&code=12612-0002&bypass=true&levelindex=0&levelid=1610983595176#abreadcrumb}{12612-0002} 
#' (01/18/2021); \href{https://www.govdata.de/dl-de/by-2-0}{dl-de/by-2-0}; 
#' processed by Eva-Maria Maier
#' 
#' @references 
#' Maier, E.-M., Stoecker, A., Fitzenberger, B., Greven, S. (2021):
#' Additive Density-on-Scalar Regression in Bayes Hilbert Spaces with an Application to Gender Economics.
#' arXiv preprint arXiv:2110.11771.
#' 
#' @examples
#' data("birthDistribution", package = "FDboost")
#' 
#' # Plot densities
#' year_col <- rainbow(70, start = 0.5, end = 1)
#' year_lty <- c(1, 2, 4, 5)
#' par(mfrow = c(1, 2))
#' funplot(1:12, birthDistribution$birth_densities[1:70, ], ylab = "densities", xlab = "month", 
#'         xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Male")
#' funplot(1:12, birthDistribution$birth_densities[71:140, ], ylab = "densities", xlab = "month", 
#'         xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Female")
#' par(mfrow = c(1, 1))
#' 
#' # fit density-on-scalar model with effects for sex and year
#' model <- FDboost(birth_densities_clr ~ 1 + bolsc(sex, df = 1) + 
#'                    bbsc(year, df = 1, differences = 1),
#'                  # use bbsc() in timeformula to ensure integrate-to-zero constraint
#'                  timeformula = ~bbsc(month, df = 4, 
#'                                      # December is followed by January of subsequent year
#'                                      cyclic = TRUE, 
#'                                      # knots = {1, ..., 12} with additional boundary knot
#'                                      # 0 (coinciding with 12) due to cyclic = TRUE
#'                                      knots = 1:11, boundary.knots = c(0, 12), 
#'                                      # degree = 1 with these knots yields identity matrix 
#'                                      # as design matrix
#'                                      degree = 1),
#'                  data = birthDistribution, offset = 0, 
#'                  control = boost_control(mstop = 1000))
#' 
#' # Plotting 'model' yields the clr-transformed effects
#' par(mfrow = c(1, 3))
#' plot(model, n1 = 12, n2 = 12)
#' 
#' # Use inverse clr transformation to get effects in Bayes Hilbert space, e.g. for intercept
#' intercept_clr <- predict(model, which = 1)[1, ]
#' intercept <- clr(intercept_clr, w = 1, inverse = TRUE)
#' funplot(1:12, intercept, xlab = "month", xaxp = c(1, 12, 11), pch = 20,
#'         main = "Intercept", ylab = expression(hat(beta)[0]), id = rep(1, 12))
#' 
#' # Same with predictions
#' predictions_clr <- predict(model)
#' predictions <- t(apply(predictions_clr, 1, clr, inverse = TRUE))
#' pred_ylim <- range(birthDistribution$birth_densities)
#' par(mfrow = c(1, 2))
#' funplot(1:12, predictions[1:70, ], ylab = "predictions", xlab = "month", ylim = pred_ylim,
#'         xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Male")
#' funplot(1:12, predictions[71:140, ], ylab = "predictions", xlab = "month", ylim = pred_ylim,
#'         xaxp = c(1, 12, 11), pch = 20, col = year_col, lty = year_lty, main = "Female")
#' par(mfrow = c(1, 1))
"birthDistribution"