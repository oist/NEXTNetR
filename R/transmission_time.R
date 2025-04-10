#' @name time_distributions
#' @title Creating time distributions
#'
#' @description
#' NEXTNetR supports a series of pre-defined distributions for transmission and
#' recover/reset time, and also allows custom distributions to be defined
#' with [userdefined_time]. 
#' 
#' @returns
#' * `exponential_time(lambda, p_infinity)`. Returns a time distribution representing
#'   an exponential distribution with rate `lambda` which in addition to finite
#'   values produces the value infinity with probability `p_infinity`.
#' * `lognormal_time(mean, var, p_infinity)`. Returns a time distribution representing
#'   a Log-normal distribution with the given `mean` and `variance`, which in addition
#'   to finite  values produces the value infinity with probability `p_infinity`.
#' * `gamma_time(mean, var, p_infinity)`. Returns a time distribution representing
#'   a Gamma distribution with the given `mean` and `variance`, which in addition
#'   to finite  values produces the value infinity with probability `p_infinity`.
#' * `weibull_time(mean, var, p_infinity)`. Returns a time distribution representing
#'   a Weibull distribution with the given `shape` and `scale` parameter, which in addition
#'   to finite  values produces the value infinity with probability `p_infinity`.
#'   The distribution has mean \eqn{b \Gamma(1 + 1/a)} and variance 
#'   \eqn{b^2 \Gamma(1 + 1/a) - \Gamma^2(1 + 1/a)}for shape \eqn{a} and scale \eqn{b}.
#' * `polynomial_rate_time(coeffs)`. Distribution with survival function
#'   \eqn{\Psi(\tau) = e^{-p(\tau)}} for a polynomial hazard rate
#'   \eqn{p = c[1] + c[2] x + c[3] x^2 + \ldots} with non-negative coefficients.
#' * `deterministic_time(tau)`. Deterministic time with fixed value `tau`.
NULL

#' @rdname time_distributions
#' @export
exponential_time <- function(lambda, p_infinity = 0.0) {
  nextnetR_exponential_time(as.double(lambda), as.double(p_infinity))
}

#' @rdname time_distributions
#' @export
lognormal_time <- function(mean, var, p_infinity = 0.0) {
  nextnetR_gamma_time(as.double(mean), as.double(var), as.double(p_infinity))
}

#' @rdname time_distributions
#' @export
gamma_time <- function(mean, var, p_infinity) {
  nextnetR_gamma_time(as.double(mean), as.double(var), as.double(p_infinity))
}

#' @rdname time_distributions
#' @export
weibull_time <- function(shape, scale, p_infinity) {
  nextnetR_gamma_time(as.double(shape), as.double(scale), as.double(p_infinity))
}

#' @rdname time_distributions
#' @export
polynomial_rate_time <- function(coeffs) {
  nextnetR_polynomial_rate_time(as.double(coeffs))
}

#' @rdname time_distributions
#' @export
deterministic_time <- function(tau) {
  nextnetR_deterministic_time(as.doubl(tau))
}

#' @title User-defined time distributions
#' 
#' @description TODO
#' 
#' @export
userdefined_time <- function(density, survivalprobability, probability_is_trinary,
                             survivalquantile, quantile_is_trinary,
                             sample, p_infinity)
{
  nextnetR_userdefined_time(as.function(density), as.function(survivalprobability),
                            as.logical(probability_is_trinary),
                            as.function(survivalquantile), as.logical(quantile_is_trinary),
                            as.function(sample), as.double(p_infinity))
}

#' @name time_functions
#' @title Functions that operate on transmission/recovery/reset time distributions
#' 
#' @description
#' A time distributions represents a family of distributions parametrized by
#' a conditioning time \eqn{t \geq 0}, a multiplicity \eqn{m > 0} over
#' a base distribution with survival function
#' \eqn{\Psi(\tau) = e^{-\lambda(\tau)}}. The parameters represent (i)
#' conditioning of the base distribution on values \eqn{\geq t} and (ii)
#' taking the *minimum* over \eqn{m} i.i.d. samples of \eqn{\Psi}. The survival
#' function of \eqn{\Psi_{t,m}} is therefore \eqn{(\Psi(\tau)/\Psi(t))^m}.
#' 
#' @param timedistribution a [time_distribution][time_distributions].
#' @param n number of samples to draw
#' @param tau non-negative time
#' @param t the conditioning time
#' @param m multiplicity
#' 
#' @seealso \code{\link{time_distributions}}
#' 
#' @returns
#' * `time_sample(n, timedistribution, t, m)`. Samples `n` values from
#'    distribution \eqn{\Psi_{t,m}}.
#'   
#' * `time_density(timedistribution, tau)`. Evaluates the density at points
#'    `tau` of distribution \eqn{\Psi_{t,m}}.
#'   
#' * `time_hazardrate(timedistribution, tau)`. Evaluates the hazardrate
#'    \eqn{\lambda(\tau)} of the *base* distribution \eqn{\Psi}. The hazard
#'    rate of \eqn{\Psi_{t,m}} is simply \eqn{m \times \lambda(\tau)}.
#'   
#' * `time_survivalprobability(timedistribution, tau, t, m)`. Evaluates the
#'    survivalfunction \eqn{\Psi_{t,m}} at points `tau`.
#'   
#' * `time_survivalquantile(timedistribution, p, t, m)`. Computes the
#'    \eqn{p}-quantiles of distribution \eqn{\Psi_{t,m}}.
NULL

#' @rdname time_functions
#' @export
time_sample <- function(n, timedistribution, t=0, m=1) {
  nextnetR_time_sample(as.integer(n), timedistribution, as.double(t), as.integer(m))
}

#' @rdname time_functions
#' @export
time_density <- function(timedistribution, tau) {
  nextnetR_time_density(timedistribution, as.double(tau))
}

#' @rdname time_functions
#' @export
time_hazardrate <- function(ttimedistributiontr, tau) {
  nextnetR_time_hazardrate(timedistribution, as.double(tau))
}

#' @rdname time_functions
#' @export
time_survivalprobability <- function(timedistribution, tau, t = 0, m = 1) {
  tau <- as.double(tau)
  t <- as.double(t)
  m <- as.integer(m)
  if (length(t) == 1)
    t <- rep(t, length(tau))
  if (length(m) == 1)
    m <- rep(m, length(tau))
  nextnetR_time_survivalprobability(timedistribution, tau, t, m)
}

#' @rdname time_functions
#' @export
time_survivalquantile <- function(timedistribution, p, t = 0, m = 1) {
  p <- as.double(p)
  t <- as.double(t)
  m <- as.integer(m)
  if (length(t) == 1)
    t <- rep(t, length(p))
  if (length(m) == 1)
    m <- rep(m, length(p))
  nextnetR_time_survivalquantile(timedistribution, p, m, t)
}

