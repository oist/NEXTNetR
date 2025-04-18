#' @name time_distributions
#' @title Creating time distributions
#'
#' @description
#' NEXTNetR supports a series of pre-defined distributions for transmission and
#' recover/reset time, and also allows custom distributions to be defined
#' with [userdefined_time].
#' 
#' Each `time_distribution` object actually represents a two-parameer family
#' of distributions, see [time_functions] for a full discussion and for functions
#' that operate on time distribution objects.
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
#'   \eqn{b^2(\Gamma(1 + 2/a) - \Gamma^2(1 + 1/a))}for shape \eqn{a} and scale \eqn{b}.
#' * `polynomial_rate_time(coeffs)`. Distribution with survival function
#'   \eqn{\Psi(\tau) = e^{-p(\tau)}} for a polynomial hazard rate
#'   \eqn{p = c[1] + c[2] x + c[3] x^2 + \ldots} with non-negative coefficients.
#' * `deterministic_time(tau)`. Deterministic time with fixed value `tau`
#' 
#' @seealso \code{\link{time_functions}}
#' 
NULL

#' @rdname time_distributions
#' @export
exponential_time <- function(lambda, p_infinity = 0.0) {
  nextnetR_exponential_time(as.double(lambda), as.double(p_infinity))
}

#' @rdname time_distributions
#' @export
lognormal_time <- function(mean, var, p_infinity = 0.0) {
  nextnetR_lognormal_time(as.double(mean), as.double(var), as.double(p_infinity))
}

#' @rdname time_distributions
#' @export
gamma_time <- function(mean, var, p_infinity = 0.0) {
  nextnetR_gamma_time(as.double(mean), as.double(var), as.double(p_infinity))
}

#' @rdname time_distributions
#' @export
weibull_time <- function(shape, scale, p_infinity = 0.0) {
  nextnetR_weibull_time(as.double(shape), as.double(scale), as.double(p_infinity))
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
#' @seealso \code{\link{time_distributions}}, \code{\link{time_functions}}
#' 
#' @export
userdefined_time <- function(sample=NULL, sample_is_trinary, survival, survival_is_trinary,
                             density=NULL, survivalquantile=NULL, quantile_is_trinary,
                             pinfinity=0.0)
{
  nextnetR_userdefined_time(
    if (!is.null(sample)) as.function(sample) else NULL,
    if (!is.null(sample)) as.logical(sample_is_trinary) else FALSE,
    as.function(survival), as.logical(survival_is_trinary),
    if (!is.null(density)) as.function(density) else NULL,
    if (!is.null(survivalquantile)) as.function(survivalquantile) else NULL,
    if (!is.null(survivalquantile)) as.logical(quantile_is_trinary) else FALSE,
    pinfinity=pinfinity)
}

#' @name time_functions
#' @title Functions that operate on transmission/recovery/reset time distributions
#' 
#' @description
#' Time distribution objects are used to specify the distribution of transmission
#' and recovery/reset times. In principle, time distribution can represent any
#' distribution on \eqn{[0,\infty]}, where the value \eqn{\infty} indicates that no
#' transmission respectively recovery/reset takes place. See [time_distributions]
#' for functions to create objects representing common time distributions such
#' as lognormal, gamma, etc. Each time distribution object actually represents
#' a two-parameter family of distributions with parameters \eqn{t} and \eqn{m},
#' see Details for a full discussion.
#' 
#' @details
#' Each time distribution object actually represents a two-parameter family of
#' distributions with parameters \eqn{t} and \eqn{m}. This family is defined in
#' terms of a base distribution with some survival function
#' \eqn{\Psi(\tau) = \exp\big(-\int_0^\tau \lambda(t) dt\big)}.
#' Note that assuming this form of the survival function does not restrict the choice
#' of base distributions; any distribution on \eqn{[0,\infty]} takes this form
#' for \eqn{\lambda(\tau) := -\Psi'(\tau) / \Psi(\tau)}. \eqn{\lambda(\tau)} is
#' called the *hazard rate* at time \eqn{\tau}.
#' 
#' In terms of hazard rate \eqn{\lambda(\tau)}, the distribution \eqn{\Psi(\tau)}
#' is naturally interpreted as the time until the first event fired by an
#' inhomogeneous Poisson process with rate function \eqn{\lambda(\tau)}. If 
#' \eqn{\int_0^\tau \lambda(t) dt = \Lambda < \infty} for all \eqn{\tau},
#' then the process does not necessarily fire and the distribution takes
#' value \eqn{\infty} with positive probability \eqn{p_\infty = e^{-\Lambda}}.
#' 
#' Parameter \eqn{m} modulates the rate function \eqn{\lambda(\tau)}, and
#' parameter \eqn{t} conditions the process to fire after time \eqn{t}. Note
#' that for conditioned distributions, \eqn{\tau} is expressed *relative*
#' to \eqn{t}, i.e. the domain of \eqn{\tau} is always \eqn{[0,\infty]}. The
#' survival function of the modulates and conditioned distribution is
#' \eqn{\Psi_{t,m}(\tau) = \exp\big(-\int_t^{t+\tau} m \lambda(t) dt\big)}.
#' 
#' For integral \eqn{m}, modulation can be interpreted as taking the minimum
#' of \eqn{m} i.i.d. copies of the unmodulated distribution. This is used by
#' some simulation algorithms to efficiently handle nodes with large numbers
#' of edges. Parameter \eqn{m} is also used when simulating on *weighted*
#' networks, where it represents the weight of an edge.
#' 
#' @param timedistribution a [time_distribution][time_distributions].
#' @param n number of samples to draw
#' @param tau non-negative time
#' @param t the conditioning time, see details
#' @param m multiplicity, see details
#' 
#' @seealso \code{\link{time_distributions}}
#' 
#' @returns
#' * `time_sample(n, timedistribution, t, m)`. Samples `n` values from
#'    distribution \eqn{\Psi_{t,m}}.
#'   
#' * `time_density(timedistribution, tau, t, m)`. Evaluates the density at points
#'    `tau` of distribution \eqn{\Psi_{t,m}}.
#'   
#' * `time_hazardrate(timedistribution, tau)`. Evaluates the hazardrate
#'    \eqn{\lambda(\tau)} of the *base* distribution \eqn{\Psi}. The hazard
#'    rate of \eqn{\Psi_{t,m}} is simply \eqn{m\lambda(\tau)}, see Details.
#'   
#' * `time_survivalprobability(timedistribution, tau, t, m)`. Evaluates the
#'    survivalfunction \eqn{\Psi_{t,m}} at points `tau`.
#'   
#' * `time_survivalquantile(timedistribution, p, t, m)`. Computes the
#'    \eqn{p}-quantiles of distribution \eqn{\Psi_{t,m}}.
NULL

#' @rdname time_functions
#' @export
time_sample <- function(n, timedistribution, t=0.0, m=1.0) {
  nextnetR_time_sample(as.integer(n), timedistribution, as.double(t), as.double(m))
}

#' @rdname time_functions
#' @export
time_density <- function(timedistribution, tau, t=0.0, m=1.0) {
  tau <- as.double(tau)
  t <- as.double(t)
  m <- as.double(m)
  if (length(t) == 1)
    t <- rep(t, length(tau))
  if (length(m) == 1)
    m <- rep(m, length(tau))
  nextnetR_time_density(timedistribution, tau, t, m)
}

#' @rdname time_functions
#' @export
time_hazardrate <- function(ttimedistributiontr, tau) {
  nextnetR_time_hazardrate(timedistribution, as.double(tau))
}

#' @rdname time_functions
#' @export
time_survivalprobability <- function(timedistribution, tau, t=0.0, m=1.0) {
  tau <- as.double(tau)
  t <- as.double(t)
  m <- as.double(m)
  if (length(t) == 1)
    t <- rep(t, length(tau))
  if (length(m) == 1)
    m <- rep(m, length(tau))
  nextnetR_time_survivalprobability(timedistribution, tau, t, m)
}

#' @rdname time_functions
#' @export
time_survivalquantile <- function(timedistribution, p, t=0.0, m=1.0) {
  p <- as.double(p)
  t <- as.double(t)
  m <- as.double(m)
  if (length(t) == 1)
    t <- rep(t, length(p))
  if (length(m) == 1)
    m <- rep(m, length(p))
  nextnetR_time_survivalquantile(timedistribution, p, m, t)
}

