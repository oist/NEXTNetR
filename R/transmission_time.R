#' @export
time_sample <- function(n, ttr, t=0, m=1) {
  episimR_time_sample(as.integer(n), ttr, as.double(t), as.integer(m))
}

#' @export
time_density <- function(ttr, taus) {
  episimR_time_density(ttr, as.double(taus))
}

#' @export
time_hazardrate <- function(ttr, taus) {
  episimR_time_hazardrate(ttr, as.double(taus))
}

#' @export
time_survivalprobability <- function(ttr, tau, t = 0, m = 1) {
  tau <- as.double(tau)
  t <- as.double(t)
  m <- as.integer(m)
  if (length(t) == 1)
    t <- rep(t, length(tau))
  if (length(m) == 1)
    m <- rep(m, length(tau))
  episimR_time_survivalprobability(ttr, tau, t, m)
}

#' @export
time_survivalquantile <- function(ttr, p, t = 0, m = 1) {
  p <- as.double(p)
  t <- as.double(t)
  m <- as.integer(m)
  if (length(t) == 1)
    t <- rep(t, length(p))
  if (length(m) == 1)
    m <- rep(m, length(p))
  episimR_time_survivalquantile(ttr, p, m, t)
}

#' @export
exponential_time <- function(lambda) {
  episimR_exponential_time(as.double(lambda))
}

#' @export
lognormal_time <- function(mean, var, pinf) {
  episimR_gamma_time(as.double(mean), as.double(var), as.double(pinf))
}

#' @export
gamma_time <- function(mean, var, pinf) {
  episimR_gamma_time(as.double(mean), as.double(var), as.double(pinf))
}

#' @export
generic_time <- function(density, survivalprobability, probability_is_trinary,
                         survivalquantile, quantile_is_trinary,
                         sample, pinfinity)
{
  episimR_generic_time(density, survivalprobability, as.logical(probability_is_trinary),
                       survivalquantile, as.logical(quantile_is_trinary),
                       sample, as.double(pinfinity))
}
