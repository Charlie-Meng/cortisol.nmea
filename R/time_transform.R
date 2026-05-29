#' Forward time transformation for asymmetric cortisol curves
#'
#' Applies the NMEA pre-peak time warp used to encode asymmetric cortisol
#' trajectories. Times at or after `mu` are left unchanged.
#'
#' @param time Numeric vector of times since waking.
#' @param alpha Numeric asymmetry parameter. Values below 1 stretch the
#'   pre-peak interval; values above 1 compress it.
#' @param mu Numeric peak-location parameter.
#'
#' @return A numeric vector on the transformed time scale.
#' @export
time_forward <- function(time, alpha, mu) {
  time_new <- alpha * time + (1 - alpha) * mu
  nochange <- which(time >= mu)
  time_new[nochange] <- time[nochange]
  time_new
}

#' Backward time transformation for asymmetric cortisol curves
#'
#' Inverts [time_forward()] for times before the peak and leaves post-peak
#' times unchanged.
#'
#' @param time_new Numeric vector on the transformed time scale.
#' @param alpha Numeric asymmetry parameter.
#' @param mu Numeric peak-location parameter.
#'
#' @return A numeric vector on the original time scale.
#' @export
time_backward <- function(time_new, alpha, mu) {
  time <- (time_new - (1 - alpha) * mu) / alpha
  nochange <- which(time_new >= mu)
  time[nochange] <- time_new[nochange]
  time
}
