#' Synthetic toy cortisol data
#'
#' A small, fully synthetic data set generated from the NMEA curve model. It is
#' intended for examples, tests, and vignettes; it contains no real participant
#' data.
#'
#' @format A data frame with 240 rows and 5 variables:
#' \describe{
#'   \item{Subject}{Synthetic subject identifier.}
#'   \item{Time}{Hours since waking.}
#'   \item{Truth}{Noise-free NMEA curve value.}
#'   \item{Sigma_t}{Observation-level noise standard deviation.}
#'   \item{Conc_Obs}{Noisy synthetic cortisol value.}
#' }
"toy_cortisol"
