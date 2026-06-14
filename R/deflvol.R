#' Corneal deflection volume at HC
#'
#' Computes corneal deflection volume (DeflVol) at the highest concavity (HC)
#' moment using a two-sphere intersection model. The default method is the
#' closed-form analytical volume; an equivalent disk-integration method is
#' also available.
#'
#' @param CR Numeric vector of corneal curvature radii.
#' @param IR Numeric vector of inverse curvature radii.
#' @param defl_amp Numeric vector of deflection amplitudes at HC. Equivalent
#'   to `DeflAmp` or `HCDeflAmp` in legacy scripts.
#' @param method `"closed"` (default) or `"integral"`.
#' @param N Number of trapezoid segments for the integral method (default 4000).
#' @return Numeric vector of deflection volumes. If inputs are in mm, volumes
#'   are in μL (1 mm³ = 1 μL). Invalid rows return `NA`; non-overlapping
#'   geometry returns 0.
#' @export
#' @examples
#' deflvol(7.8, 6.4, 0.93)
#' deflvol(c(7.8, 7.5), c(6.4, 6.2), c(0.93, 0.90))
#' # with dplyr: data |> dplyr::mutate(DeflVol = deflvol(CR, IR, DeflAmp))
deflvol <- function(
  CR,
  IR,
  defl_amp,
  method = c("closed", "integral"),
  N = 4000L
) {
  method <- match.arg(method)
  if (method == "closed") {
    mapply(
      compute_deflvol_closed_one,
      CR, IR, defl_amp,
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  } else {
    mapply(
      compute_deflvol_integral_one,
      CR, IR, defl_amp,
      MoreArgs = list(N = N),
      SIMPLIFY = TRUE,
      USE.NAMES = FALSE
    )
  }
}
