compute_deflvol_closed_one <- function(CR, IR, defl_amp) {
  if (!validate_inputs(CR, IR, defl_amp)) {
    return(NA_real_)
  }

  d <- CR + IR - defl_amp

  if (!is.finite(d) || d <= 0) {
    return(NA_real_)
  }

  if (d >= CR + IR) {
    return(0)
  }

  if (d <= abs(CR - IR)) {
    return(4 * pi * min(CR, IR)^3 / 3)
  }

  pi * (CR + IR - d)^2 *
    (d^2 + 2 * d * (CR + IR) - 3 * (CR - IR)^2) / (12 * d)
}
