DEFAULT_N <- 4000L

compute_deflvol_integral_one <- function(CR, IR, defl_amp, N = DEFAULT_N) {
  if (!validate_inputs(CR, IR, defl_amp)) {
    return(NA_real_)
  }
  if (!is.finite(N) || N < 2) {
    stop("N must be an integer >= 2.")
  }

  d <- CR + IR - defl_amp

  z_inv <- 0
  z_cor <- d

  z_min <- max(z_cor - CR, z_inv - IR)
  z_max <- min(z_cor + CR, z_inv + IR)

  if (!is.finite(z_min) || !is.finite(z_max) || z_min >= z_max) {
    return(0)
  }

  z <- seq(z_min, z_max, length.out = as.integer(N))

  rho_cor2 <- CR^2 - (z - z_cor)^2
  rho_inv2 <- IR^2 - (z - z_inv)^2

  rho_sq <- pmin(pmax(rho_cor2, 0), pmax(rho_inv2, 0))

  pi * sum(diff(z) * (head(rho_sq, -1) + tail(rho_sq, -1)) / 2)
}
