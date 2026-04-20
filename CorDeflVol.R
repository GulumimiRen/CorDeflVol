# ==============================================================
# DeflVol calculation
# The input data.csv must contain columns:
#   - CR: cornea sphere radius
#   - IR: inverse sphere radius
#   - DeflAmp: deflection amplitude parameter (HC moment)
# Ensure consistent units among CR, IR, and DeflAmp (e.g., mm).
# Output:
#   - df$DeflVol (volume, mm^3 if inputs are mm)
# ==============================================================

library(readr)
library(dplyr)

df <- read_csv("data.csv")

compute_deflvol_one <- function(CR, IR, DeflAmp, N = 4000) {
  if (any(is.na(c(CR, IR, DeflAmp)))) return(NA_real_)
  if (CR <= 0 || IR <= 0) return(NA_real_)

  # Sphere-center distance:
  # d = z_cornea - z_inverse = CR + IR - DeflAmp
  d <- CR + IR - DeflAmp

  # Choose a convenient coordinate system:
  # z_inverse = 0, z_cornea = d
  z_inv <- 0
  z_cor <- d

  # Integration bounds from the overlap of the two spheres' valid ranges
  z_min <- max(z_cor - CR, z_inv - IR)
  z_max <- min(z_cor + CR, z_inv + IR)

  if (!is.finite(z_min) || !is.finite(z_max) || z_min >= z_max) return(0)

  # Disk integration (numerical quadrature via trapezoid rule)
  z <- seq(z_min, z_max, length.out = N)

  # Cross-sectional radii with numerical safety (clamp negatives to 0)
  rho_cor2 <- CR^2 - (z - z_cor)^2
  rho_inv2 <- IR^2 - (z - z_inv)^2
  rho_cor <- sqrt(pmax(rho_cor2, 0))
  rho_inv <- sqrt(pmax(rho_inv2, 0))

  # Overlapping-region effective radius
  rho <- pmin(rho_cor, rho_inv)
  y <- rho^2

  # Integral of y dz using the trapezoidal rule
  dz <- diff(z)
  integral_y <- sum((y[-1] + y[-length(y)]) / 2 * dz)

  # DeflVol = pi * integral[rho(z)^2 dz]
  DeflVol <- pi * integral_y
  return(DeflVol)
}

# Compute DeflVol for all rows
df$DeflVol <- mapply(
  compute_deflvol_one,
  df$CR, df$IR, df$DeflAmp,
  MoreArgs = list(N = 4000)
)

# Optional: save results
# write_csv(df, "data_with_DeflVol.csv")
