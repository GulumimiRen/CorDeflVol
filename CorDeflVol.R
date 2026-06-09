# ==============================================================
# CorDeflVol: DeflVol calculation with two equivalent algorithms
#
# The input data.csv must contain columns:
#   - CR: cornea sphere radius
#   - IR: inverse sphere radius
#   - DeflAmp: deflection amplitude at the HC moment
#
# Ensure consistent units among CR, IR, and DeflAmp (e.g., mm).
# If inputs are in mm, DeflVol is returned in mm^3 (= microlitre).
#
# Output columns:
#   - df$DeflVol_integral: numerical disk integration
#   - df$DeflVol_closed:   closed-form two-sphere intersection
#   - df$DeflVol:          default reported value, equal to closed form
#   - df$DeflVol_abs_error
#   - df$DeflVol_rel_error
# ==============================================================

library(readr)

DEFAULT_N <- 4000L

validate_inputs <- function(CR, IR, DeflAmp) {
  if (any(is.na(c(CR, IR, DeflAmp)))) {
    return(FALSE)
  }
  is.finite(CR) && is.finite(IR) && is.finite(DeflAmp) && CR > 0 && IR > 0
}

# Numerical disk integration: DeflVol = pi * integral[rho(z)^2 dz].
compute_deflvol_integral_one <- function(CR, IR, DeflAmp, N = DEFAULT_N) {
  if (!validate_inputs(CR, IR, DeflAmp)) {
    return(NA_real_)
  }
  if (!is.finite(N) || N < 2) {
    stop("N must be an integer >= 2.")
  }

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

  if (!is.finite(z_min) || !is.finite(z_max) || z_min >= z_max) {
    return(0)
  }

  # Disk integration (numerical quadrature via trapezoid rule)
  z <- seq(z_min, z_max, length.out = as.integer(N))

  # Cross-sectional radii squared, with numerical safety (clamp negatives to 0)
  rho_cor2 <- CR^2 - (z - z_cor)^2
  rho_inv2 <- IR^2 - (z - z_inv)^2

  # Overlap effective radius squared: min(rho_cor^2, rho_inv^2) equals min(rho_cor, rho_inv)^2
  rho_sq <- pmin(pmax(rho_cor2, 0), pmax(rho_inv2, 0))

  # DeflVol = pi * integral[rho(z)^2 dz] via the trapezoidal rule
  pi * sum(diff(z) * (head(rho_sq, -1) + tail(rho_sq, -1)) / 2)
}

# Closed-form volume of intersection of two spheres.
compute_deflvol_closed_one <- function(CR, IR, DeflAmp) {
  if (!validate_inputs(CR, IR, DeflAmp)) {
    return(NA_real_)
  }

  # Sphere-center distance: d = CR + IR - DeflAmp
  d <- CR + IR - DeflAmp

  # Invalid geometry: non-finite or non-positive center distance
  if (!is.finite(d) || d <= 0) {
    return(NA_real_)
  }

  # Separated spheres: no overlap volume
  if (d >= CR + IR) {
    return(0)
  }

  # One sphere fully contained in the other: volume of the smaller sphere
  if (d <= abs(CR - IR)) {
    return(4 * pi * min(CR, IR)^3 / 3)
  }

  # Partial overlap: standard two-sphere intersection formula
  pi * (CR + IR - d)^2 *
    (d^2 + 2 * d * (CR + IR) - 3 * (CR - IR)^2) / (12 * d)
}

compute_deflvol_one <- function(
  CR,
  IR,
  DeflAmp,
  method = c("closed", "integral"),
  N = DEFAULT_N
) {
  method <- match.arg(method)
  if (method == "closed") {
    compute_deflvol_closed_one(CR, IR, DeflAmp)
  } else {
    compute_deflvol_integral_one(CR, IR, DeflAmp, N = N)
  }
}

compute_deflvol_dual_one <- function(CR, IR, DeflAmp, N = DEFAULT_N) {
  deflvol_integral <- compute_deflvol_integral_one(CR, IR, DeflAmp, N = N)
  deflvol_closed <- compute_deflvol_closed_one(CR, IR, DeflAmp)
  abs_error <- abs(deflvol_integral - deflvol_closed)
  rel_error <- abs_error / deflvol_closed

  c(
    DeflVol_integral = deflvol_integral,
    DeflVol_closed = deflvol_closed,
    DeflVol_abs_error = abs_error,
    DeflVol_rel_error = rel_error
  )
}

compute_deflvol_dual <- function(df, N = DEFAULT_N) {
  required <- c("CR", "IR", "DeflAmp")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required column(s): ", paste(missing, collapse = ", "))
  }

  values <- t(mapply(
    compute_deflvol_dual_one,
    df$CR, df$IR, df$DeflAmp,
    MoreArgs = list(N = N)
  ))
  values <- as.data.frame(values)

  df$DeflVol_integral <- values$DeflVol_integral
  df$DeflVol_closed <- values$DeflVol_closed
  df$DeflVol <- df$DeflVol_closed
  df$DeflVol_abs_error <- values$DeflVol_abs_error
  df$DeflVol_rel_error <- values$DeflVol_rel_error
  df
}

if (file.exists("data.csv")) {
  df <- read_csv("data.csv", show_col_types = FALSE)
  df <- compute_deflvol_dual(df, N = DEFAULT_N)
} else {
  message("No data.csv found; DeflVol functions are loaded but df was not created.")
}

# Optional: save results
# write_csv(df, "data_with_DeflVol.csv")
