compute_deflvol_dual_one <- function(CR, IR, defl_amp, N = DEFAULT_N) {
  deflvol_integral <- compute_deflvol_integral_one(CR, IR, defl_amp, N = N)
  deflvol_closed <- compute_deflvol_closed_one(CR, IR, defl_amp)
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
