validate_inputs <- function(CR, IR, defl_amp) {
  if (any(is.na(c(CR, IR, defl_amp)))) {
    return(FALSE)
  }
  is.finite(CR) && is.finite(IR) && is.finite(defl_amp) && CR > 0 && IR > 0
}
