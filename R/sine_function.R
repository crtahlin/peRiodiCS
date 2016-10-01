#' @title Sine function used in the simulations (for simulating probability)
#' 
#' @export
sine_function <- function(x, par1sin=1, par2sin=0.25, par3sin=0.25, par4sin=1, par5sin=0) {
  result <- (par1sin + sin(x * 2 * pi * par4sin + (par5sin * 2 * pi) )) * par2sin + par3sin 
  return(result)
}