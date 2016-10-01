#' @title Sine function used in the simulations (for simulating probability)
#' 
#' @export
sine_function <- function(par1sin, par2sin, par3sin, par4sin, par5sin) {
  result <- (par1sin + sin(x * 2 * pi * par4sin + (par5sin * 2 * pi) )) * par2sin + par3sin 
  return(result)
}