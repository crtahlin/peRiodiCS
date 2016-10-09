#' @title Sine function used in the simulations (for simulating probability)
#' 
#' @export
sine_function <- function(x,
                          par1sin = 1,   par2sin=0.25,    par3sin=0.25,  par4sin=1,        par5sin=0,
                          par1trend = 1, par2trend = 0.1, par3trend = 0, par4trend = 1/10, par5trend = 0,
                          add_trend = TRUE,
                          maxvalue = 0.95,
                          minvalue = 0.05) {
  result <- ((par1sin   + sin(x * 2 * pi * par4sin   + (par5sin * 2 * pi)   )) * par2sin    + par3sin ) +
            add_trend *
            ((par1trend + sin(x * 2 * pi * par4trend + (par5trend * 2 * pi) )) * par2trend  + par3trend)
  # keep results inside som threshold (probability has to be between [0,1])
  result[result > maxvalue] <- maxvalue
  result[result < minvalue] <- minvalue
  return(result)
}