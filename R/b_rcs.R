#' @title Basis for restricted cubic splines
#' 
#' @description 
#' Function that derives the restricted cubic splines 
#' for a value/vector of values, given the knots; 
#' obtains exactly the same results as the rcs function included in the rms package.
#' 
#' @param x numerical vector
#' @param knots vector specifying the knot locations
#' @param inclx logical, if TRUE returns also the x vector
#' 
b_rcs <- function(x, knots, inclx=FALSE) {
  num.knots = length(knots)
  tk = knots[num.knots]
  tkmin1 = knots[num.knots-1]
  
  res <- lapply(1:(num.knots-2), 
                function(i){
                  tj <- knots[i]  
                  pmax((x-tj)^3, 0) -
                    pmax((x-tkmin1)^3, 0)*(tk-tj)/(tk-tkmin1) +
                    pmax((x-tk)^3, 0)*(tkmin1-tj)/(tk-tkmin1)
                })
  
  res <- matrix(unlist(res), ncol=num.knots-2)
  
  if (inclx) { res <- cbind(x, res) }

  # return result
  return(res)
  
} # end b.rcs