#' @title Derive first derivatives of RCS
#' 
#' @description 
#' function that derives the first derivative of the restricted cubic splines 
#' for a value/vector of values, given the knots
#' 
#' @param x vector of values
#' @param knots vector of knot locations
#' 
#' @export
b_rcs_prime <- function (x, knots) {
  num.knots <- length(knots)
  tk <- knots[num.knots]
  tkmin1 <- knots[num.knots - 1]

  res <- lapply( (1:(num.knots - 2) ),
                function(i){
                  tj <- knots[i]  
                  ifelse((x - tj)^3 > 0, 3*(x - tj)^2, 0) -
                    ifelse((x - tkmin1)^3 > 0, 3*(x - tkmin1)^2*(tk - tj)/(tk - tkmin1), 0) +
                    ifelse((x - tk)^3 > 0, 3*(x - tk)^2*(tkmin1 - tj)/(tk - tkmin1), 0)
                })
  
  res <- matrix(unlist(res), ncol = num.knots - 2)
  
  # return result
  return(res)
  
} # end b_rcs_prime