#' @title Generate design matrix for periodic restricted cubic spline
#' 
#' @description 
#' Generate design matrix for periodic restricted cubic spline.
#' 
#' @param x numerical x values to transform to new basis
#' @param knots vector with locations of the knots of the spline
#' @param nk number of knots, used only if the knots are not specified, overridden otherwise
#' @param xmin value of the (theoretical) minimum of x
#' @param xmax value of the (theoretical) maximum of x
#' 
#' @export
rcs_per <- function(x,
                    knots = NULL,
                    nk = 5,
                    xmin = min(x, na.rm=TRUE),
                    xmax = max(x, na.rm=TRUE)){
  
  # derive the knot locations as in the rms package, if they are not provided
  if(is.null(knots)) {
    knots <- rcspline.eval(x, nk = nk, knots.only = TRUE)
  }
  
  # check if the number of knots if at least 4
  nk <- length(knots)
  if (nk < 4) stop("To use the periodic RCS you must specify at least 4 knots")
  
  b.x.all <- b_rcs(x, knots) # matrix with the expansion of the splines, n*(k-2)
  b.xmax.all <- b_rcs(xmax, knots) # value of the spline for x=xmax, vector 1*(k-2)
  
  b.prime.x.all <- b_rcs_prime(x, knots) # matrix with the first derivative of the expansion of the splines, n*(k-2)
  b.prime.xmax.all <- b_rcs_prime(xmax, knots) # vector with the value of the first derivative for x=xmax, vector 1*(k-2)
  
  # terms to add in gamma1
  gamma1.c1 <- (xmin-xmax) / b.xmax.all[nk-2] # constant
  gamma.f1 <- gamma1.c1 * b.x.all[, nk-2] # vector n*1
  
  gamma.f2 <- gamma1.c1 * b.prime.xmax.all[nk-2]
  
  denom.beta.kMinus3 <- b.prime.xmax.all[nk-3] - b.xmax.all[nk-3] / b.xmax.all[nk-2] * b.prime.xmax.all[nk-2]
  
  beta.j <- matrix(NA, ncol=nk-3, nrow = length(x))
  
  for(j in 1:(nk-3)) {beta.j[,j] <- b.x.all[,j] - b.xmax.all[j] / b.xmax.all[nk-2] * b.x.all[,nk-2]}
  
  if(nk > 4) {
    beta.c2.j <- matrix(NA, ncol = nk-4, nrow = length(x))
    for(j in 1:(nk-4)) {beta.c2.j[,j] <- b.prime.xmax.all[,j] - b.xmax.all[j] / b.xmax.all[nk-2] * b.prime.xmax.all[nk-2]}
  }
  
  num.beta.kMinus3 <- beta.j[, nk-3]
  
  # construct dataframe with the result
  if (nk > 4) {
    result <- 
      cbind(x+gamma.f1-gamma.f2/denom.beta.kMinus3*num.beta.kMinus3,
               beta.j[,1:(nk-4)]-beta.c2.j[,1:(nk-4)]*num.beta.kMinus3/denom.beta.kMinus3
            )
    } else {
      result <- x+gamma.f1-gamma.f2/denom.beta.kMinus3*num.beta.kMinus3
    }
  
  # saves the knot locations as attributes
  attr(result, "knots") <- knots
  
  # return result
  return(result)
  
} # end rcs_per