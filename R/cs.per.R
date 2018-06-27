#' @title Periodic cubic splines design matrix
#' 
#' @description implementation of the formulas of zhang for any support of x
#' 
#' @param x numeric vector
#' @param knots vector with knot locations
#' @param nk number of knots, overridden if knots are provided
#' @param xmax value of the end of the cycle
#' @param mmin value of the beginning of the cycle
#' 
#' @export
cs.per=function(x, knots = NULL,
                nk = 5,
                xmax = max(x, na.rm=TRUE),
                xmin = min(x, na.rm=TRUE)){
# browser()
  # if they are not provided, derive the knots using function from the rms package
  if( is.null(knots) ) {
    knots <- rcspline.eval(x, nk = nk, knots.only = TRUE)
  }
  
  nk <- length(knots) # set number of knots
  
  rcs.out <- matrix(NA, ncol = nk, nrow = length(x)) # prepare an empty design matrix
  
  # do a loop across all columns of the design matrix, calculating the columns
  for( j in 1:nk ){
    
    a_j <- (-1 / (xmax - xmin)) * (  
      ((xmax^2 + xmin^2 + 4 * xmin * xmax) / 2)  * (xmax - knots[j]) -
        ((3 * (xmax + xmin) / 2) * (xmax - knots[j])^2) + (xmax - knots[j])^3 ) # a_j
  
    b_j <- ((3 * (xmax + xmin) * (xmax - knots[j])) / (2*(xmax - xmin))) - 
      ((3 * (xmax - knots[j])^2) / (2 * (xmax - xmin))) # b_j
    
    c_j <- ( - (xmax - knots[j])/(xmax - xmin)) # c_j 
    # c_j <- ( (xmax - knots[j])/(xmax - xmin)) # c_j # test without a minus sign
    
    
    rcs.out[, j] <- (a_j * x) + (b_j * (x^2)) + (c_j * (x^3)) + ifelse( (x - knots[j] > 0), (x - knots[j])^3, 0) # value of column

  }
  
  # saves the knot locations as attributes
  attr(rcs.out, "knots") <- knots
  
  # return result
  return(rcs.out)
  
}#end b.rcs.cubic.2c