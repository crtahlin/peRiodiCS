#' @title Generate design matrix for periodic cubic splines
#' 
#' @description 
#' Generate design matrix for periodic cubic splines.
#' 
#' @param x numerical x values to transform to new basis
#' @param knots vector with locations of the knots of the spline
#' @param nk number of knots, used only if the knots are not specified, overridden otherwise
#' @param xmax value of the (theoretical) minimum of x
#' @param xmin value of the (theoretical) maximum of x
#' 
#' @import Hmisc
#' @importFrom Hmisc rcspline.eval
#' 
#' @examples 
#' # load example data; see help("viral_east_mediteranean")
#' data("viral_east_mediteranean")
#' 
#' # calculate location of knots to use
#' Knots <- 
#'  Hmisc::rcspline.eval(x = viral_east_mediteranean$EpiWeek,
#'                       nk = 5, knots.only = TRUE)
#'
#' # model viral infections vs weeks
#' model <- glm(RSV ~ cs_per(EpiWeek, knots = Knots), data = viral_east_mediteranean)
#'
#' # plot model (with many points, to make it smooth)
#' plot_per_mod(Model = model, XvarName = "EpiWeek")
#' 
#' @export
cs_per <- function(x,
                   knots = NULL,
                   nk = 5,          
                   xmax = max(x, na.rm=TRUE),
                   xmin = min(x, na.rm=TRUE)){
  
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
  
}# end cs.per