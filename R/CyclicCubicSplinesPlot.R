#' @title Function to plot a curve for a cyclic cubic spline using truncated power basis.
#' 
#' @description
#' function for generating the curve predicted by Cyclic Cubic Spline using 
#' Truncated Power Series 
#' 
#' @param data The data to plot. (Should have at most ~15k observations, 
#' due to memory issues on 32bit systems?)
#' @param curve OBSOLETE TRUE or FALSE. Plots either a curve (TRUE) or 
#' a set of points (FALSE).
#' @param add If TRUE, adds the curve to an existing plot. 
#' If FALSE, it plots a new plot.
#' 
#' @export
cyclicCubicSplineTPBCurve <- function (data,
                                       add = FALSE,
                                       se = FALSE) {
  # real data points
  xReal <- data$CyclicTime
  yReal <- data$Value
  
  # presume 9 equaly spaced knots (10 segments)
  # TODO: use the sam knots in cyclicCubicSplineTPB()
  knots <- seq(0.1, 0.9, by=0.1)
  
##############################################################################
# helper functions for calculating the predicted x values basis

  # construct matrix of values to predict from
  constructXcoordinates <- function (x) {
    dummyData <- matrix(nrow=length(x), ncol=13) # matrix of data 
    for(i in 1:length(x)) {
      dummyData[i,1] <- 1
      dummyData[i,2] <- x[i]
      dummyData[i,3] <- x[i]^2
      dummyData[i,4] <- x[i]^3
      dummyData[i,5] <- tpower(x=x[i], t=knots[1], 3)
      dummyData[i,6] <- tpower(x=x[i], t=knots[2], 3)
      dummyData[i,7] <- tpower(x=x[i], t=knots[3], 3)
      dummyData[i,8] <- tpower(x=x[i], t=knots[4], 3)
      dummyData[i,9] <- tpower(x=x[i], t=knots[5], 3)
      dummyData[i,10] <- tpower(x=x[i], t=knots[6], 3)
      dummyData[i,11] <- tpower(x=x[i], t=knots[7], 3)
      dummyData[i,12] <- tpower(x=x[i], t=knots[8], 3)
      dummyData[i,13] <- tpower(x=x[i], t=knots[9], 3)
    } 
    return(dummyData)
  }
  
  # function to predict values at x based on parameters estimated from real data
  prediction <- 
    function (x, parameters=cyclicCubicSplineTPB(data)) {
    constructXcoordinates(x)%*%parameters
    }
  
  # id add is FALSE, a new plot should be generated 
  if (add==FALSE) {plot(data$Value ~ data$CyclicTime, ylim=c(30,200))}
  
  # add a curve predicted for x from 0 to 1
  curve(prediction, from=0, to=1, add=TRUE, lwd=2, col="red")
   
  # add std. error bands if se = TRUE
  if (se == TRUE) {
    # TODO : check logic bellow
    # variance of residuals around the graph; residual / degrees of freedom
    sigmaSquared <- sum((data$CyclicTime - prediction(data$CyclicTime))^2)/
      (dim(data)[1]-10)  # degrees of freedom
    # pointwise standard error
    se <- sqrt(sigmaSquared * diag(cyclicCubicSplineTPB(data, coef=TRUE)))
    # width of the error bands (t distribution, 95% confidence)
    bandWidth <- qt(.975,(dim(data)[1]-10))*se
    # plot the error bands - order the x values first, to make lines visible
    sortedTime <- sort(data$CyclicTime)
    points(sortedTime, prediction(sortedTime) + bandWidth, col="blue", type="l")
    points(sortedTime, prediction(sortedTime) - bandWidth, col="blue", type="l")
    axis(side=4)
  }
    
}