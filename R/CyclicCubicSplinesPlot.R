#' @title Function to plot a curve for a cyclic cubic spline using truncated power basis.
#' 
#' @description
#' function for generating the curve predicted by Cyclic Cubic Spline using Truncated Power Series 
#' @param data The data to plot. (Should have at most ~15k observations, due to memory issues on 32bit systems?)
#' @param curve TRUE or FALSE. Plots either a curve (TRUE) or a set of points (FALSE).
#' @param add If TRUE, adds the curve to an existing plot. If FALSE, it plots a new plot.
#' @export
cyclicCubicSplineTPBCurve <- function (data, curve=FALSE, add=FALSE) {
  # real data points
  xReal <- data$CyclicTime
  yReal <- data$Value
  
  # presume 9 equaly spaced knots (10 segments)
  # TODO: use the sam knots in cyclicCubicSplineTPB()
  knots <- seq(0.1, 0.9, by=0.1)
  
  # predicted data points
  x <- seq(0, 1, by=0.001)
  testData <- matrix(nrow=1001, ncol=13)
  # TODO: rewrite this loop using some version of constructXcoordinates() below?
  for(i in 1:length(x)) {
    testData[i,1] <- 1
    testData[i,2] <- x[i]
    testData[i,3] <- x[i]^2
    testData[i,4] <- x[i]^3
    testData[i,5] <- tpower(x=x[i], t=knots[1], 3)
    testData[i,6] <- tpower(x=x[i], t=knots[2], 3)
    testData[i,7] <- tpower(x=x[i], t=knots[3], 3)
    testData[i,8] <- tpower(x=x[i], t=knots[4], 3)
    testData[i,9] <- tpower(x=x[i], t=knots[5], 3)
    testData[i,10] <- tpower(x=x[i], t=knots[6], 3)
    testData[i,11] <- tpower(x=x[i], t=knots[7], 3)
    testData[i,12] <- tpower(x=x[i], t=knots[8], 3)
    testData[i,13] <- tpower(x=x[i], t=knots[9], 3)
  }
  
  ##############################################################################
  # helper function for calculating the predicted x values basis
  constructXcoordinates <- function (x) {
    testData <- matrix(nrow=length(x), ncol=13)
    for(i in 1:length(x)) {
      testData[i,1] <- 1
      testData[i,2] <- x[i]
      testData[i,3] <- x[i]^2
      testData[i,4] <- x[i]^3
      testData[i,5] <- tpower(x=x[i], t=knots[1], 3)
      testData[i,6] <- tpower(x=x[i], t=knots[2], 3)
      testData[i,7] <- tpower(x=x[i], t=knots[3], 3)
      testData[i,8] <- tpower(x=x[i], t=knots[4], 3)
      testData[i,9] <- tpower(x=x[i], t=knots[5], 3)
      testData[i,10] <- tpower(x=x[i], t=knots[6], 3)
      testData[i,11] <- tpower(x=x[i], t=knots[7], 3)
      testData[i,12] <- tpower(x=x[i], t=knots[8], 3)
      testData[i,13] <- tpower(x=x[i], t=knots[9], 3)
    } 
    return(testData)
  }
  
  prediction <- function (x, parameters=cyclicCubicSplineTPB(data)) {constructXcoordinates(x)%*%parameters}
  
  
  if (add==FALSE) {plot(data$Value ~ data$CyclicTime, ylim=c(30,200))}
  if (curve==TRUE) {
    curve(prediction, from=0, to=1, add=TRUE, lwd=4, col="red")
    
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
  } else {
    yPredicted <- testData %*% cyclicCubicSplineTPB(data)
    points(yPredicted~x, col="blue")}
}