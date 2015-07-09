#' @title Function for generating the curve estimated by cyclic fractional polynomial.
#' 
#' @param data Data to estimate on.
#' @param curve If TRUE, a curve is plotted instead of points.
#' @param add If TRUE, the plot is added to an existing plot.
#' 
#' @export
cyclicFractionalPolynomialCurve <- function (data,
                                             curve=TRUE,
                                             add=FALSE) {
  
  # real data points
  xReal <- data$CyclicTime
  yReal <- data$Value
  
  # predicted data points
  x <- seq(0,1, by=0.001)
  testData <- matrix(nrow=1001, ncol=7)
  for(i in 1:length(x)) {
    testData[i,1] <- 1
    testData[i,2] <- x[i]
    testData[i,3] <- x[i]^2
    testData[i,4] <- x[i]^3
    testData[i,5] <- (x[i]+0.5)^(-2)
    testData[i,6] <- (x[i]+0.5)^(-1)
    testData[i,7] <- (x[i]+0.5)^(-1/2)
  }
  
  
  # helper function for calculating the predicted x values basis
  constructXcoordinates <- function (x) {
    testData <- matrix(nrow=length(x), ncol=7)
    for(i in 1:length(x)) {
      testData[i,1] <- 1
      testData[i,2] <- x[i]
      testData[i,3] <- x[i]^2
      testData[i,4] <- x[i]^3
      testData[i,5] <- (x[i]+0.5)^(-2)
      testData[i,6] <- (x[i]+0.5)^(-1)
      testData[i,7] <- (x[i]+0.5)^(-1/2)
    } 
    return(testData)
  }
  
  # function to calculate predicte y from x
  prediction <- function (x, parameters=cyclicFractionalPolynomial(data)) {constructXcoordinates(x)%*%parameters}
  
  # plot the real points
  if (add==FALSE) {plot(data$Value ~ data$CyclicTime)}
  #browser()
  if (curve==TRUE) {
    curve(prediction, from=0, to=1, add=TRUE, lwd=4, col="blue")
  } else {
    yPredicted <- testData %*% cyclicFractionalPolynomial(data)
    points(yPredicted~x, col="blue")}
  
}