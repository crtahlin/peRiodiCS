# Cyclic cubic spline using truncated power basis
# Author: Crt Ahlin, crt.ahlin@numbersinlife.com
# Date: 2014

#' @title Function to calculate coefficients for a cyclic cubic spline using truncated power basis.
#' @description
#' The coefficients for a cyclic function constructed using cubic splines with truncated power basis.
#' 
#' @export
cyclicCubicSplineTPB <- function (data, coef=FALSE)
{
  # presume cyclic time data with domain [0,1]
  # presume 9 equaly spaced knots (10 segments)
  knots <- seq(0.1, 0.9, by=0.1)
  
  # presume a column named CyclicTime in data contains the time data with
  # domain [0,1] (the x values)
  # load time data into separate variable, for clarity
  x <- data$CyclicTime
  
  # construct coeficients of the design matrix (beta[1:3] not needed as they are dependent on others; see article)
  const <- rep(x=1,times=dim(data)[1])
  beta4 <- tpower(x, t=knots[1], 3) - 0.9*x^3 - 3/2*((0.9)^2 - (0.9))*x^2 + (-(0.9)^3 + 3/2*(0.9)^2 - 1/2*(0.9))*x
  beta5 <- tpower(x, t=knots[2], 3) - 0.8*x^3 - 3/2*((0.8)^2 - (0.8))*x^2 + (-(0.8)^3 + 3/2*(0.8)^2 - 1/2*(0.8))*x
  beta6 <- tpower(x, t=knots[3], 3) - 0.7*x^3 - 3/2*((0.7)^2 - (0.7))*x^2 + (-(0.7)^3 + 3/2*(0.7)^2 - 1/2*(0.7))*x
  beta7 <- tpower(x, t=knots[4], 3) - 0.6*x^3 - 3/2*((0.6)^2 - (0.6))*x^2 + (-(0.6)^3 + 3/2*(0.6)^2 - 1/2*(0.6))*x
  beta8 <- tpower(x, t=knots[5], 3) - 0.5*x^3 - 3/2*((0.5)^2 - (0.5))*x^2 + (-(0.5)^3 + 3/2*(0.5)^2 - 1/2*(0.5))*x
  beta9 <- tpower(x, t=knots[6], 3) - 0.4*x^3 - 3/2*((0.4)^2 - (0.4))*x^2 + (-(0.4)^3 + 3/2*(0.4)^2 - 1/2*(0.4))*x
  beta10 <- tpower(x, t=knots[7], 3) - 0.3*x^3 - 3/2*((0.3)^2 - (0.3))*x^2 + (-(0.3)^3 + 3/2*(0.3)^2 - 1/2*(0.3))*x
  beta11 <- tpower(x, t=knots[8], 3) - 0.2*x^3 - 3/2*((0.2)^2 - (0.2))*x^2 + (-(0.2)^3 + 3/2*(0.2)^2 - 1/2*(0.2))*x
  beta12 <- tpower(x, t=knots[9], 3) - 0.1*x^3 - 3/2*((0.1)^2 - (0.1))*x^2 + (-(0.1)^3 + 3/2*(0.1)^2 - 1/2*(0.1))*x
  
  # construct the design matrix from the coeficients
  designMatrix <- cbind(const, beta4, beta5, beta6, beta7, beta8, beta9, beta10, beta11, beta12)
  
  # solve the least squares system
  # assuming the data has a column named Value containing the y values
  solution <- solve(t(designMatrix)%*%designMatrix) %*% t(designMatrix) %*% data$Value
  
  # calculate the dependent coefficients
  beta3 <- - tempsum(solution[-1,], knots, 1)
  beta2 <- (3/2) * ( - tempsum(solution[-1,], knots, 2) - beta3 )
  beta1 <- - beta2 - beta3 - tempsum(solution[-1,], knots, 3)
  
  # merge coefficients into final solution (inserting beta[1:3] between the rest of them)
  finalSolution <- rbind(solution[1, , drop=FALSE], beta1, beta2, beta3, solution[-1, , drop=FALSE] )
  
  if (coef==FALSE) {
    # return the final solution
    return(finalSolution)
  }
  # NOTE: R seems to only support matrix size (element count) of up tu 2^31-1 (=2147483647)
  # See: https://stat.ethz.ch/pipermail/r-help/2007-June/133238.html
  # So the largest number of observations supported (with 10 coeeficients to estimate) would be
  # sqrt(2147483647) = 46340.95 or 46340? Although it seems to break at 20.000 already...
  # Perhaps the bigmemory package would solve the problem? 
  # See: http://cran.r-project.org/web/packages/bigmemory/bigmemory.pdf
  if (coef==TRUE) {
    return(designMatrix %*% solve(t(designMatrix)%*%designMatrix) %*% t(designMatrix))
  }
  
}

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

######## Helper functions ##########

# truncated p-th power function
tpower <- function(x, t, p) {
  (x - t) ^ p * (x > t)}

# function to help in calculating the fixed coefficients
tempsum <- function (solution, knots, power) {
  sum <- 0
  for (i in 1:length(knots)) {
    sum <- sum + solution[i]*((1 - knots[i])^power)
  }
  return(unname(sum))
}


