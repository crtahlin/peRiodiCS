# Cyclic cubic spline using truncated power basis
# Author: Crt Ahlin, crt.ahlin@numbersinlife.com
# Date: 2014

#' @title Function to calculate coefficients for a cyclic cubic spline using 
#' truncated power basis.
#' @description
#' The coefficients for a cyclic function constructed using cubic splines
#' with truncated power basis.
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