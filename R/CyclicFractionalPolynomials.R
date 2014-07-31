# Cyclic fractional polynomials
# Author: Crt Ahlin, crt.ahlin@numbersinlife.com
# Date: 2014

#' @title Function to calculate coefficients for a cyclic fractional polynomial.
#' @description
#' The coefficients for a cyclic function constructed using cyclic fractional polynomial.
#' @param data The data to estimate on.
#' @export
cyclicFractionalPolynomial <- function(data) {
  # Ta verzija deluje. Očitno je Korenski člen nekaj pokvaril lepo obnašanje -
  # mogoče kaj v zvezi s tem, da postane pri odvodu ^(-1/2), ta člen pa je preveč
  # podoben mnogim ostalim? Pri tej funkciji so tudi enaki samo prvi odvodi,
  # potrebno poskusiti še z drugimi odvodi... Ampak malo kasneje. Sicer tudi tako
  # ne izgleda slabo.
  
  
  # Umaknem člen z ^(1/2), ker se bo mogoče funkcija lepše obnašala
  # f2[x_] := 
  #   b0 + b1*x + b2*x^2 + b3*x^3 + b4*(x + 0.5)^(-2) + 
  #   b5*(x + 0.5)^(-1) + b6*(x + 0.5)^(-1/2) 
  
  # Solve[{f[0] == f[1], f'[0] == f'[1]}, {b0, b1, b2, b3, b4, b5, b6}]
  
  # {{b5 -> 0. - 1.89556 b1 - 3.87973 b2 - 4.87182 b3 - 8.54569 b4, 
  #   b6 -> 0. + 5.90148 b1 + 10.3276 b2 + 12.5406 b3 + 13.1144 b4}}
  
  # creating the design matrix
  # create the coefficients of the design matrix
  # assuming the data have columns named "CyclicTime" and "Values"
  b0 <- rep(1, times=dim(data)[1]) # no restrictions
  b1 <- data$CyclicTime +
    (-1.89556)  * (data$CyclicTime + 0.5)^(-1) +   # b5 restriction
    (5.90148)  * (data$CyclicTime + 0.5)^(-1/2)     # b6 restriction
  b2 <- data$CyclicTime^2 +
    (- 3.87973)  * (data$CyclicTime + 0.5)^(-1) +   # b5 
    (10.3276)  * (data$CyclicTime + 0.5)^(-1/2)     # b6
  b3 <- data$CyclicTime^3 +
    (- 4.87182)  * (data$CyclicTime + 0.5)^(-1) +   # b5 
    (12.5406)  * (data$CyclicTime + 0.5)^(-1/2)     # b6
  b4 <- (data$CyclicTime + 0.5)^(-2) +
    (- 8.54569 )   * (data$CyclicTime + 0.5)^(-1) +   # b65
    (13.1144 ) * (data$CyclicTime + 0.5)^(-1/2)     # b6
  # create the design matrix
  designMatrix <- cbind(b0, b1, b2, b3, b4)
  # solve via least squares
  solution <- solve(t(designMatrix)%*%designMatrix) %*% t(designMatrix) %*% data$Value
  
  # calculate the dependent coeficients
  b5 <- unname((- 1.89556)*solution["b1",] - 
                 3.87973*solution["b2",]  - 
                 4.87182*solution["b3",]  - 
                 8.54569*solution["b4",] ) 
  b6 <- unname((5.90148)*solution["b1",] +
                 10.3276*solution["b2",] + 
                 12.5406*solution["b3",] + 
                 13.1144*solution["b4",] )
  # and create the matrix of the final solution
  finalSolution <- rbind(solution, b5, b6)
  
  # return the final solution
  return(finalSolution)
}


#' @title Function for generating the curve estimated by cyclic fractional polynomial.
#' 
#' @param data Data to estimate on.
#' @param curve If TRUE, a curve is plotted instead of points.
#' @param add If TRUE, the plot is added to an existing plot.
#' @export
cyclicFractionalPolynomialCurve <- function (data, curve=FALSE, add=FALSE) {
  
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