#' @title Function to calculate coefficients for a cyclic fractional polynomial.
#' 
#' @description
#' Returns the  coefficients for a cyclic function constructed using cyclic 
#' fractional polynomial.
#' 
#' @param data The data to estimate on. Assumes a column named CyclicTime.
#' 
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