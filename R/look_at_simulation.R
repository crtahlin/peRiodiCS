#' @title Prints ot results for a particular simulation
#' 
#' @description Prints out results for a particular simulation
#'  (if they exist) saved ina a list.
#'  
#' @export
look_at_simulation <- function(par4sin, par5sin, results_list, file = NULL) {
  
  hash <- digest(paste0(par4sin,",", par5sin)) # calculate hash value to find the result
  results <- results_list[[hash]] # extract results from a list
  
  # calculate average from data saved in my.res
  averages <- apply(X= results$my.res, FUN = function(x) {mean(x)}, MARGIN = 2)
  
  # make the averages into a data frame with named rows/columns
  attach(as.list(averages))
  averages_dt <- data.frame(
    rcs = c(brier.rcs.train, brier.rcs, AUCTrainEst.rcs, AUCNewEst.rcs, cal.rcs.1, cal.rcs.2, lrt.p.value.rcs.train, score.p.value.rcs.train),
    rcs.per = c(brier.rcs.per.train, brier.rcs.per, AUCTrainEst.rcs.per, AUCNewEst.rcs.per, cal.rcs.per.1, cal.rcs.per.2, lrt.p.value.rcs.per.train, score.p.value.rcs.per.train),
    cs.per = c(brier.cs.per.train, brier.cs.per, AUCTrainEst.cs.per, AUCNewEst.cs.per, cal.cs.per.1, cal.cs.per.2, lrt.p.value.cs.per.train, score.p.value.cs.per.train),
    row.names = c("Brier - train", "Brier - test", "AUC - train", "AUC - test", "Calibration intercept", "Calibration slope", "LRT P value", "Score P value")
  ) 
  
  # plot sine curve used as the originial probability function
  curve(sine_function(x, par4sin = results$par4sin, par5sin = results$par5sin), from = 0, to = 1, main = "Probability function used in simulation")
  # print out the averages
  print(averages_dt)
  # print out used parameters of the sine function (TODO: add the rest)
  print(paste("par4sin = ",results$par4sin,"; par5sin =", results$par5sin))
}