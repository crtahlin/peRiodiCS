#' @title Ploting function for periodic curves model
#' 
#' @description 
#' Plots graph of periodic curves with confidence intervals.
#' Data should be included in the model.
#' 
#' @param Model The built model 
#' @param Xvarname Name of the x variable in the dataset (column name)
#' @param Ylab Label on vertical (y) axis
#' @param Xlab Label on horizontal (x) axis
#' @param Xlim Limits of x axis
#' @param Ylim Limits of y axis
#' @param Knots Locations of knots of the splines
#' @param Title Title of the plot
#' @param Vlines Where to plot vertical lines
#' @param Hlines Where to plot horizontal lines
#' @param Cex.lab Character expansion (aka "size of font") for the labels
#' @param Cex.main Character expansion for main text
#' @param Cex.axis Character expansion for the axis text
#' @param Axes Plot axes
#' @param Add Add to existing plot
#' @param Col Color of the plotted lines
#' @param PlotCI Plot confidence intervals
#' @param Smooth Make the Xaxis values equidistant (and the curve smoother)
#' @param Xmin The min X of data to be predicted (if Smooth)
#' @param Xmax The max X of data to be predicted (if Smooth)
#' @param xLocation If smooth FALSE, the location of the x term in model$x[, xLocation]
#' 
#' @export
Plot.per.mod <- function(Model,
                         XvarName,
                         # my.var=NULL,
                         Ylab="Response",
                         Xlab="Covariate",
                         Ylim=NULL,
                         Xlim=NULL,
                         Xmin=NULL,
                         Xmax=NULL,
                         Knots=NULL,
                         Title=NULL,
                         Vlines=NULL,
                         Hlines=NULL,
                         Cex.lab=NULL,
                         Cex.main=NULL,
                         Cex.axis=NULL, 
                         Axes=TRUE,
                         Add=FALSE,
                         Col="black",
                         PlotCI=TRUE,
                         Smooth=FALSE,
                         xLocation=2) {
  
  # NOTE: depends on how the model was built - if NOT by data=some_data, but by direct reference to
  # some_data$variable ~ someothervariable
  # then logic should be different

  # if model does not contain data, stop execution
  stopifnot(!is.null(Model$data))
  
  # load the X variable from the data in the model
  Xvar <- na.omit(Model$data[[XvarName]])
  Xvar <- Xvar[order(Xvar)]
  Xvar <- unique(Xvar)
  
  # determine where to plot vertical lines
  Intervals <- c(0.0001, 0.001, .01, .1, 1, 10, 100, 1000, 10000, 100000)
  my.range <- diff(range(Xvar, na.rm=TRUE))
  if(is.null(Vlines)){ # if vertical line are not user defined
    By.v <- Intervals[which(my.range/Intervals < 10)[1]]
    Vlines <- seq(-100000, 10000, by=By.v) # TODO: do this more generaly?
  }
  
  # determine where to plot horizontal lines
  Prediction <- predict(Model, type="response", se=TRUE) # determine possible y
  my.range.y <- diff(range(Prediction, na.rm=T))
  if(is.null(Hlines)){ # if horizontal lines are not defined
    By.h<-Intervals[which(my.range.y/Intervals < 10)[1]]
    Hlines <- seq(-100000, 10000, by=By.h)
  }
  # make new X values to make predictions smoother
  if (Smooth == TRUE) { # if using prediction on equidistant x value interval
    if (is.null(Xmin)) {Xmin <- min(Model$data[[XvarName]], na.rm = TRUE)}
    if (is.null(Xmax)) {Xmax <- max(Model$data[[XvarName]], na.rm = TRUE)}
    Xvar <- seq(Xmin, Xmax, length.out = dim(Model$data)[1] ) # make a sequence of 1000 Xses to plot smoothly
    NewData <- data.frame((Xvar))
    colnames(NewData) <- c(XvarName)
    Prediction <- predict.glm(Model, type = "response", se = TRUE, newdata = NewData) 
  } else {
    # to work with rcs, this has to be done without using new data ...
    Xvar <- Model$x[, xLocation]
    Prediction <- predict(Model, type = "response", se = TRUE) 
    }
  
  # plot main curve
  if(Add==FALSE) {
    matplot( rbind(Xvar, Xvar, Xvar), # make a new plot
             rbind(Prediction$fit, Prediction$fit-1.96*Prediction$se, Prediction$fit+1.96*Prediction$se),
             pch=rep(1, length(Prediction$fit)*3),
             type="n",
             xlab=Xlab, ylab=Ylab, xlim=Xlim, ylim=Ylim,
             main=Title, cex.lab=Cex.lab,
             axes=Axes, cex.main=Cex.main, cex.axis=Cex.axis, col=Col) 
    } 
  # plot the lines in any case
  lines( Xvar[order(Xvar)], Prediction$fit[order(Xvar)], col=Col) # just add line to existing plot
      
  # plot confidence intervals
  if (PlotCI==TRUE) {lines( Xvar[order(Xvar)], (Prediction$fit+1.96*Prediction$se)[order(Xvar)], lty=2, col=Col)}
  if (PlotCI==TRUE) {lines( Xvar[order(Xvar)], (Prediction$fit-1.96*Prediction$se)[order(Xvar)], lty=2, col=Col)}
  
  # plot horizontal and vertical lines
  if(!is.null(Hlines) | !is.null(Vlines)) abline(h=Hlines, v= Vlines, lty=3, col="light grey") else grid()
  if(!is.null(Knots)) axis(1, at=Knots, line=2, cex.axis=.65)
  
}