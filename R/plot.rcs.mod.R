#' @title Ploting function for periodic curves
#' 
#' @description 
#' Plots graph of periodic curves with confidence intervals.
#' 
#' @param Model The built model 
#' 
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
#' @param 
#' 
#' @export
Plot.rcs.mod<-function(Model,
                       my.var,
                       Ylab="Response",
                       Xlab="Covariate",
                       Ylim=NULL,
                       Xlim=NULL,
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
                       predict=FALSE,
                       sequence.var=NULL){
  
  my.intervals<-c(0.0001, 0.001, .01, .1, 1, 10, 100, 1000, 10000, 100000)
  my.range<- diff(range(my.var, na.rm=T))
  
  if(is.null(Vlines)){
    my.by<-my.intervals[which(my.range/my.intervals<10)[1]]
    my.vlines=seq(-100000, 10000, by=my.by)
    
  }
  
  
  my.pred<-predict(Model, type="response", se=T)
  
  my.range.y=diff(range(my.pred, na.rm=T))
  
  if(is.null(Hlines)){
    my.by.h<-my.intervals[which(my.range/my.intervals<10)[1]]
    my.hlines=seq(-100000, 10000, by=my.by.h)
    
  }
  
  if (predict == TRUE) {
    my.pred <- predict(Model, type = "response", se = TRUE, newdata = sequence.var) 
  }
  
  #browser()
  if(Add==FALSE)
    matplot( rbind(my.var, my.var, my.var),
             rbind(my.pred$fit, my.pred$fit-1.96*my.pred$se, my.pred$fit+1.96*my.pred$se),
             pch=rep(1, length(my.pred$fit)*3),
             type="n", xlab=Xlab, ylab=Ylab, xlim=Xlim, ylim=Ylim,
             main=Title, cex.lab=Cex.lab,
             axes=Axes, cex.main=Cex.main, cex.axis=Cex.axis, col=Col)
  #else matlines( rbind(my.var, my.var, my.var), rbind(my.pred$fit, my.pred$fit-1.96*my.pred$se, my.pred$fit+1.96*my.pred$se), pch=rep(1, length(my.pred$fit)*3), col=1, type="n")
  lines( my.var[order(my.var)], my.pred$fit[order(my.var)], col=Col)
  if (PlotCI==TRUE) {lines( my.var[order(my.var)], (my.pred$fit+1.96*my.pred$se)[order(my.var)], lty=2, col=Col)}
  if (PlotCI==TRUE) {lines( my.var[order(my.var)], (my.pred$fit-1.96*my.pred$se)[order(my.var)],  lty=2, col=Col)}
  
  
  # abline(h=seq(0,1,by=.1), v= Vlines  , lty=3, col="light grey")
  if(!is.null(Hlines) | !is.null(Vlines)) abline(h=Hlines, v= Vlines  , lty=3, col="light grey") else grid()
  if(!is.null(Knots)) axis(1, at=Knots, line=2, cex.axis=.65)
  
}