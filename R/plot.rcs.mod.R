#' @title Just some ploting function
#' 
#' @description Plots graph in a certain, nice looking way ;)
#' 
#' @export
Plot.rcs.mod<-function(My.mod,
                       my.var,
                       my.ylab="Response",
                       my.xlab="Covariate",
                       my.xlim=NULL,
                       my.ylim=NULL,
                       my.knots=NULL,
                       my.title=NULL,
                       vlines=NULL,
                       hlines=NULL,
                       my.cex.lab=NULL,
                       my.axes=TRUE,
                       my.cex.main=NULL,
                       my.cex.axis=NULL){
  
  my.intervals<-c(0.0001, 0.001, .01, .1, 1, 10, 100, 1000, 10000, 100000)
  my.range<- diff(range(my.var, na.rm=T))
  
  if(is.null(vlines)){
    my.by<-my.intervals[which(my.range/my.intervals<10)[1]]
    my.vlines=seq(-100000, 10000, by=my.by)
    
  } 
  
  my.pred<-predict(My.mod, type="response", se=T)
  
  my.range.y=diff(range(my.pred, na.rm=T))
  
  if(is.null(hlines)){
    my.by.h<-my.intervals[which(my.range/my.intervals<10)[1]]
    my.hlines=seq(-100000, 10000, by=my.by.h)
    
  } 
  
  
  matplot( rbind(my.var, my.var, my.var), rbind(my.pred$fit, my.pred$fit-1.96*my.pred$se, my.pred$fit+1.96*my.pred$se), pch=rep(1, length(my.pred$fit)*3), col=1, type="n", xlab=my.xlab, ylab=my.ylab, xlim=my.xlim, ylim=my.ylim, main=my.title, cex.lab=my.cex.lab, axes=my.axes, cex.main=my.cex.main, cex.axis=my.cex.axis)
  lines( my.var[order(my.var)], my.pred$fit[order(my.var)])
  lines( my.var[order(my.var)], (my.pred$fit+1.96*my.pred$se)[order(my.var)], lty=2)
  lines( my.var[order(my.var)], (my.pred$fit-1.96*my.pred$se)[order(my.var)],  lty=2)
  
  
  # abline(h=seq(0,1,by=.1), v= vlines  , lty=3, col="light grey")
  if(!is.null(hlines) | !is.null(vlines)) abline(h=hlines, v= vlines  , lty=3, col="light grey") else grid()
  if(!is.null(my.knots)) axis(1, at=my.knots, line=2, cex.axis=.65)
  
}