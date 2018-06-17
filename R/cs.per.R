#' @title Periodic cubic splines design matrix
#' 
#' @description implementation of the formulas of zhang for any support of x
#' 
#' @param x numeric vector
#' @param knots vector with knot locations
#' @param nk number of knots, overridden if knots are provided
#' @param xmax value of the end of the cycle
#' @param mmin value of the beginning of the cycle
#' 
#' @export
cs.per=function(x, knots=NULL, nk=5, xmax=max(x, na.rm=TRUE), xmin=min(x, na.rm=TRUE)){

  #derive the knots as in the rms package, if not provided
  if(is.null(knots)) {
    knots=rcspline.eval(x, nk=nk, knots.only=TRUE)
  }
  
  nk=length(knots)  
  
  rcs.out=matrix(NA, ncol=nk, nrow=length(x))   
  
  for(j in 1:length(knots)){
    
    my.aj=-1/(xmax-xmin)*(  
      (xmax^2+xmin^2+4*xmin*xmax)/2*(xmax-knots[j])-
        3*(xmax+xmin)/2*(xmax-knots[j])^2+(xmax-knots[j])^3)
  
    my.bj=3*(xmax+xmin)*(xmax-knots[j])/(2*(xmax-xmin))-3*(xmax-knots[j])^2/(2*(xmax-xmin))
    
    my.cj=-(xmax-knots[j])/(xmax-xmin)
    
    rcs.out[,j]=my.aj*x+my.bj*(x^2)+my.cj*(x^3)+ifelse(x-knots[j]>0, (x-knots[j])^3, 0)

  }
  
  #saves the knots
  attr(rcs.out, "knots") <- knots
  
  #return result
  return(rcs.out)
  
}#end b.rcs.cubic.2c