#' @title Basis for restricted cubic splines
#' 
#' @description 
#' Function that derives the restricted cubic splines 
#' for a value/vector of values, given the knots; 
#' obtains exactly the same results as the rcs function included in the rms package.
#' 
#' @param x numerical vector
#' @param knots vector specifying the knots
#' @param inclx logical, if TRUE returns also the x vector
#' 
#' @export
b.rcs=function(x, knots, inclx=FALSE) {

  num.knots = length(knots)
  tk = knots[num.knots]
  tkmin1 = knots[num.knots-1]
  
  my.res=lapply(1:(num.knots-2), 
                function(i){
                  tj=knots[i]  
                  pmax((x-tj)^3, 0)-
                    pmax((x-tkmin1)^3,0)*(tk-tj)/(tk-tkmin1)+
                    pmax((x-tk)^3,0)*(tkmin1-tj)/(tk-tkmin1)
                })
  
  my.res=matrix(unlist(my.res), ncol=num.knots-2)
  
  if (inclx) my.res=cbind(x, my.res)

  # return result
  return(my.res)
  
}#end b.rcs


#check: 
#
#x=rnorm(1000)
#library(rms)
#xx = rcspline.eval(x, nk=6, norm=0, inclx=FALSE)
#knots=attr(xx, "knots")
#max(abs(xx-b.rcs(x, knots)))

#example
#matplot(x,b.rcs(x, knots))
