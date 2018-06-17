#' @title RCS.PER, one restriction
#' @description 
#' Experimental: 
#' periodic restricted cubic spline, with only the constraint on equality of estimated values at beginning/end of period
#' 
#' @export
rcs.per.1c=function(x, knots, xmin, xmax){
  #x: numerical variable
  #vector with the knots of the spline
  #xmin: value of the (theoretical) minimum of x
  #xmax: value of the (theoretical) maximum of x
  
  nk=length(knots)
  
  b.x.all=b.rcs(x, knots) #matrix with the expansion of the splines, n*(k-2)
  b.xmax.all=b.rcs(xmax, knots) #value of the spline for x=xmax, vector 1*(k-2)
  
  
  #terms to add in gamma1
  gamma1.c1=(xmin-xmax)/b.xmax.all[nk-2] #constant
  gamma.f1=gamma1.c1*b.x.all[,nk-2] #vector n*1
  
  #betaj.f1=b.x.all[,c(1:(nk-3))]-b.xmax.all[,c(1:(nk-3))]/b.xmax.all[nk-2]*b.x.all[,nk-2]
  
  beta.j=matrix(NA, ncol=nk-3, nrow=length(x))
  for(j in 1:(nk-3)) {beta.j[,j]=b.xmax.all[,j]/b.xmax.all[nk-2]*b.x.all[,nk-2]}
  
  
  
  cbind(x+gamma.f1,
        b.x.all[,1:(nk-3)]-beta.j
  )
}#end f.get.design.rcs.c1
