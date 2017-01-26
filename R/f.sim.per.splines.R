#' @title Simulation function
#' 
#' @param B number of simulations
#' @param n number of points in training sample?
#' @param n.new number of points in test sample?
#' @param nk number of knots
#' @param knots vector with locations of knots
#' @param quatiles.cs quantiles at which to place the knots by default?
#' @param par1sin constant to add to "shift" the sine function to generate values in acceptable range (defaults to 1, to have function values between 0 and 2)
#' @param par2sin factor to multiply the sine function by (<1, "squashing" the function min-max difference; defaults to 0.25, to have mix&max of the output differ by 0.5 from original 2)
#' @param par3sin constant to add to additionaly "shift" the whole function to wanted level (defaults to 0.25, to raise the whole function above 0)
#' @param par4sin the number of sine cycles per one of our periods (defaults to 1 cycle)
#' @param par5sin offset of the sine curve in multiples of 2*Pi (defaults to 0); with 0.25 we would get cosine curve
#' @param tmax the end of the interval we are using, e.g. 3 for "3 years"; useful if studying trends (defaults to 1)
#' @param add_trend boolean. whether to simulate a trend in probability changeing over years. Trend also uses a sine function. defaults to FALSE.
#' @param par1trend constant to add to "shift" the trend part. (defaults to 1, to have function values between 0 and 2)
#' @param par2trend factor to multiply the trend sine function (<1, squashin it; defaults to 0.1, to have min&max differ by 0.2)
#' @param par3trend constant to add to additionaly shift the whole trend to wanted level (defaults to 0, for no additional shifting, e.g. start at 0)
#' @param par4trend the number of trend sine cycles per one of our periods (cycles) (defaults to 1/10, for repeating after 10 cycles)
#' @param par5trend offset of the trend sine curve in multiples of 2*Pi (defaults to 0, for a sine curve)
#' @param max_prob_value the maximum generated allowed probability value (defaults to 0.95)
#' @param min_prob_value the minimum allowed generated probability value (defaults to 0.05)
#' @param set_seed sets the seed - the same for each simulation "batch" (defaults to 1234)
#' 
#' @export
f.sim.per.splines=function(B=100,
                           n=100,
                           n.new=100,  
                           nk=5,
                           knots=NULL, 
                           quantiles.cs = NULL, # c(0.05, 0.25, 0.50, 0.75, 0.95),
                           #parameters for the generation of the periodic data
                           par1sin=1,
                           par2sin=0.25,
                           par3sin=0.25,
                           par4sin=1,
                           par5sin=0,
                           # number of "cycles" to take
                           tmax = 1,
                           # parameters for trend
                           add_trend = FALSE,
                           par1trend = 1,
                           par2trend = 0.1,
                           par3trend = 0,
                           par4trend = 1/10,
                           par5trend = 0,
                           # min/max allowed probability
                           max_prob_value = 0.95,
                           min_prob_value = 0.05,
                           set_seed = 1234){
  
  # set the seed ####
  set.seed(set_seed)
  
  ############### initialization of the outputs ####
  num.res = 3 + 3 + 2 + 4*2 + 2*3 + 3 + 3 # number of result columns to initialize
  #4*2: AUC, 2*3: calibration, 3: LRT for training models, 3: score test for training models
  my.res = matrix(NA, nrow=B, ncol=num.res)
  
  # add names to columns
  colnames(my.res) <- 
    c("brier.rcs.train", "brier.rcs.per.train", "brier.cs.per.train",
      "brier.rcs", "brier.rcs.per", "brier.cs.per",
      "AUCTrainEst.rcs", "AUCTrainEst.rcs.per", "AUCTrainEst.cs.per",
      "AUCNewEst.rcs", "AUCNewEst.rcs.per", "AUCNewEst.cs.per",
      "cal.rcs.1", "cal.rcs.per.1", "cal.cs.per.1",
      "cal.rcs.2", "cal.rcs.per.2", "cal.cs.per.2",
      "lrt.p.value.rcs.train", "lrt.p.value.rcs.per.train", "lrt.p.value.cs.per.train",
      "score.p.value.rcs.train", "score.p.value.rcs.per.train", "score.p.value.cs.per.train",
      "AUCTrainMax", "AUCNewMax", "prop.events.train", "prop.events.test")
  
  # theoretical min and max of the covariate
  min.th.x=0
  max.th.x=1
  
  # points where the coverage is tested ####
  x.test.points = seq(min.th.x + 0.025, max.th.x - 0.025, by=.025)
  
  num.test.points = length(x.test.points)
  
  #matrices where there the coverage probabilities are stored
  
  #predicted probabilities ####
  prob.coverage.cs.per=prob.coverage.rcs=prob.coverage.rcs.per=matrix(NA, ncol=num.test.points*4, nrow=B)
  
  prob.coverage.cs.per.train=prob.coverage.rcs.train=prob.coverage.rcs.per.train=matrix(NA, ncol=num.test.points*4, nrow=B)
  
  #predicted lp ####
  lp.coverage.cs.per=lp.coverage.rcs=lp.coverage.rcs.per=matrix(NA, ncol=num.test.points*4, nrow=B)
  
  lp.coverage.cs.per.train=lp.coverage.rcs.train=lp.coverage.rcs.per.train=matrix(NA, ncol=num.test.points*4, nrow=B)
  
  # structures to hold models
  models_rcs <- list()
  models_rcs.per <- list()
  models_cs.per <- list()
  
  # define probability function in one place
  probability_function <- 
    function(x,
             par1sin. = par1sin,
             par2sin. = par2sin,
             par3sin. = par3sin,
             par4sin. = par4sin,
             par5sin. = par5sin,
             add_trend. = add_trend,
             par1trend. = par1trend,
             par2trend. = par2trend,
             par3trend. = par3trend,
             par4trend. = par4trend,
             par5trend. = par5trend,
             max_prob_value. = max_prob_value,
             min_prob_value. = min_prob_value) {
      sine_function(x,
                    par1sin=par1sin.,
                    par2sin=par2sin.,
                    par3sin=par3sin.,
                    par4sin=par4sin.,
                    par5sin=par5sin.,
                    add_trend = add_trend.,
                    par1trend = par1trend.,
                    par2trend = par2trend.,
                    par3trend = par3trend.,
                    par4trend = par4trend.,
                    par5trend = par5trend.,
                    max_prob_value = max_prob_value.,
                    min_prob_value = min_prob_value.
                  )
    }
  
  ############### beginning of the iterations
  for(b in 1:B) {
    
    # generate x on the interval allown using uniform distribution
    x = runif(n, min = 0, max = tmax) # we generate 100 x values
    
    # simulation of the probabilities
    x.transf.2 = x.transf = probability_function(x) 
    
    # but from now on, we have to transform the X to keep it inside [0,1] using modulus division, e.g. (x %% 1)
    x <- x %% 1
    
    # simulation of the binary events
    y = ifelse(runif(n) < x.transf.2, 1, 0) # event happens if generated random number [0,1] smaller than simulated probability
    # TODO: why not use rbinom() ?
  
    # transformation on the linear predictor scale
    x.transf.2.lp = log(x.transf.2 / (1 - x.transf.2))
    
    # generation of the design matrices ####
    # classic RCS
    x.rcs = rcspline.eval(x, nk = nk, inclx = TRUE)
    #saving the knots used for the splines, uses the default from rms package
    my.knots = attr(x.rcs, "knots")
    
    #periodic RCS
    x.rcs.per = rcs.per(x, nk = nk, xmin = min.th.x, xmax = max.th.x)
    
    # TODO!: check that my.knot are correct for each of the functions! and that each returns the same knot locations
    # TODO!: set the default quantiles.cs to NULL? otherwise they will have knots at different locations than other variants
    # DONE : set to NULL for the simulation function
    
    # periodic cubic splines
    if(is.null(quantiles.cs)){ 
      x.cs.per = cs.per(x, nk=nk, xmin=min.th.x, xmax=max.th.x)
      my.knots.cs = my.knots
    } else {
      my.knots.cs = as.numeric(quantile(ecdf(x), quantiles.cs))
      x.cs.per = cs.per(x, knots=my.knots.cs, nk=NULL, xmin=min.th.x, xmax=max.th.x)
    }
    
    # save each of the knots
    my.knots.rcs = attr(x.rcs, "knots")
    my.knots.rcs.per = attr(x.rcs.per, "knots")
    my.knots.cs.per = attr(x.cs.per, "knots")
    # TODO: maybe later explicitly use each of these in the test set?
    
    # TODO!: save these models to results? to plot later
    # estimate models
    mod.rcs.per = glm(y~x.rcs.per, family="binomial") # model binary responses with different versions of splines
    mod.rcs = glm(y~x.rcs, family="binomial")
    mod.cs.per = glm(y~x.cs.per, family="binomial")
    
    # save models to results placeholder
    models_rcs[[b]] <- mod.rcs
    models_rcs.per[[b]] <- mod.rcs.per
    models_cs.per[[b]] <- mod.cs.per
    
    # Brier's scores
    brier.rcs.per.train = mean((y - mod.rcs.per$fitted.values)^2) # Brier's score on fitted values (square of diff between real and fitted)
    brier.rcs.train = mean((y-mod.rcs$fitted.values)^2)
    brier.cs.per.train = mean((y-mod.cs.per$fitted.values)^2)
    
    # calculate proportion of events happening
    prop.events.train = mean(y)
    
    ################### predicted probabilities, training ####
    
    # estimated linear predictor with se
    pred.rcs.per.train = predict(mod.rcs.per, se.fit = TRUE) # on the scale of linear predictors
    pred.rcs.train = predict(mod.rcs, se.fit=TRUE)
    pred.cs.per.train = predict(mod.cs.per, se.fit=TRUE)
    
    # saving the linear predictors
    lp.rcs.per.train = pred.rcs.per.train$fit
    lp.rcs.train = pred.rcs.train$fit
    lp.cs.per.train = pred.cs.per.train$fit
    
    # deriving the estimated probabilities
    p.rcs.per.train = 1/(1+exp(-lp.rcs.per.train))
    p.rcs.train = 1/(1+exp(-lp.rcs.train))
    p.cs.per.train = 1/(1+exp(-lp.cs.per.train))
    
    #estimated probabilities predictor with se
    p.pred.rcs.per.train = predict(mod.rcs.per, se.fit = TRUE, type="response")
    p.pred.rcs.train = predict(mod.rcs, se.fit=TRUE, type="response")
    p.pred.cs.per.train = predict(mod.cs.per, se.fit=TRUE, type="response")
    
    ######## coverage #############
    ########### coverage predicted probabilities #######
    
    my.index=numeric(num.test.points)
    
    for(i in 1:num.test.points){
      my.index[i]=which.min(abs(x-x.test.points[i]))}
    # calculate if the real (simulated) point was covered by the confidence interval (and save some other info) for each of the variants 
    prob.coverage.rcs.per.train[b,]=c(p.pred.rcs.per.train$fit[my.index]-1.96*p.pred.rcs.per.train$se.fit[my.index]<=x.transf.2[my.index] & p.pred.rcs.per.train$fit[my.index]+1.96*p.pred.rcs.per.train$se.fit[my.index]>=x.transf.2[my.index],
                                      
                                      p.pred.rcs.per.train$fit[my.index]-1.96*p.pred.rcs.per.train$se.fit[my.index], 
                                      p.pred.rcs.per.train$fit[my.index],
                                      p.pred.rcs.per.train$fit[my.index]+1.96*p.pred.rcs.per.train$se.fit[my.index])
    
    prob.coverage.rcs.train[b,]=c(p.pred.rcs.train$fit[my.index]-1.96*p.pred.rcs.train$se.fit[my.index]<=x.transf.2[my.index] & p.pred.rcs.train$fit[my.index]+1.96*p.pred.rcs.train$se.fit[my.index]>=x.transf.2[my.index],
                                  
                                  p.pred.rcs.train$fit[my.index]-1.96*p.pred.rcs.train$se.fit[my.index], 
                                  p.pred.rcs.train$fit[my.index],
                                  p.pred.rcs.train$fit[my.index]+1.96*p.pred.rcs.train$se.fit[my.index])
    
    
    prob.coverage.cs.per.train[b,]=c(p.pred.cs.per.train$fit[my.index]-1.96*p.pred.cs.per.train$se.fit[my.index]<=x.transf.2[my.index] & p.pred.cs.per.train$fit[my.index]+1.96*p.pred.cs.per.train$se.fit[my.index]>=x.transf.2[my.index],
                                     
                                     p.pred.cs.per.train$fit[my.index]-1.96*p.pred.cs.per.train$se.fit[my.index], 
                                     p.pred.cs.per.train$fit[my.index],
                                     p.pred.cs.per.train$fit[my.index]+1.96*p.pred.cs.per.train$se.fit[my.index])
    
    ################ end predicted probabilities training
    ##### predicted lp coverage, training 
    
    lp.coverage.rcs.per.train[b,]=c(pred.rcs.per.train$fit[my.index]-1.96*pred.rcs.per.train$se.fit[my.index]<=x.transf.2.lp[my.index] & pred.rcs.per.train$fit[my.index]+1.96*pred.rcs.per.train$se.fit[my.index]>=x.transf.2.lp[my.index],
                                    
                                    pred.rcs.per.train$fit[my.index]-1.96*pred.rcs.per.train$se.fit[my.index], 
                                    pred.rcs.per.train$fit[my.index],
                                    pred.rcs.per.train$fit[my.index]+1.96*pred.rcs.per.train$se.fit[my.index])
    
    
    lp.coverage.rcs.train[b,]=c(pred.rcs.train$fit[my.index]-1.96*pred.rcs.train$se.fit[my.index]<=x.transf.2.lp[my.index] & pred.rcs.train$fit[my.index]+1.96*pred.rcs.train$se.fit[my.index]>=x.transf.2.lp[my.index],
                                
                                pred.rcs.train$fit[my.index]-1.96*pred.rcs.train$se.fit[my.index], 
                                pred.rcs.train$fit[my.index],
                                pred.rcs.train$fit[my.index]+1.96*pred.rcs.train$se.fit[my.index])
    
    
    lp.coverage.cs.per.train[b,]=c(pred.cs.per.train$fit[my.index]-1.96*pred.cs.per.train$se.fit[my.index]<=x.transf.2.lp[my.index] & pred.cs.per.train$fit[my.index]+1.96*pred.cs.per.train$se.fit[my.index]>=x.transf.2.lp[my.index],
                                   
                                   pred.cs.per.train$fit[my.index]-1.96*pred.cs.per.train$se.fit[my.index], 
                                   pred.cs.per.train$fit[my.index],
                                   pred.cs.per.train$fit[my.index]+1.96*pred.cs.per.train$se.fit[my.index])
    
    # end predicted lp coverage, training
    
    # AUC max on test data
    AUCTrainMax = as.numeric(performance(prediction(x.transf.2, y), "auc")@y.values)
    
    #AUC estimated on test data
    AUCTrainEst.cs.per = as.numeric(performance(prediction(p.pred.cs.per.train$fit, y), "auc")@y.values)
    
    AUCTrainEst.rcs.per=as.numeric(performance(prediction(p.pred.rcs.per.train$fit, y), "auc")@y.values)
    
    AUCTrainEst.rcs=as.numeric(performance(prediction(p.pred.rcs.train$fit, y), "auc")@y.values)
    
    
    ########## lrt for variable selection
    lrt.p.value.cs.per.train = 1-pchisq(-1*(mod.cs.per$deviance-mod.cs.per$null.deviance), mod.cs.per$df.null-mod.cs.per$df.residual)
    
    lrt.p.value.rcs.per.train = 1-pchisq(-1*(mod.rcs.per$deviance-mod.rcs.per$null.deviance), mod.rcs.per$df.null-mod.rcs.per$df.residual)
    
    lrt.p.value.rcs.train = 1-pchisq(-1*(mod.rcs$deviance-mod.rcs$null.deviance), mod.rcs$df.null-mod.rcs$df.residual)
    
    ########## score for variable selection
    score.p.value.cs.per.train=anova(mod.cs.per, test="Rao")[2,6]
    
    score.p.value.rcs.per.train=anova(mod.rcs.per, test="Rao")[2,6]
    
    score.p.value.rcs.train=anova(mod.rcs, test="Rao")[2,6]
    
    #4debug
    x.old=x
    
    ############## new data (test set) #######################
    
    x=runif(n.new, min = 0, max = tmax)
    
    x.transf.2=x.transf = probability_function(x) 
    x <- x %% 1 # merge into one cycle
    
    ynew=ifelse(runif(n)<x.transf, 1, 0)
    x.transf.2.lp=log(x.transf.2/(1-x.transf.2))
    
    #using the knots specified in the training data
    x.rcs.per = rcs.per(x, knots=my.knots, xmin=min.th.x, xmax=max.th.x) 
    x.rcs = rcspline.eval(x, knots=my.knots, inclx=TRUE)
    x.cs.per = cs.per(x, knots=my.knots.cs, nk=NULL, xmax=max.th.x, xmin=min.th.x)
    
    # construct data frames to hold data to use in test models   
    new.data.rcs.per = data.frame(ynew=ynew, x.rcs.per)
    dimnames(new.data.rcs.per)[[2]][-1]=seq(1:dim(x.rcs.per)[2])
    
    new.data.rcs = data.frame(ynew=ynew, x.rcs)
    dimnames(new.data.rcs)[[2]][-1]=seq(1:dim(x.rcs)[2])
    
    new.data.cs.per=data.frame(ynew=ynew, x.cs.per)
    dimnames(new.data.cs.per)[[2]][-1]=seq(1:dim(x.cs.per)[2])
    
    #estimated linear predictor with se
    pred.rcs.per = predict(mod.rcs.per, newdata=new.data.rcs.per, se.fit = TRUE)
    pred.rcs = predict(mod.rcs, newdata=new.data.rcs, se.fit=TRUE)
    pred.cs.per = predict(mod.cs.per, newdata=new.data.cs.per, se.fit=TRUE)
    
    #saving the linear predictors
    lp.rcs.per = pred.rcs.per$fit
    lp.rcs = pred.rcs$fit
    lp.cs.per = pred.cs.per$fit
    
    #deriving the estimated probabilities
    p.rcs.per = 1/(1+exp(-lp.rcs.per))
    p.rcs = 1/(1+exp(-lp.rcs))
    p.cs.per = 1/(1+exp(-lp.cs.per))
    
    #estimated probabilities predictor with se
    p.pred.rcs.per = predict(mod.rcs.per, newdata=new.data.rcs.per, se.fit = TRUE, type="response")
    p.pred.rcs = predict(mod.rcs, newdata=new.data.rcs, se.fit=TRUE, type="response")
    p.pred.cs.per = predict(mod.cs.per, newdata=new.data.cs.per, se.fit=TRUE, type="response")
      
    ############ brier's score on new data ###############
    brier.rcs.per = mean((ynew-p.rcs.per)^2)
    brier.rcs = mean((ynew-p.rcs)^2)
    brier.cs.per = mean((ynew-p.cs.per)^2)
    
    #proportion of events on new data
    prop.events.test=mean(ynew)
    
    ############ Prediction CI coverage score on new data ###############
    
    my.index=numeric(num.test.points)
    
    for(i in 1:num.test.points){
      my.index[i]=which.min(abs(x-x.test.points[i]))}
    
    prob.coverage.rcs.per[b,]=c(p.pred.rcs.per$fit[my.index]-1.96*p.pred.rcs.per$se.fit[my.index]<=x.transf.2[my.index] & p.pred.rcs.per$fit[my.index]+1.96*p.pred.rcs.per$se.fit[my.index]>=x.transf.2[my.index],
                                
                                p.pred.rcs.per$fit[my.index]-1.96*p.pred.rcs.per$se.fit[my.index], 
                                p.pred.rcs.per$fit[my.index],
                                p.pred.rcs.per$fit[my.index]+1.96*p.pred.rcs.per$se.fit[my.index])
    
    
    prob.coverage.rcs[b,]=c(p.pred.rcs$fit[my.index]-1.96*p.pred.rcs$se.fit[my.index]<=x.transf.2[my.index] & p.pred.rcs$fit[my.index]+1.96*p.pred.rcs$se.fit[my.index]>=x.transf.2[my.index],
                            
                            p.pred.rcs$fit[my.index]-1.96*p.pred.rcs$se.fit[my.index], 
                            p.pred.rcs$fit[my.index],
                            p.pred.rcs$fit[my.index]+1.96*p.pred.rcs$se.fit[my.index])
    
    
    prob.coverage.cs.per[b,]=c(p.pred.cs.per$fit[my.index]-1.96*p.pred.cs.per$se.fit[my.index]<=x.transf.2[my.index] & p.pred.cs.per$fit[my.index]+1.96*p.pred.cs.per$se.fit[my.index]>=x.transf.2[my.index],
                               
                               p.pred.cs.per$fit[my.index]-1.96*p.pred.cs.per$se.fit[my.index], 
                               p.pred.cs.per$fit[my.index],
                               p.pred.cs.per$fit[my.index]+1.96*p.pred.cs.per$se.fit[my.index])
    
    
    ############# end coverage of predicted probabilities
    
    
    #### coverage of linear predictor, new data
    
    ##### predicted lp coverage, training 
   lp.coverage.rcs.per[b,]=c(pred.rcs.per$fit[my.index]-1.96*pred.rcs.per$se.fit[my.index]<=x.transf.2.lp[my.index] & pred.rcs.per$fit[my.index]+1.96*pred.rcs.per$se.fit[my.index]>=x.transf.2.lp[my.index],
                              
                              pred.rcs.per$fit[my.index]-1.96*pred.rcs.per$se.fit[my.index], 
                              pred.rcs.per$fit[my.index],
                              pred.rcs.per$fit[my.index]+1.96*pred.rcs.per$se.fit[my.index])
    
    
    lp.coverage.rcs[b,]=c(pred.rcs$fit[my.index]-1.96*pred.rcs$se.fit[my.index]<=x.transf.2.lp[my.index] & pred.rcs$fit[my.index]+1.96*pred.rcs$se.fit[my.index]>=x.transf.2.lp[my.index],
                          
                          pred.rcs$fit[my.index]-1.96*pred.rcs$se.fit[my.index], 
                          pred.rcs$fit[my.index],
                          pred.rcs$fit[my.index]+1.96*pred.rcs$se.fit[my.index])
    
    
    lp.coverage.cs.per[b,]=c(pred.cs.per$fit[my.index]-1.96*pred.cs.per$se.fit[my.index]<=x.transf.2.lp[my.index] & pred.cs.per$fit[my.index]+1.96*pred.cs.per$se.fit[my.index]>=x.transf.2.lp[my.index],
                             
                             pred.cs.per$fit[my.index]-1.96*pred.cs.per$se.fit[my.index], 
                             pred.cs.per$fit[my.index],
                             pred.cs.per$fit[my.index]+1.96*pred.cs.per$se.fit[my.index])
    
    ############## end of coverage of linear predictor, new data
    
   # TODO: are these three lines duplicated from above (not needed)? 
    brier.rcs.per=mean((ynew-p.rcs.per)^2)
    brier.rcs=mean((ynew-p.rcs)^2)
    brier.cs.per=mean((ynew-p.cs.per)^2)
    
    #AUC max on test data
    AUCNewMax = as.numeric(performance(prediction(x.transf.2, ynew), "auc")@y.values)
    
    #AUC estimated on test data
    AUCNewEst.cs.per = as.numeric(performance(prediction(p.pred.cs.per$fit, ynew), "auc")@y.values)
    
    AUCNewEst.rcs.per = as.numeric(performance(prediction(p.pred.rcs.per$fit, ynew), "auc")@y.values)
    
    AUCNewEst.rcs = as.numeric(performance(prediction(p.pred.rcs$fit, ynew), "auc")@y.values)
    
    #calibration on new data
    
    #intercept and slope, 6 new parameters
    cal.cs.per=glm(ynew~pred.cs.per$fit, family="binomial")$coef
    cal.rcs.per=glm(ynew~pred.rcs.per$fit, family="binomial")$coef
    cal.rcs=glm(ynew~pred.rcs$fit, family="binomial")$coef
    
    my.res[b,]=c(brier.rcs.train, brier.rcs.per.train, brier.cs.per.train, 
                 brier.rcs, brier.rcs.per, brier.cs.per, 
                 AUCTrainEst.rcs, AUCTrainEst.rcs.per, AUCTrainEst.cs.per, 
                 AUCNewEst.rcs, AUCNewEst.rcs.per, AUCNewEst.cs.per,  
                 cal.rcs[1], cal.rcs.per[1], cal.cs.per[1],
                 cal.rcs[2], cal.rcs.per[2], cal.cs.per[2],
                 lrt.p.value.rcs.train, lrt.p.value.rcs.per.train, lrt.p.value.cs.per.train,
                 score.p.value.rcs.train, score.p.value.rcs.per.train, score.p.value.cs.per.train,
                 AUCTrainMax, AUCNewMax, prop.events.train, prop.events.test)
    
    
    
  }#end for b
  
  # TODOs
  # -pojdi skozi izvajanje, preglej rezultate
  # -dodaj možnost generiranja podatkov z zamaknjeno sin funkcijo
  # -dodaj možnost izvajanja z non-sin funkcijo?
  # -dodaj možnost izvajanja 
  
  return(list(my.res=my.res, # results
              # coverages 
              prob.coverage.rcs=prob.coverage.rcs,
              prob.coverage.rcs.per=prob.coverage.rcs.per,
              prob.coverage.cs.per=prob.coverage.cs.per,
              lp.coverage.rcs=lp.coverage.rcs,
              lp.coverage.rcs.per=lp.coverage.rcs.per,
              lp.coverage.cs.per=lp.coverage.cs.per,
              prob.coverage.rcs.train=prob.coverage.rcs.train,
              prob.coverage.rcs.per.train=prob.coverage.rcs.per.train,
              prob.coverage.cs.per.train=prob.coverage.cs.per.train,
              lp.coverage.rcs.train=lp.coverage.rcs.train,
              lp.coverage.rcs.per.train=lp.coverage.rcs.per.train,
              lp.coverage.cs.per.train=lp.coverage.cs.per.train,
              # points at which coverage was tested
              x.test.points=x.test.points, num.test.points=num.test.points, 
              #saving the parameters of the simulation
              B=B, n=n, n.new=n.new, nk=nk, knots=my.knots, knots.cs=my.knots.cs, 
              quantiles.cs=quantiles.cs,
              #parameters for the generation of the periodic data
              par1sin=par1sin, par2sin=par2sin, par3sin=par3sin, par4sin=par4sin, par5sin=par5sin, prop.events.train=prop.events.train, prop.events.test=prop.events.test,
              x=x, x.transf.2.lp=x.transf.2.lp,
              tmax = tmax,
              # trend parameters
              add_trend = add_trend,
              par1trend = par1trend,
              par2trend = par2trend,
              par3trend = par3trend,
              par4trend = par4trend,
              par5trend = par5trend,
              max_prob_value = max_prob_value,
              min_prob_value = min_prob_value,
              # seed
              set_seed = set_seed,
              # knots info explicitly
              my.knots.rcs = my.knots.rcs,
              my.knots.rcs.per = my.knots.rcs.per,
              my.knots.cs.per = my.knots.cs.per #,
              # estimated models for each simulation 
              # NOTE: remove temporarily, due to large size of results
              # models_rcs = models_rcs,
              # models_rcs.per = models_rcs.per,
              # models_cs.per = models_cs.per
  ))
  
  
  
}#end function sim