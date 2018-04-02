# test adding tmax ####

# a simluation without using tmax
set.seed(1234)
result1 <- f.sim.per.splines()
save(result1, file = "result1.Rdata")
str(result1)


# another, but specify prob fnuction at beginning
set.seed(1234)
result2 <- f.sim.per.splines()
save(result2, file = "result2.Rdata")
str(result2)

identical(result1, result2)
# juhu!

# nuw just add that tmax is used in runif, keep it at 1
set.seed(1234)
result3 <- f.sim.per.splines(tmax = 1)
save(result3, file = "result3.Rdata")
str(result3)

identical(result1, result3)

# now try collapsing data, but keep tmax =1 
set.seed(1234)
result4 <- f.sim.per.splines(tmax = 1)
save(result4, file = "result4.Rdata")
str(result4)

identical(result1, result4)
# juhu!

# now try if collapsing works
set.seed(1234)
result5 <- f.sim.per.splines(tmax = 2)
save(result5, file = "result5.Rdata")
str(result5)

identical(result1, result5)
# off course, it should not be identical ... and it isn't

apply(X= result5$my.res, FUN = function(x) {mean(x)}, MARGIN = 2)
apply(X= result1$my.res, FUN = function(x) {mean(x)}, MARGIN = 2)
# these are somewhat different... will try a larger number of simulations (without setting seed?)
result1.large <- f.sim.per.splines(B=5000, tmax = 1)
result5.large <- f.sim.per.splines(B=5000, tmax = 2)
result1.large.averages <- apply(X= result1.large$my.res, FUN = function(x) {format(mean(x), digits = 5) }, MARGIN = 2)
result5.large.averages <- apply(X= result5.large$my.res, FUN = function(x) {format(mean(x), digits = 5) }, MARGIN = 2)
# they seem close enough ...
save(result5.large, file = "result5.large.Rdata")
save(result1.large, file = "result1.large.Rdata")
result6.large <- f.sim.per.splines(B=5000, tmax = 3) # še merganje 3 ciklov
result6.large.averages <- apply(X= result6.large$my.res, FUN = function(x) {format(mean(x), digits = 5) }, MARGIN = 2)
save(result6.large, file = "result6.large.Rdata")
# check difference
as.numeric(result5.large.averages) - as.numeric(result1.large.averages)
as.numeric(result6.large.averages) - as.numeric(result1.large.averages)

# check after adding  support for trends (but trend turned off)
set.seed(1234)
result7 <- f.sim.per.splines()
save(result7, file = "result7.Rdata")
str(result7)

identical(result1$my.res, result7$my.res)
# the main result still the same

# check with trend turned on
set.seed(1234)
result8 <- f.sim.per.splines(add_trend = TRUE)
save(result8, file = "result8.Rdata")
str(result8)

identical(result1$my.res, result8$my.res)
# these are no longer identical, and they shouldn't be

# use default trend, make it 3 years
set.seed(1234)
result9 <- f.sim.per.splines(add_trend = TRUE, tmax = 3)
save(result9, file = "result9.Rdata")
str(result9)
result9.averages <- apply(X= result9$my.res, FUN = function(x) {format(mean(x), digits = 5) }, MARGIN = 2)
result9.averages
result9.large <- f.sim.per.splines(add_trend = TRUE, tmax = 3, B=5000)


