# example code to export a design design matrix to use with periodic splines
# TODO: install package from CRAN ; currently GitHub

library(devtools) #load library for installing from GitHub and other development functionalities
install_github("crtahlin/peRiodic") # install package for periodic splines


library(peRiodiCS) # load package for periodic splines

# you can use your data, but in this example we will use data available in the package
??peRiodiCS::viral_east_mediteranean # for more info about the dataset
data("viral_east_mediteranean") # load the data set
viral_east_mediteranean$EpiWeek # contains the week of measurement (1-53)
summary(viral_east_mediteranean$EpiWeek) # summary of the data

# there are some data missing, keep only rows with existing data
viral_nonNAData <- na.omit(viral_east_mediteranean[, c("RSV","EpiWeek")])
viral_nonNAData$EpiWeek <- as.numeric(viral_nonNAData$EpiWeek)
str(viral_nonNAData)


#############################
# RCS - restricted cubic splines
design.matrix.rcs <- rcspline.eval(viral_nonNAData$EpiWeek, nk = 5, inclx = TRUE) # 5 knots, includes X values (has to be explicitly specified)
print(design.matrix.rcs) # prints the matrix on screen
write.csv2(design.matrix.rcs, file = "design.matrix.rcs.csv") # write to CSV file

# build a model, just for example
Knots <- rcspline.eval(x = viral_nonNAData$EpiWeek, nk = 5, knots.only = TRUE)
mod.rcs <- glm(RSV ~ rcs(EpiWeek, knots = Knots, inclx = TRUE), family = "binomial", data = viral_nonNAData, x = TRUE)
Plot.per.mod(mod.rcs, XvarName = "EpiWeek", Smooth = FALSE, xLocation = 2)
Plot.per.mod(mod.rcs, XvarName = "EpiWeek", Smooth = TRUE) # doesn't work
mod.rcs.2 <- glm(RSV ~ rcspline.eval(EpiWeek, knots = Knots, inclx = TRUE), family = "binomial", data = viral_nonNAData, x = TRUE)
Plot.per.mod(mod.rcs.2, XvarName = "EpiWeek", Smooth = TRUE)

#############################
# RCS.PER - periodic restricted cubic spline
design.matrix.rcs.per <- rcs.per(viral_nonNAData$EpiWeek, nk = 5) # note X is not included in the design matrix
print(design.matrix.rcs.per)
write.csv2(design.matrix.rcs.per, file = "design.matrix.rcs.per.csv") # write to CSV file

# build a model, just for example
Knots <- rcspline.eval(x = viral_nonNAData$EpiWeek, nk = 5, knots.only = TRUE)
mod.rcs.per <- glm(RSV ~ rcs.per(EpiWeek, knots = Knots), family = "binomial", data = viral_nonNAData) 
Plot.per.mod(mod.rcs.per, XvarName = "EpiWeek", Smooth = TRUE)

#############################
# CS.PER - periodic cubic splines
design.matrix.cs.per <- cs.per(viral_nonNAData$EpiWeek, nk = 5, xmin = 1, xmax = 53)
print(design.matrix.cs.per)
write.csv2(design.matrix.cs.per, file = "design.matrix.cs.per.csv") # write to CSV file

# build a model, just for example
Knots <- rcspline.eval(x = viral_nonNAData$EpiWeek, nk = 5, knots.only = TRUE)
mod.cs.per <- glm(RSV ~ cs.per(EpiWeek, knots = Knots), data = viral_nonNAData, family = "binomial", x = TRUE, y = TRUE )
Plot.per.mod(mod.cs.per, XvarName = "EpiWeek", Smooth = TRUE)

