# Import data ####
# AIRS data frame
load("./extdata//AIRSData.Rdata")
head(AIRSData)

# load the library
library(cyclicCurvesAlfa)

# cyclic cubuc splines
cyclicCubicSplineTPB(AIRSData[sample(1:1699655, size=15000),], coef = TRUE)
cyclicCubicSplineTPBCurve(AIRSData[sample(1:1699655, size=200000),])
cyclicCubicSplineTPBCurve(AIRSData[sample(1:1699655, size=10000),], curve = TRUE)

# cyclic fractional polynomials
cyclicFractionalPolynomial(AIRSData)
cyclicFractionalPolynomialCurve(AIRSData[sample(1:1699655, size=10000),], curve=TRUE, add=FALSE)
