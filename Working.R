# Import data ####
# AIRS data frame
load("./extdata//AIRSData.Rdata")
head(AIRSData)

# load the library
library(cyclicCurvesAlfa)

# cyclic cubuc splines
cyclicCubicSplineTPB(AIRSData[sample(1:1699655, size=1500),], coef = TRUE)
cyclicCubicSplineTPB(AIRSData[sample(1:1699655, size=1500),], coef = FALSE)
cyclicCubicSplineTPBCurve(AIRSData[sample(1:1699655, size=2000),])
cyclicCubicSplineTPBCurve(AIRSData[sample(1:1699655, size=1000),], curve = TRUE)

# cyclic fractional polynomials
cyclicFractionalPolynomial(AIRSData)
cyclicFractionalPolynomialCurve(AIRSData[sample(1:1699655, size=10000),], curve=TRUE, add=FALSE)
