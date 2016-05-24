# Import data ####
# AIRS data frame
load("./extdata//AIRSData.Rdata")
head(AIRSData)

# load the library
library(cyclicCurvesAlfa)

# cyclic cubic splines
CCS1 <- cyclicCubicSplineTPB(AIRSData[sample(1:1699655, size=1500),], coef = TRUE)
str(CCS1) # "solution" matrix
CCS2 <- cyclicCubicSplineTPB(AIRSData[sample(1:1699655, size=1500),], coef = FALSE)
str(CCS2) # regression coeficients 
# curve
cyclicCubicSplineTPBCurve(AIRSData[sample(1:1699655, size=2000),])
# curve with confidence intervals
cyclicCubicSplineTPBCurve(AIRSData[sample(1:1699655, size=1000),], se = TRUE)

# cyclic fractional polynomials
cyclicFractionalPolynomial(AIRSData) # regression coeficients
# curve
cyclicFractionalPolynomialCurve(AIRSData[sample(1:1699655, size=10000),], curve=TRUE, add=FALSE)
