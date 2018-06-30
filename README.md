[![Travis-CI Build Status](https://travis-ci.org/crtahlin/peRiodiCS.svg?branch=master)](https://travis-ci.org/crtahlin/peRiodiCS)


peRiodiCS
===========

Code for periodic versions of cubic splines and restricted cubic splines.
The functions transform a numeric (time) variable `x` into a new basis, 
to model with either a periodic variant of cubic splines or
periodic restricted cubic spline. 

Install from GitHub with:

```
library(devtools)
install_github("crtahlin/peRiodiCS")
```

Example models for different variants:

```
# assuming y is a binary response vector and x is numerical vector with
# values to transform to new basis

# RCS - restricted cubic splines (non-periodic, just for reference)
mod.rcs <- glm(y ~ rcs(x), family = "binomial", data = your_data)

# periodic RCS
mod.rcs.per <- glm(y ~ rcs_per(x, xmin = 0, xmax = 1, nk = 5),
                                  family = "binomial", data = your_data)
                                  
# periodic CS (cubic spline)
mod.cs.per <- glm(y ~ cs_per(x, xmin = 0, xmax = 1, nk = 5),
                                  family = "binomial", data = your_data)

```