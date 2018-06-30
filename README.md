peRiodic
===========

Code for the periodic functions, under development.

Install with :

```
library(devtools)
install_github("crtahlin/peRiodiCS")
```

Example models for different variants:

```
# assuming y is the binary response vector and x is numerical vector with
# values to transform to new basis

# RCS - restricted cubic splines
mod.rcs <- glm(y ~ rcs(x), family = "binomial", data = your_data)

# periodic RCS
mod.rcs.per <- glm(y ~ rcs_per(x, xmin = 0, xmax = 1, nk = 5),
                                  family = "binomial", data = your_data)
                                  
# periodic CS (cubic spline)
mod.cs.per <- glm(y ~ cs_per(x, xmin = 0, xmax = 1, nk = 5),
                                  family = "binomial", data = your_data)

```