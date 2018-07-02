### CS PER --------------------------------------

# load example data; see help("viral_east_mediteranean")
data("viral_east_mediteranean")

# calculate location of knots to use
Knots <- 
  Hmisc::rcspline.eval(x = viral_east_mediteranean$EpiWeek,
                       nk = 5, knots.only = TRUE)

# model viral infections vs weeks
model <- glm(RSV ~ cs_per(EpiWeek, knots = Knots), data = viral_east_mediteranean)

# plot model (with many points, to make it smooth)
plot_per_mod(Model = model, XvarName = "EpiWeek", Smooth = TRUE)


### RCS -------------------------------------------

# load example data; see help("viral_east_mediteranean")
data("viral_east_mediteranean")

# calculate location of knots to use
Knots <- 
  Hmisc::rcspline.eval(x = viral_east_mediteranean$EpiWeek,
                       nk = 5, knots.only = TRUE)

# model viral infections vs weeks
model <- glm(RSV ~ Hmisc::rcspline.eval(EpiWeek, inclx = TRUE, knots = Knots), data = viral_east_mediteranean)

# plot model (with many points, to make it smooth)
plot_per_mod(Model = model, XvarName = "EpiWeek", Smooth = TRUE)



################ TEST EXAMPLE IN README

# load example data; see help("viral_east_mediteranean")
data("viral_east_mediteranean")

# calculate location of knots to use
Knots <- Hmisc::rcspline.eval(x = viral_east_mediteranean$EpiWeek,
                              nk = 5, knots.only = TRUE)

### RCS - restricted cubic splines (non-periodic, just for reference) #########
# model viral infections vs weeks
model_rcs <- glm(RSV ~ Hmisc::rcspline.eval(EpiWeek, inclx = TRUE, knots = Knots), data = viral_east_mediteranean)

# plot model (with many points, to make it smooth)
plot_per_mod(Model = model_rcs, XvarName = "EpiWeek", Smooth = TRUE)

### periodic RCS ##############################################################
# model viral infections vs weeks
model_rcs_per <- glm(RSV ~ rcs_per(EpiWeek, knots = Knots), data = viral_east_mediteranean)

# plot model (with many points, to make it smooth)
plot_per_mod(Model = model_rcs_per, XvarName = "EpiWeek", Smooth = TRUE)

### periodic CS (cubic spline) ###############################################
# model viral infections vs weeks
model_cs_per <- glm(RSV ~ cs_per(EpiWeek, knots = Knots), data = viral_east_mediteranean)

# plot model (with many points, to make it smooth)
plot_per_mod(Model = model_cs_per, XvarName = "EpiWeek", Smooth = TRUE)
