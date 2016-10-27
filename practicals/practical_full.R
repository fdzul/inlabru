#'---
#' title: "INLAbru practical"
#' author: Fabian E. Bachl, Yuan Yuan, Lindesay Scott-Hayward, Janine Illian, David Borchers, Finn Lindgren
#' output:
#'   html_document:
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'       smooth_scroll: false
#'---

#'## Intrduction
#'
#' This practical is also available as [R code](practical_full.R) and [R markdown document](practical_full.Rnw).
#'


#'## Initial settings

#' Load the package
#' 
#+results="hide",message=FALSE

  devtools::load_all("/home/fbachl/devel/git/essmod/iDistance")

#' Special settings for the practical
#' 
#'* Use emprical Bayes for inference (faster but not the real deal)
#'* Compute DIC and WAIC criteria
#'* Settings to make predictions possible
#'
#+results="hide",message=FALSE
  
  init.tutorial()
#'


#' Random fields in one dimension
#'-------------------------------------------------------------------------------------

#'### The INLA approach

#' Build a mesh

  x = seq(0, 50, by = 1)
  mesh <- inla.mesh.1d(x, degree = 1)

#' Invent an intensity in one dimension

  lambda = function(x) { 5*exp(sin(0.2*x))}

#' Sample from intensity

  pts = sample.lgcp(mesh, log(lambda(x)))

#' Plot samples

  ggplot(data=data.frame(x=x,y=lambda(x))) + geom_line(aes(x=x,y=y)) + 
    geom_point(data = pts, aes(x=x), y = 2, shape = "|")

#' Make a 1D SPDE model
#' 
#'* Note: `g()` is the same as inla `f()` except that the mesh has to be stated explicitly.
#'* Using `f()` does not work if meshes are involved.
#'* SPDE parameters go into the inla.spde2.matern() bit. 

  mdl = x ~ g(x, model = inla.spde2.matern(mesh), mesh = mesh)

#' Run LGCP inference
#' 
#+results="hide"

  r = lgcp(points = pts, model = mdl)

#' A very concise summary

  summary(r)

#' Predict intensity
#' The n parameter sets the number of samples used for the prediction. Using more gives
#' better predictions but is also slower.

  intensity = predict(r, x ~ exp(x + Intercept), n = 5000)

#' Plot the ground truth intensity as well as our estimate

  plot(intensity) + geom_line(data=data.frame(x=x,y=lambda(x)), aes(x=x,y=y), linetype = 2, alpha = 0.5)

#'### Comparison to GAMs

#' Bin x-coordinates of data and compute count per bin

  breaks = seq(0,50, length.out = 101)
  binwidth = breaks[2] - breaks[1]
  hst = hist(pts$x, breaks = breaks, plot = FALSE)
  gamdata = data.frame(x = hst$mid, count = hst$count, exposure = binwidth)

#' GAM inference
#' 
#+results="hide" 
  
  library(mgcv)
  r.gam = gam(count ~ s(x), data = gamdata, family = poisson())

#' Predict

  newdata = data.frame(x = seq(0,50,length.out = 100))
  intensity.gam = exp(predict(r.gam, newdata = newdata))

#' Compare to INLA
#'
#'* The step function (dotted line) represents the histogram of the binned data
#'* The blue line shows the GAM fit

  plot(intensity) + 
    geom_line(data=data.frame(x=x,y=lambda(x)), aes(x=x,y=y), linetype = 2, alpha = 0.5) +
    geom_line(data = data.frame(x=newdata$x, y = intensity.gam), aes(x,y), color = "blue") +
    geom_step(data=gamdata, aes(x=x-0.5*exposure,y=count/exposure), linetype = 3, alpha = 0.5)

#'* While the confidence band around the INLA prediction reflects the variability of the histogram,
#'  a quantitative comparison is discuraged. More on this during the goodness-of-fit practical!
#'* The GAM result depends strongly on the number of breaks. This is quite interesting!

#'### TO DO: 
#'
#'* Adjust SPDE parameters to correspond to GAM smoother

#'


#' Random fields in two dimensions
#'-------------------------------------------------------------------------------------

#' Let's create a know intensity function. The samples below cover a domain of 
#' x in [-5,5] and y [-5,5], so chose an intensity that makes sense on that domain.

  lambda = function(loc) { 1000 * dnorm(coordinates(loc)[,1], sd = 2) * dnorm(coordinates(loc)[,2], sd = 2)}
  
#' This is a little helper function that will generate samples
  
  sc = toy.completesample(lambda)

#' Plot points and detections

  ggplot() + gg.mesh(sc$mesh) + gg.point(sc$points)

#' Run `lgcp`
#' 
#+results="hide"

  r = lgcp(sc$points, mesh = sc$mesh)

#' Summary

  summary(r)

#' Predict field

  intensity = predict(r, coordinates ~ exp(Intercept + spde))
  plot(intensity)

#' Predict abundance

  abundance = predict(r, ~ exp(Intercept + spde), integrate = "coordinates")
  abundance
  plot(abundance)

#' Set up model manually if you want to adjust the hyper parameters
#' 
#'* The first argument of `g()` only serves as a name of the SPDE model (see predict call later on)
#'* Intercept is added automatically unless the right hand side contains -1
#'* To do: use SPDE parameterization that reproduces GAM smoother (via pcprior)

  spde.mdl = inla.spde2.matern(sc$mesh, prior.variance.nominal = 10, theta.prior.prec = 0.01)
  mdl = coordinates ~ g(spat, model = spde.mdl, mesh = sc$mesh)

#' Run `lgcp`
#' 
#+results="hide"
  
  r = lgcp(sc$points, model = mdl)

#' Predict intensity field

  intensity = predict(r, coordinates ~ exp(Intercept + spat))
  plot(intensity)



#' TO DO: 
#'
#'* Adjust SPDE parameters to correspond to GAM smoother
#'* Use GAM smoother for comparison
#'


#' Sampling strategies
#'-------------------------------------------------------------------------------------
require(maptools) # Needed for toy data generator
lambda = function(loc) { 1000 * dnorm(coordinates(loc)[,1], sd = 2) * dnorm(coordinates(loc)[,2], sd = 2)}

#'### Point sampling
#'
#' The toy data generator for plot sampling

  sp = toy.pointsample(lambda, area = 0.2^2, n = 300)

#' Plot samplers and detections
  ggplot() + gg.mesh(sp$mesh) + 
    gg.point(sp$samplers, size = 5, color = "black", alpha = 0.1) +
    gg.point(sp$points, size = 2, color = "red")
  

#' LGCP inference
#' Since the domain was observed through plot sampling we need to add the `samplers` 
#' argument to our `lgcp` call:

  r = lgcp(points = sp$points, samplers = sp$samplers, mesh = sp$mesh)
  
#' Predict spatial intensity

  spintens = predict(r, coordinates ~ exp(spde + Intercept))
  plot(spintens) +
    gg.point(sp$samplers, size = 5, color = "black", alpha = 0.2) +
    gg.point(sp$points, size = 2, color = "red")
  
#' Predict total intensity
  
  sptintens = predict(r, ~ exp(spde + Intercept), integrate = "coordinates")
  sptintens
  

#'### Line sampling
#'
#' Sample from SpatialLines

  sl = toy.linesample(lambda, n = 5, width = 0.1)

#' Plot lines and detections

  ggplot() + gg.mesh(sl$mesh) + 
    gg.segment(sl$samplers, size = 3, alpha = 0.2) +
    gg.point(sl$points, size = 1.5, color = "red")
  
#' LGCP inference

  r = lgcp(points = sl$points, samplers = sl$lines, mesh = sl$mesh)

#' Predict spatial intensity
  
  slintens = predict(r, coordinates ~ exp(spde + Intercept))
  plot(slintens) +
    gg.segment(sl$samplers, size = 3, alpha = 0.2) +
    gg.point(sl$points, size = 1.5, color = "red")
  
#' Predict total intensity

  sltintens = predict(r, ~ exp(spde + Intercept), integrate = "coordinates")  
  sltintens

#'### Areal sampling
#' 
#' Generate toy data

  sa = toy.polysample(lambda, n = 30, seed = 3)

#' Plot polgygons and detections

  ggplot() + 
    gg.mesh(sa$mesh) + 
    gg.polygon(sa$samplers) + 
    gg.point(sa$points, size = 1, color = sa$points$sampler)

#' Run LGCP

  r = lgcp(points = sa$points, samplers = sa$samplers, mesh = sa$mesh)

#' Predict spatial intensity

  saintens = predict(r, coordinates ~ exp(spde + Intercept))
  
#' Plot spatial intensity

  plot(saintens) +
    gg.polygon(sa$samplers) +
    gg.point(sa$points, size = 1, color = "red")
  
#' Predict total intensity
  
  satintens = predict(r, ~ exp(spde + Intercept), integrate = "coordinates")  
  satintens

  
#' Compare total intensities for sampling methods

  plot(sptintens,  satintens, sltintens)
  
#' Predict counts (=abundance) and compare with total intensity posterior
#' 
#'* Difference mostly due to low sample size
#'* But that depends on the scale of the problem

  sacnt = predict(r, ~ rpois(rep(1,length(spde)), exp(spde + Intercept)), integrate = "coordinates")  
  plot(sacnt, satintens)
  
#'### Count data
#'
#' The toy data generator for point sampling stores the number of points 
#' found within each sampling area. We can use this to run a regression
#' on counts.

  sp = toy.pointsample(lambda, area = 0.2^2, n = 300)
  
#' Plot samplers and the number of points found by each sampler

  ggplot() + gg.mesh(sp$mesh) +
    gg.point(sp$samplers, size = 5, color = "black", alpha = 0.2) +
    geom_point(data = as.data.frame(sp$samplers), aes(x,y, size = n), color = "red")
  
#' Run Poisson regression
#' 
#' The two terms of the formula's left hand side define the number of counts modeled
#' and the exposure, which is the area of each point sample. This data is available from the
#' `sp$samplers` data frame. The `.` on the right hand side indicates that the default model
#' should be use, i.e. SPDE + Intercept.
#' 
#+results="hide"
  
  r = poiss(sp$samplers, n + weight ~ ., mesh = sp$mesh)
 
#' Plot intensity
#' 
#' TO DO: prediction method for poisson regression results not yet added to package

  plot.spatial(r) +
    gg.point(sp$samplers, size = 5, color = "black", alpha = 0.2) +
    geom_point(data = as.data.frame(sp$samplers), aes(x,y, size = n), color = "red")
  