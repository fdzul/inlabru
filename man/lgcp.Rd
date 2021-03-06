% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bru.inference.R
\name{lgcp}
\alias{lgcp}
\title{Log Gaussian Cox process (LGCP) inference using INLA}
\usage{
lgcp(components, data, samplers = NULL, domain = NULL, ips = NULL,
  formula = . ~ ., E = 1, options = list())
}
\arguments{
\item{components}{A formula describing the latent components}

\item{data}{A data frame or SpatialPoints[DataFrame] object}

\item{samplers}{A data frame or Spatial[Points/Lines/Polygons]DataFrame objects}

\item{domain}{Named list of domain definitions}

\item{ips}{Integration points (overrides \code{samplers})}

\item{formula}{If NULL, the linear combination implied by the \code{components} is used as a predictor for the point location intensity. If a (possibly non-linear) expression is provided the respective Taylor approximation is used as a predictor. Multiple runs if INLA are then required for a better approximation of the posterior.}

\item{E}{Single numeric used rescale all integration weights by a fixed factor}

\item{options}{See \link{bru.options}}
}
\value{
An \link{bru} object
}
\description{
This function performs inference on a LGCP observed via points residing possibly multiple dimensions. 
These dimensions are defined via the left hand side of the formula provided via the model parameter.
The left hand side determines the intensity function that is assumed to drive the LGCP. This may include
effects that lead to a thinning (filtering) of the point process. By default, the log intensity is assumed
to be a linear combination of the effects defined by the formula's RHS. More sofisticated models, e.g.
non-linear thinning, can be achieved by using the predictor argument. The latter requires multiple runs
of INLA for improving the required approximation of the predictor. In many applications
the LGCP is only observed through subsets of the dimensions the process is living in. For example, spatial
point realizations may only be known in sub-areas of the modeled space. These observed subsets of the LGCP
domain are called samplers and can be provided via the respective parameter. If samplers is NULL it is
assumed that all of the LGCP's dimensions have been observed completely.
}
\examples{

\dontrun{

Load the Gorilla data
data(gorillas)

# Use tutorial setting and thus empirical Bayes for faster inference
init.tutorial()

# Plot the Gorilla nests, the mesh and the survey boundary
ggplot() + 
  gg(gorillas$mesh) + 
  gg(gorillas$nests) + 
  gg(gorillas$boundary) + 
  coord_fixed()

# Define SPDE prior
matern <- inla.spde2.pcmatern(gorillas$mesh, 
                              prior.sigma = c(0.1, 0.01), 
                              prior.range = c(5, 0.01))

# Define domain of the LGCP as well as the model components (spatial SPDE effect and Intercept)
cmp <- coordinates ~ mySmooth(map = coordinates, model = matern) + Intercept

# Fit the model
fit <- lgcp(cmp, gorillas$nests, samplers = gorillas$boundary)

# Predict the spatial intensity surface
lambda <- predict(fit, pixels(gorillas$mesh), ~ exp(mySmooth + Intercept))

# Plot the intensity
ggplot() + 
  gg(lambda) +
  gg(gorillas$mesh) + 
  gg(gorillas$nests) + 
  gg(gorillas$boundary) + 
  coord_fixed()

}

}
