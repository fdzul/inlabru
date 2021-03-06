% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mesh.R
\name{is.inside}
\alias{is.inside}
\title{Query if a point is inside the mesh boundary}
\usage{
is.inside(mesh, loc, mesh.coords = NULL)
}
\arguments{
\item{mesh}{an inla.mesh object.}

\item{loc}{Points in space stored either as data.frame, a two-column matrix of x and y coordinates or a SpatialPoints object.}

\item{mesh.coords}{Coordinate names of the mesh. Use only if loc is a data.frame with respective column names.}
}
\value{
Single column matrix of Boolean values indicating if a point is inside the mesh.
}
\description{
Query if a point is inside the mesh boundary
}
\examples{
\dontrun{
# Load Gorilla data

data("gorillas")

# Check if all Gorilla nests are inside the mesh

all(is.inside(gorillas$mesh, gorillas$nests))

# Also works for locations not stored as SpatialPoints object

loc = coordinates(gorillas$nests)
all(is.inside(gorillas$mesh, loc))
}

}
\author{
Fabian E. Bachl <\email{bachlfab@gmail.com}>
}
\keyword{internal}
