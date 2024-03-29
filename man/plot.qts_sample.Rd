% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qts-sample-class.R
\name{plot.qts_sample}
\alias{plot.qts_sample}
\title{Plot for \code{\link{qts_sample}} objects}
\usage{
\method{plot}{qts_sample}(x, memberships = NULL, highlighted = NULL, with_animation = FALSE, ...)
}
\arguments{
\item{x}{An object of class \code{\link{qts_sample}}.}

\item{memberships}{A vector coercible as factor specifying a group membership
for each QTS in the sample. Defaults to \code{NULL}, in which case no grouping
structure is displayed.}

\item{highlighted}{A boolean vector specifying whether each QTS in the sample
should be hightlighted. Defaults to \code{NULL}, in which case no QTS is
hightlighted w.r.t. the others.}

\item{with_animation}{A boolean value specifying whether to create a an
animated plot or a static \link[ggplot2:ggplot]{ggplot2::ggplot} object. Defaults to \code{FALSE}
which will create a static plot.}

\item{...}{Further arguments to be passed to methods.}
}
\value{
No return value, called for side effects.
}
\description{
This function creates a visualization of a sample of QTS \strong{without}
returning the corresponding \link[ggplot2:ggplot]{ggplot2::ggplot} object
}
\examples{
plot(vespa64$igp)
}
