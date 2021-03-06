% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/teststatistics.R
\name{get.crit.FB}
\alias{get.crit.FB}
\title{Critical values for the fixed-b unit root test}
\usage{
get.crit.FB(b)
}
\arguments{
\item{b}{Relative blocklength.}
}
\value{
A list containing the following components:
\item{crit.values}{A vector of critical values for different significance levels.}
\item{rel.blocklength}{The relative blocklength for which the critical values are reported. If the critical values for the input value of b are not available, the nearest available value is considered.}
}
\description{
Provides critical values for the fixed-b unit root test for a given relative blocklength for different significance levels.
}
\examples{
get.crit.FB(0.2)
get.crit.FB(0.48)
}
