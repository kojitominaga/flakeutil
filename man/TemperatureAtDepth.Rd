\name{TemperatureAtDepth}
\alias{TemperatureAtDepth}
\title{Calculate Water Temperature at a Given Depth}
\description{
  Calculates water temperature at a given depth based on a results from
  \code{runflake}.
}
\usage{
TemperatureAtDepth(flakeresult, parameters, z)
}
\arguments{
  \item{flakeresult}{the \code{data.frame} object from
    \code{flaker::runflake}}  
  \item{parameters}{the named numeric \code{vector}, which was a part of
    the input to \code{flaker::runflake}}
  \item{z}{the depth in metre at which temperature is calculated}
}
\value{
  a numeric vector of length \code{nrow(flakeresult)} in Kelvin
}


