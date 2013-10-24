\name{TemperatureAtDepths}
\alias{TemperatureAtDepths}
\title{Calculate Water Temperature at a Given Series of Depths}
\description{
  Calculates water temperature at given depths based on a results from
  \code{runflake}.
}
\usage{
TemperatureAtDepths(flakeresult, parameters, z)
}
\arguments{
  \item{flakeresult}{the \code{data.frame} object from
    \code{flaker::runflake}}  
  \item{parameters}{the named numeric \code{vector}, which was a part of
    the input to \code{flaker::runflake}}
  \item{z}{the depths in metre at which temperature is calculated}
}
\value{
  a numeric matrix of size (no. time steps x depths) in Kelvin
}


