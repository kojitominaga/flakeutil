\name{hurr2hurs}
\alias{hurr2hurs}
\title{Calculate specific humidity kg kg-1}
\description{
  Calculates air specific humidity from relative humidity, 
  air temperature and air pressure.
}
\usage{
hurr2hurs(hurr, ta, ps)
}
\arguments{
  \item{hurr}{\code{numeric} \code{vector} for relative humidity in percent}
  \item{ta}{\code{numeric} \code{vector} for air temperture in K}
  \item{ps}{\code{numeric} \code{vector} for air pressure in Pa}
}
\value{\code{numeric} \code{vector} for specific humidity in kg kg-1}}
\references{\url{http://www.cactus2000.de/js/calchum.pdf}}


