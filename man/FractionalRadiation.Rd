\name{FractionalRadiation}
\alias{FractionalRadiation}
\title{24-hourize (or more or less) the radiation based on Julian day}
\description{
  Calculates fractional allocation of the radiation for the sub-daily
  timesteps based on the day of year and geographical location, returns a
  vector. The values add up to 1.
}
\usage{
FractionalRadiation(d, lat, lon, n = 24)
}
\arguments{
  \item{d}{a \code{Date} object, on which sub-daily allocation of
    radiation is calculated}
  \item{lat}{latitude in degrees, positive north}
  \item{lon}{longitude in degrees, positive east}
  \item{n}{an integer for the number of division in a day, defaults to 24}
}
\value{
  a numeric \code{vector} of length \code{n} that adds up to 1.
}
%% \seealso{
%%   \code{\link{sunAngle}}, \code{\link{as.Date}}.
%% }

