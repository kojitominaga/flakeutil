\name{HourlyAirTemperature}
\alias{HourlyAirTemperature}
\title{Dirty 24-hourization of Air Temperature from Daily Statistics}
\description{
  A very unscientific guess of what might have happened to air
  temperature during the 24 hours, based on daily mean TA, daily maximum
  TA, and daily minimum TA. It also vaguely considers the starting TA
  that was "inherited" from the previous day, and known information
  about the minimum TA and maximum TA on the following day.
}
\usage{
HourlyAirTemperature(d, lat, lon,
                     startT, minT, maxT, meanT, minTn,
                     maxTn, n = 24)
}
\arguments{
  \item{d}{a \code{Date} object, on which sub-daily allocation of
    radiation is calculated}
  \item{lat}{latitude in degrees, positive north}
  \item{lon}{longitude in degrees, positive east}
  \item{startT}{the starting air temperature in Kelvin that was
    inherited from the previous day}
  \item{minT}{the daily mininum air temperature in Kelvin of the 
    concerned day}
  \item{maxT}{the daily maximum air temperature in Kelvin of the
    concerned day}
  \item{meanT}{the daily mean air temperature in Kelvin of the concerned
    day}
  \item{minTn}{the daily mininum air temperature in Kelvin of the
    following day}
  \item{maxTn}{the daily maximum air temperature in Kelvin of the
    following day}
  \item{n}{an integer for the number of division in a day, defaults to
    24}
}
\value{
  a numeric \code{vector} of length \code{n}, crude guesses for
  sub-daily temperature in Kelvin
}


  
