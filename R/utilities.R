PlotSolarRadiation <- function(flakeresult) {
  Iatm <- flakeresult[['I_atm_flk_out']]
  Isnow <- flakeresult[['I_snow_flk_out']]
  Iice <- flakeresult[['I_ice_flk_out']]
  Iw <- flakeresult[['I_w_flk_out']]
  Ih <- flakeresult[['I_h_flk_out']]
  Ibot <- flakeresult[['I_bot_flk_out']]
  Asnow <- flakeresult[['I_refl_snow']]
  Aice <- flakeresult[['I_refl_ice']]
  Aw <- flakeresult[['I_refl_w']]
  AttenSnow <- Isnow - Iice - Aice
  AttenIce <- Iice - Iw - Aw
  AttenML <- Iw - Ih
  AttenTC <- Ih - Ibot
  print(round(apply(data.frame(Asnow, Aice, Aw,
                               AttenSnow, AttenIce, AttenML, AttenTC,
                               Ibot),
                    2, summary
                    ),
              digits = 1
              )
        )
  x <- c(1:nrow(flakeresult), rev(1:nrow(flakeresult)))
  par(bg = 'darkgray')
  plot(0, 0, type = 'n',
       xlim = range(x), ylim = range(c(Iatm, 0)),
       xlab = 'time step', ylab = 'fate of incoming global radiation (W m-2)',
       ...
       )
  abline(h = 0, col = 'gray')
  l0 <- rep(0, times = nrow(flakeresult))
  l1 <- l0 + Ibot
  l2 <- l1 + AttenTC
  l3 <- l2 + AttenML
  l4 <- l3 + Aw
  l5 <- l4 + AttenIce
  l6 <- l5 + Aice
  l7 <- l6 + AttenSnow
  l8 <- l7 + Asnow
  polygon(x, c(l0, rev(l1)), border = NA, col = 'gray')      # Ibot
  polygon(x, c(l1, rev(l2)), border = NA, col = 'darkblue')  # AttenTC
  polygon(x, c(l2, rev(l3)), border = NA, col = 'blue')      # AttenML
  polygon(x, c(l3, rev(l4)), border = NA, col = 'lightblue') # Awater
  polygon(x, c(l4, rev(l5)), border = NA, col = 'red')       # AttenIce
  polygon(x, c(l5, rev(l6)), border = NA, col = 'pink')      # Aice
  polygon(x, c(l6, rev(l7)), border = NA, col = 'lightgray') # AttenSnow
  polygon(x, c(l7, rev(l8)), border = NA, col = 'white')     # Asnow
  legend('topleft',
         legend = c(
           'reflected on snow',
           'reflected on ice',
           'reflected on water',
           'absorbed in snow',
           'absorbed in ice',
           'absorbed in mixed layer',
           'absorbed in thermocline (hypolimnion)',
           'absorbed in sediment'),
         fill = c('white', 'pink', 'lightblue', 'lightgray', 'red', 'blue',
           'darkblue', 'gray'),
         border = NA,
         ncol = 2,
         bty = 'n',
         cex = 0.8
         )
  
  plot(0, 0, type = 'n',
       xlim = range(x), ylim = c(0, 1.25),
       xlab = 'time step',
       ylab = 'fate of incoming global radiation (fraction)',
       ...
       )
  abline(h = 0, col = 'gray')
  l0 <- rep(0, times = nrow(flakeresult))
  l1 <- l0 + Ibot
  l2 <- l1 + AttenTC
  l3 <- l2 + AttenML
  l4 <- l3 + Aw
  l5 <- l4 + AttenIce
  l6 <- l5 + Aice
  l7 <- l6 + AttenSnow
  l8 <- l7 + Asnow
  l0 <- l0 / Iatm
  l1 <- l1 / Iatm
  l2 <- l2 / Iatm
  l3 <- l3 / Iatm
  l4 <- l4 / Iatm
  l5 <- l5 / Iatm
  l6 <- l6 / Iatm
  l7 <- l7 / Iatm
  l8 <- l8 / Iatm
  polygon(x, c(l0, rev(l1)), border = NA, col = 'gray')      # Ibot
  polygon(x, c(l1, rev(l2)), border = NA, col = 'darkblue')  # AttenTC
  polygon(x, c(l2, rev(l3)), border = NA, col = 'blue')      # AttenML
  polygon(x, c(l3, rev(l4)), border = NA, col = 'lightblue') # Awater
  polygon(x, c(l4, rev(l5)), border = NA, col = 'red')       # AttenIce
  polygon(x, c(l5, rev(l6)), border = NA, col = 'pink')      # Aice
  polygon(x, c(l6, rev(l7)), border = NA, col = 'lightgray') # AttenSnow
  polygon(x, c(l7, rev(l8)), border = NA, col = 'white')     # Asnow
  legend('topleft',
         legend = c(
           'reflected on snow',
           'reflected on ice',
           'reflected on water',
           'absorbed in snow',
           'absorbed in ice',
           'absorbed in mixed layer',
           'absorbed in thermocline (hypolimnion)',
           'absorbed in sediment'),
         fill = c('white', 'pink', 'lightblue', 'lightgray', 'red', 'blue',
           'darkblue', 'gray'),
         border = NA,
         ncol = 2,
         bty = 'n',
         cex = 0.8
         )
}

FractionalRadiation <- function(d, lat, lon, n = 24) {
  require(oce)
  midnight <- as.POSIXct(floor(as.numeric(as.POSIXct(d)) / 86400) * 86400,
                         origin = '1970-01-01', tz = 'UTC'
                         )
  ## bounds <- seq(from = midnight,
  ##               to = midnight + 60 * 60 * 24,
  ##               length.out = n + 1)
  middletimes <- seq(from = midnight + 60 * 60 * 24 / n / 2,
                     to = midnight + 60 * 60 * 24 * (n - 0.5) / n,
                     length.out = n)
  intensity <- sin(pmax(sunAngle(middletimes, lat, lon)[['altitude']], 0)
                   / 180 * pi)
  ## below fixes for the polar winter without sun
  if (sum(intensity) == 0) {
    intensity <- rep(1 / length(intensity), times = length(intensity))
  }
  fraction <- intensity / sum(intensity)
  return(fraction)
}

## fraction <- fractionradiation(as.Date('2013-12-09'), 59, 11)
## meanI <- 200 # W m-2
## hourlyI <- meanI * n * fraction

TemperatureAtDepth <- function(flakeresult, parameters, z) {
  ## flakeresult should be data.frame with certain names.
  ## parameters should be named vector. need only for lake depth.
  ## z is the depth at which temperature is wanted.
  ## See runflake.R for details
  Ts <- flakeresult[['T_wML']] # mixed layer (top) temperature
  Tb <- flakeresult[['T_bot']] # bottom temperature 
  h <- flakeresult[['h_ML']] # mixed layer depth
  C <- flakeresult[['C_T']] # shape factor
  D <- parameters[['depth_w']]
  zeta <- (z - h) / (D - h)
  ## eq. 55
  c1 <- 40 / 3
  c2 <- 20 / 3
  c3 <- 5 / 3
  c4 <- 10 / 3
  is.in.ML <- z <= h
  T <- ifelse(is.in.ML,
              Ts,
              zeta * (c1 * C - c2 + zeta * (18 - 30 * C +
                zeta * (20 * C - 12 + zeta * (c3 - c4 * C)))) *
              (Tb - Ts) + Ts
              )
}

TemperatureAtDepths <- function(flakeresult, parameters, z) {
  ## flakeresult should be data.frame with certain names.
  ## parameters should be named vector. need only for lake depth.
  ## z is the depth at which temperature is wanted. NOW VECTOR
  ## See runflake.R for details
  ## value: matrix first dimension time step, second dimension z
  Ts <- flakeresult[['T_wML']] # mixed layer (top) temperature
  Tb <- flakeresult[['T_bot']] # bottom temperature 
  h <- flakeresult[['h_ML']] # mixed layer depth
  C <- flakeresult[['C_T']] # shape factor
  D <- parameters[['depth_w']]
  zeta <- apply(matrix(1:length(z), ncol = length(z)),
                2,
                function(zi) (z[zi] - h) / (D - h)
                )
  ## eq. 55
  c1 <- 40 / 3
  c2 <- 20 / 3
  c3 <- 5 / 3
  c4 <- 10 / 3
  is.in.ML <- apply(matrix(1:length(z), ncol = length(z)),
                    2,
                    function(zi) z[zi] <= h
                    )
  T <- apply(matrix(1:length(z), ncol = length(z)),
             2,
             function(zi) {
               ifelse(is.in.ML[ , zi],
                      Ts,
                      zeta[ , zi] * (c1 * C - c2 + zeta[ , zi] *
                                     (18 - 30 * C + zeta[ , zi] *
                                      (20 * C - 12 + zeta[ , zi] *
                                       (c3 - c4 * C)
                                       )
                                      )
                                     ) *
                      (Tb - Ts) + Ts
                      )
             }
             )
}
                
HourlyAirTemperature <- function(d, lat, lon,
                                 startT, minT, maxT, meanT, minTn, maxTn,
                                 n = 24) {
  ## decided to handle errors in running scripts
  ## !! TODO (see commit details on 2013-10-27)
  require(oce)
  threshold <- ifelse(meanT < 273, 0.38, 0.30)
  if (!(n = 24)) {
    stop('not sure if it works with n being not 24')
  }
  T <- numeric(n)
  midnight <- as.POSIXct(floor(as.numeric(as.POSIXct(d)) / 86400) * 86400,
                         origin = '1970-01-01', tz = 'UTC'
                         )
  ## bounds <- seq(from = midnight,
  ##               to = midnight + 60 * 60 * 24,
  ##               length.out = n + 1)
  middletimes <- seq(from = midnight + 60 * 60 * 24 / n / 2,
                     to = midnight + 60 * 60 * 24 * (n - 0.5) / n,
                     length.out = n)
  intensity <- sin(pmax(sunAngle(middletimes, lat, lon)[['altitude']], 0)
                   / 180 * pi)
  sunupi <- min(which(intensity > 0))
  timestart <- 0 # hour
  timeminT <- 0.5 + 24 / n * (sunupi - 1) + 1
  ## hour, an hour later first middle time
  timemaxT <- 14 + 24 / n / 2
  ## hour, the middle time after 14:00
  mt <- seq(from = 24 / n / 2,
            to = 24 - 24 / n / 2,
            length.out = n)
  if (abs(maxT - minTn) <= 0.5) {
    if (startT < minT) {
      ## print('case 1.1')
      ## temperature increased throughout the 24 hours, probably.
      ## assume triangle shape at 12:30
      noonT <- 2 * meanT - 0.5 * minT - 0.5 * maxT
      for (tsi in 1:(n / 2)){
        T[tsi] <- minT + (noonT - minT) * (mt[tsi] - 1) / 11.5
      }
      for (tsi in (n / 2 + 1):n) {
        T[tsi] <- noonT + (maxT - noonT) * (mt[tsi] - 12.5) / 11.5
      }
      ## adjustment at the middle 12 hours
      T[7:18] <- T[7:18] - (sum(T) - meanT * n) / 12
    } else {
      ## print('case 1.2')
      ## temperature goes down first and then up (end increasing)
      tsl1 <- mt <= timeminT
      tsii1 <- c(1:n)[tsl1]
      for (tsi in tsii1) {
        T[tsi] <- startT + (startT - minT) * sin(pi / 4) / (1 - sin(pi / 4)) -
          sin(mt[tsi] / timeminT * (pi / 4) + pi / 4) *
            (startT - minT) / (1 - sin(pi / 4))
      }
      rest <- meanT * n - sum(T)
      ## let x be noonT (unknown) minus minT
      ## then fx(time), gx(time) are scaling functions s.t.
      ## fx * x + gx * (maxT - x - minT) +  minT = T
      fx <- numeric()
      gx <- numeric()
      for (tsi in tsii1) {
        ## zero because it's already taken into account befor calculating rest
        fx[tsi] <- 0
        gx[tsi] <- 0
      }
      tsl2 <- (!tsl1) & (mt < 12)
      tsii2 <- c(1:n)[tsl2]
      for (tsi in tsii2) {
        fx[tsi] <- 1 +
          sin((mt[tsi] - timeminT) / (12 - timeminT) * pi / 2 - pi / 2)
        gx[tsi] <- 0
      }
      tsl3 <- (!tsl1) & (!tsl2)
      tsii3 <- c(1:n)[tsl3]
      for (tsi in tsii3) {
        fx[tsi] <- 1
        gx[tsi] <- (mt[tsi] - 12) / 11.5
      }
      ## now
      ## sum(fx) * x + sum(gx) * (maxT - x - minT) + minT * sum(tsl2 | tsl3)
      ##  == rest
      x <- (rest - minT * sum(tsl2 | tsl3) - sum(gx) * (maxT - minT)) /
        (sum(fx) - sum(gx))
      for (tsi in tsii2) {
        T[tsi] <- fx[tsi] * x + gx[tsi] * (maxT - x - minT) + minT
      }
      for (tsi in tsii3) {
        T[tsi] <- fx[tsi] * x + gx[tsi] * (maxT - x - minT) + minT
      }
    }
  } else if (abs(minT - maxTn) <= 0.5) {
    if (startT > maxT) {
      ## print('case 2.1')
      ## temperature decreased throughout the 24 hours, probably.
      ## assume triangle shape at 12:30
      noonT <- 2 * meanT - 0.5 * minT - 0.5 * maxT
      for (tsi in 1:(n / 2)){
        T[tsi] <- maxT - (maxT - noonT) * (mt[tsi] - 0.5) / 11.5
      }
      for (tsi in (n / 2 + 1):n) {
        T[tsi] <- noonT - (noonT - minT) * (mt[tsi] - 12) / 11.5
      }
      ## adjustment at the middle 12 hours
      T[7:18] <- T[7:18] - (sum(T) - meanT * n) / 12
    } else {
      ## print('case 2.2')
      ## temperature goes up first and then down (and decreasing)
      ## assume timemaxT at 10:00 (need to make room to adjust in the afternoon)
      timemaxT <- 10
      tsl1 <- mt <= timemaxT
      tsii1 <- c(1:n)[tsl1]
      for (tsi in tsii1) {
        T[tsi] <- startT + (maxT - startT) *
          sin((mt[tsi] / timemaxT) * (pi / 2))
      }
      rest <- meanT * n - sum(T)
      ## let x be T at 18:00 (unknown) minus minT
      ## then fx(time), gx(time) are scaling functions s.t.
      ## fx * x + gx * (maxT - x - minT) +  minT = T
      fx <- numeric()
      gx <- numeric()
      for (tsi in tsii1) {
        ## zero because it's already taken into account befor calculating rest
        fx[tsi] <- 0
        gx[tsi] <- 0
      }
      tsl2 <- (!tsl1) & (mt < 18)
      tsii2 <- c(1:n)[tsl2]
      for (tsi in tsii2) {
        fx[tsi] <- 1 
        gx[tsi] <- sin((mt[tsi] - timemaxT) / (18 - timemaxT) * (pi / 2) +
                       pi / 2)
      }
      tsl3 <- (!tsl1) & (!tsl2)
      tsii3 <- c(1:n)[tsl3]
      for (tsi in tsii3) {
        fx[tsi] <- 1 - (mt[tsi] - 18) / 5.5
        gx[tsi] <- 0
      }
      ## now
      ## sum(fx) * x + sum(gx) * (maxT - x - minT) + minT * sum(tsl2 | tsl3)
      ##  == rest
      x <- (rest - minT * sum(tsl2 | tsl3) - sum(gx) * (maxT - minT)) /
        (sum(fx) - sum(gx))
      for (tsi in tsii2) {
        T[tsi] <- fx[tsi] * x + gx[tsi] * (maxT - x - minT) + minT
      }
      for (tsi in tsii3) {
        T[tsi] <- fx[tsi] * x + gx[tsi] * (maxT - x - minT) + minT
      }
    }
  } else {
    if (((meanT - minT) / (maxT - minT)) < threshold) {
      ## print('case 3.1.1')
      ## temperature increased throughout the 24 hours, probably.
      ## assume triangle shape at 12:30
      noonT <- 2 * meanT - 0.5 * minT - 0.5 * maxT
      for (tsi in 1:(n / 2)){
        T[tsi] <- minT + (noonT - minT) * (mt[tsi] - 1) / 11.5
      }
      for (tsi in (n / 2 + 1):n) {
        T[tsi] <- noonT + (maxT - noonT) * (mt[tsi] - 12.5) / 11.5
      }
      ## adjustment at the middle 12 hours
      T[7:18] <- T[7:18] - (sum(T) - meanT * n) / 12
    } else if (((maxT - meanT) / (maxT - minT)) < threshold) {
      ## print('case 3.1.2')
      ## temperature decreased throughout the 24 hours, probably.
      ## assume triangle shape at 12:30
      noonT <- 2 * meanT - 0.5 * minT - 0.5 * maxT
      for (tsi in 1:(n / 2)){
        T[tsi] <- maxT - (maxT - noonT) * (mt[tsi] - 0.5) / 11.5
      }
      for (tsi in (n / 2 + 1):n) {
        T[tsi] <- noonT - (noonT - minT) * (mt[tsi] - 12) / 11.5
      }
      ## adjustment at the middle 12 hours
      T[7:18] <- T[7:18] - (sum(T) - meanT * n) / 12
    } else {
      if (startT < minT) {
        ## print('case 3.2.1')
        startT <- minT
        timeminT <- 5
      } else if (startT > maxT) {
        ## print('case 3.2.2')
        startT <- maxT
        timeminT <- 10
      } else if ((((meanT - minT) / (maxT - minT)) < (threshold + 0.1)) &
                 (((startT - maxT) / (maxT - minT)) < (threshold + 0.1))) {
        ## print('case 3.2.3')
        startT <- meanT
        timeminT <- 5
      } else if ((((maxT - meanT) / (maxT - minT)) < (threshold + 0.1)) &
                 (((startT - minT) / (maxT - minT)) < (threshold + 0.1))) {
        ## print('case 3.2.4')
        startT <- meanT
        timeminT <- 5
      } else {
        ## print('case 3.2.5')
        ## the most ordinary case, down, up and down
      }
      tsl1 <- mt <= timeminT
      tsii1 <- c(1:n)[tsl1]
      for (tsi in tsii1) {
        T[tsi] <- startT + (startT - minT) * sin(pi / 4) / (1 - sin(pi / 4)) -
          sin(mt[tsi] / timeminT * (pi / 4) + pi / 4) *
            (startT - minT) / (1 - sin(pi / 4))
      }
      tsl2 <- (!tsl1) & (mt <= timemaxT)
      tsii2 <- c(1:n)[tsl2]
      for (tsi in tsii2) {
        T[tsi] <-
          (minT + maxT) / 2 +
            sin((mt[tsi] - timeminT) / (timemaxT - timeminT) * pi - pi / 2) *
              (maxT - minT) / 2
      }
      
      tsl3 <- !(tsl1 | tsl2)
      tsii3 <- c(1:n)[tsl3]
      ## calculate the rest of the temperature sum
      rest <- meanT * n - sum(T)
      if (rest < 0) {
        stop("")
      }
      ## !!! check that the min stays minimum and max stays max
      ## calculate how to distribute the rest
      frac <- rep(0, length = 24) # this would be the fraction above the endT
      for (tsi in tsii3) {
        frac[tsi] <- sin((mt[tsi] - timemaxT) / (24 - timemaxT) * (pi * 3 / 4) +
                         pi / 2)
      }
      endT <- (rest - maxT * sum(frac)) / (sum(tsl3) - sum(frac))
      for (tsi in tsii3) {
        T[tsi] <- endT + (maxT - endT) * frac[tsi]
      }
    }
  }
  ## last check
  allowance <- 2
  if (T[n] < (minT - allowance)) {
    Tnew <- ifelse(T < minT, minT, T)
    sumdiff <- sum(Tnew) - sum(T)
    Tnew[7:18] <- Tnew[7:18] - c(0.02, 0.03, 0.06, 0.09, 0.13, 0.17,
                                 0.17, 0.13, 0.09, 0.06, 0.03, 0.02) * sumdiff
    Told <- T
    T <- Tnew
  } else if (T[n] > (maxT + allowance)) {
    Tnew <- ifelse(T > maxT, maxT, T)
    sumdiff <- sum(T) - sum(Tnew)
    Tnew[7:18] <- Tnew[7:18] + c(0.02, 0.03, 0.06, 0.09, 0.13, 0.17,
                                 0.17, 0.13, 0.09, 0.06, 0.03, 0.02) * sumdiff
    T <- Tnew
  }
  return(T)
}

LAIndices <- function(flakeresult, parameters, bthA, bthD) {
  require(rLakeAnalyzer)
  if (length(bthD) != length(bthA)) stop('bthA and D must have the same length')
  if (max(bthD) != parameters[['depth_w']]) stop('max(bthD) must be the depth')
  if (min(bthD) != 0) stop('min(bthD) must be zero')
  if (rev(bthA)[1] != 0) stop('lake bottom must have zero area')
  cat('This function uses rLakeAnalyzer functions.\n')
  cat('See the documentation for this function\n')
  cat('and references for rLakeAnalyzer\n')
  depths <- bthD
  cat('... calculating temperature at depths\n')
  tzm <- TemperatureAtDepths(data.frame(flakeresult), parameters, depths)
  cat('... calculating depth broundaries for metalimnion\n')
  md <- apply(tzm - 273.15, 1, function(x) meta.depths(x, depths))
  cat('... calculating thermocline depth\n')
  td <- apply(tzm - 273.15, 1, function(x) thermo.depth(x, depths))
  cat("... calculating Schmidt's stability\n")
  ss <- apply(tzm - 273.15, 1,
              function(x) schmidt.stability(x, depths, bthA, depths))
  cat('... calculating mean densities for epi- and hypolimnion\n') 
  ed <- numeric()
  hd <- numeric()
  for (ti in 1:nrow(flakeresult)) {
    ed[ti] <- layer.density(0, md[1, ti], tzm[ti, ] - 273.15,
                            depths, bthA, depths)
    hd[ti] <- layer.density(md[2, ti],
                            parameters[['depth_w']], tzm[ti, ] - 273.15,
                            depths, bthA, depths)
  }
  cat('... calculating u*, friction velocity\n')
  us <- numeric()
  uain <- data.frame(flakeresult)[['U_a_in']]
  for (ti in 1:nrow(flakeresult)) {
    us[ti] <- uStar(uain[ti], 10, ed[ti]) 
  }
  cat('... calculating lake number\n')
  ln <- lake.number(bthA, depths, us, ss, md[1, ], md[2, ], hd)
  cat('... calculating Wedderburn number\n')
  wn <- wedderburn.number(hd - ed, md[1, ], us, bthA[1], hd)
  out <- data.frame(ThermoclineDepth = td,
                    MixedLayerDepth = md[1, ],
                    SchmidtStability = ss,
                    WedderburnNumber = wn,
                    LakeNumber = ln)
}
  

