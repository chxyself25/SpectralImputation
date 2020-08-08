# interpolation on spatial direction, same footprint
# used for imputation for missing values
# first interpolation in spatial direction, then interpolation in wavelength

SpatialFill <- function(y, long, lat) {
  na.idx <- which(is.na(y))
  nna.idx <- which(!is.na(y))
  intp <- interp(x=long[nna.idx], y=lat[nna.idx], z=y[nna.idx], xo=long[na.idx], yo=lat[na.idx],
                 linear = TRUE, duplicate = 'mean', extrap = TRUE)
  y[na.idx] <- diag(intp$z)
  return(y)
}

WavesFill <- function(y, waves) {
  na.idx <- which(is.na(y))
  nna.idx <- which(!is.na(y))
  intp <- approx(x = waves[nna.idx], y = y[nna.idx], xout = waves[na.idx], method = 'linear', rule = 2)
  y[na.idx] <- intp$y
  return(y)
}

NAFill <- function(ymat, long, lat, waves) {
  # spatial bivariate interpolation
  ymat1 <- apply(ymat, 2, function(yy) {
    if(any(is.na(yy))) {
      SpatialFill(yy, long, lat)
    } else {yy}
  })
  if (any(is.na(ymat1))) {
    # interpolation in wavelength direction
    ymat2 <- t(apply(ymat1, 1, function(yy) {
      if (any(is.na(yy))) {
        WavesFill(yy, waves)
      } else {yy}
    }))
  } else {
    ymat2 <- ymat1
  }
  return(ymat2)
}

