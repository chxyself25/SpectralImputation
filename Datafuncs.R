## functions for organizing dataset, prediction on xi, imputation and land fraction estimation

library(ncdf4)

## function for extracting variables from data downloaded in NASA data center
# file is the directory of h5 file, radiance output is devided by 10^19
GetdataNASA <- function(file, save.path = NULL) {
  options(scipen = 999)
  nc1 <- nc_open(file)
  lat <- ncvar_get(nc1, "RetrievalGeometry/retrieval_latitude")
  long <- ncvar_get(nc1, "RetrievalGeometry/retrieval_longitude")
  land_fraction <- ncvar_get(nc1, "RetrievalGeometry/retrieval_land_fraction")
  sound_id <- as.character(ncvar_get(nc1, "RetrievalHeader/sounding_id"))
  ftprint <- as.numeric(substr(sound_id, start = 16, stop = 16)) # the last digit is the footprint
  track <- substr(sound_id, start = 9, stop = 15)
  f <- strsplit(file, split = "/")[[1]]
  orbit <- substr(strsplit(f[length(f)], split = "_")[[1]][3], 1, 5)
  #solar_az <- ncvar_get(nc1, "RetrievalGeometry/retrieval_solar_azimuth")
  #solar_zen <- ncvar_get(nc1, "RetrievalGeometry/retrieval_solar_zenith")
  locs <- data.frame(id = sound_id, orbit = orbit, long = long, lat = lat, ftprint = ftprint, 
                     mix = land_fraction, track = track, stringsAsFactors = FALSE)
  row.names(locs) <- NULL
  meas_rad <- t(ncvar_get(nc1, "SpectralParameters/measured_radiance"))
  waves <- t(ncvar_get(nc1, "SpectralParameters/sample_indexes"))
  wavelength <- t(ncvar_get(nc1, "SpectralParameters/wavelength"))
  num_colors <- ncvar_get(nc1, "SpectralParameters/num_colors_per_band")[1, ] # only O2 band
  meas_rad_o2 <- foreach(i = 1:nrow(meas_rad)) %dopar% {
    meas_rad[i, 1:num_colors[i]]/(10^19)
  }
  waves_o2 <- foreach(i = 1:nrow(waves)) %dopar% {
    waves[i, 1:num_colors[i]]
  }
  wavelength_o2 <- foreach(i = 1:nrow(wavelength)) %dopar% {
    wavelength[i, 1:num_colors[i]]
  }
  options(scipen = 0)
  if (!is.null(save.path)) {
    saveRDS(locs, file = paste0(save.path, "/nasa_", orbit, "_locs.rds"))
    saveRDS(meas_rad_o2, file = paste0(save.path, "/nasa_", orbit, "_rads.rds"))
    saveRDS(waves_o2, file = paste0(save.path, "/nasa_", orbit, "_waves.rds"))
  }
  return(list(locations = locs, radiance = meas_rad_o2, waves = waves_o2, wavelength = wavelength_o2))
}

# library(doParallel); registerDoParallel(cores = 8)
# file <- "./pacific_data/oco2_L2DiaGL_14793a_170413_B8100r_170910085430.h5"
# file <- "./pacific_data/oco2_L2DiaGL_15985a_170703_B8100r_171018175833.h5"
# file <- "./pacific_data/oco2_L2DiaGL_17123a_170920_B8100r_171013011336.h5"
# res <- GetdataNASA(file, save = "./pacific_data/")


## function for getting variables from Coastal provided by Jon in JPL
# not fully tested, under construction
Getdata <- function(file, num.threshold = NULL, per.threshold = 1) {
  options(scipen = 999)
  nc1 <- nc_open(file)
  lat <- ncvar_get(nc1, "latitude")
  long <- ncvar_get(nc1, "longitude")
  id <- as.character(ncvar_get(nc1, "sounding_id"))
  measured <- ncvar_get(nc1, "measured_radiance")
  wavelength <- ncvar_get(nc1, "wavelength")
  ftprint <- ncvar_get(nc1, "footprint")
  orbit <- ncvar_get(nc1, "orbit")
  land_frac <- ncvar_get(nc1, "land_fraction")
  locs <- data.frame(long = long, lat = lat, id = id, ftprint = ftprint, orbit = orbit, mix = land_frac, stringsAsFactors = FALSE)
  rads <- t(measured)[,1:1016] #only O2 band
  wavelength <- t(wavelength)[,1:1016]
  #preprocessing data
  #function that delete radiance with missing observations
  nacols <- na_cols(rads)
  if (!is.null(num.threshold)) {
    cutoff <- num.threshold
  } else if (!is.null(per.threshold)) {
    cutoff <- per.threshold*nrow(rads)
  }
  rads <- rads[,which(nacols >= cutoff), drop = FALSE]
  wavelength <- wavelength[,which(nacols >= cutoff), drop = FALSE]
  waves <- (1:1016)[which(nacols >= cutoff)]
  narows <- na_rows(rads)
  if (!is.null(narows)) {
    rads <- rads[-narows, ]
    locs <- locs[-narows, ]
    wavelength <- wavelength[-narows,]
  }
  # make sure it is ordered by id(time) and footprint
  rads <- rads[order(locs$id, locs$ftprint),]
  wavelength <- wavelength[order(locs$id, locs$ftprint),]
  locs <- locs[order(locs$id, locs$ftprint),]
  options(scipen = 0)
  return(list(locs = locs, rads = rads/(10^19), waves = waves, wavelength = wavelength))
}

#function that delete radiance with missing observations
na_cols <- function(mat) {
  na_col <- c()
  n <- ncol(mat)
  apply(mat, 2, function(x) {sum(!is.na(x))})
}
#function that delete observations with no radiance
na_rows <- function(mat) {
  na_row <- c()
  m <- nrow(mat)
  na_row <- (1:m)[apply(mat, 1, function(i) {all(is.na(i))})]
  if (length(na_row) == 0) {
    return(NULL)
  } else {
    return(na_row)
  }
}

## function for organizing coastal data (only), and determine the mixing regions, used in mixing problem
# data is the output from Getdata, upper.bd is the cutoff for the range of mean radiance
# return the bds containing information about mixing regions: list of latitude of boundaries in mixing region
Organizedata <- function(data, upper.bd) {
  locs <- data$locs
  rads <- data$rads
  # preprocessing: remove observations with mean beyond upper bound
  out.idx <- rowMeans(rads, na.rm = TRUE) > upper.bd
  locs <- locs[!out.idx, ]
  rads <- rads[!out.idx, ]
  mixes <- which(locs$mix < 100 & locs$mix > 0)
  # identify the mixing regions: upper and lower latitude
  locs1 <- locs[order(locs$lat, decreasing = FALSE),]
  mixes1 <- with(locs1, which(mix < 100 & mix > 0))
  bds <- list()
  s.idx <- min(mixes1) # the first index must be a start of a mixing region
  e.idx <- NULL
  i <- 1 # the first mixing region
  for (m in mixes1[-1]) {
    if (is.null(s.idx)) { # everytime s.idx is set as NULL, it is after a stable region
      s.idx <- m
    } else {
      # see if it is the end index
      e.landfrac <- ifelse(locs1$mix[s.idx-1] == 0, 100, 0) # cover type of the end index
      if ((m+16 < nrow(locs1)) && all(locs1$mix[(m+1):(m+16)] == e.landfrac)) {
        e.idx <- m
      } else if (all(locs1$mix[(m+1):(m+16)] == locs1$mix[s.idx-1])) { # pass the region
        s.idx <- NULL
      }
    }
    if (!is.null(s.idx) & !is.null(e.idx)) {
      bds[[i]] <- c(locs1$lat[s.idx], locs1$lat[e.idx])
      i <- i+1
      # start over to find next one
      s.idx <- NULL 
      e.idx <- NULL 
    }
  }
  return(list(locs = locs, rads = rads, mixes = mixes, bds = bds, reg.num = length(bds)))
}

