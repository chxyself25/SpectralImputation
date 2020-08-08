## preparation for simulation on land fraction estimation: function and fpca results
## verify that kriging involved method is better than the simplest interpolation using radiance from nearest footprints
library(fdapace)
library(gstat)
library(sp)
library(geosphere)
library(MASS)
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("./SpecFPCA/", trace = FALSE)
source("./Datafuncs.R")
# source("./model_v3.5/Imputefuncs_v3.5.R")

# ##################### data preparation########################
# ## 5216 orbit as the reference
# seg <- Getdata(file = "./coastal_data/crossing_7001333_oco2_radiance.h5")
# waves <- seg$waves
# locs <- seg$locs; rads <- seg$rads
# ## locations for simulation: latitude and longitude
# lat.diff <- mean(diff(locs$lat))
# long.diff <- mean(diff(locs$long))
# lat <- seq(35, 36, by = lat.diff)
# long <- seq(24, 23, by = long.diff)
# n <- length(lat)
# # the middle pixel is the mixing location
# locs <- data.frame(lat = lat, long = long[1:n],
#                    mix = c(rep(0, floor((n-1)/2)), NA, rep(1, ceiling((n-1)/2))))
# mix.idx <- which(is.na(locs$mix))
# saveRDS(locs, file = "./sim_data/simulation_locs.rds")
# 
# ## mean function and eigen functions for water and land area
# data <- Organizedata(seg, upper.bd = 30)
# locs <- data$locs; rads <- data$rads
# mix <- data$bds[[1]] + c(-0.005, 0.005)
# L1 <- mix[1]; L2 <- mix[2]
# area1 <- with(locs, which(lat <= L1 & lat >= (L1 - 0.5))) # water area
# area2 <- with(locs, which(lat >= L2 & lat <= (L2 + 0.5))) # land area
# aream <- with(locs, which(lat < L2 & lat > L1))
# area1 <- intersect(area1, with(locs, which(mix == 0)))
# area2 <- intersect(area2, with(locs, which(mix == 100)))
# # area1: water area
# rads1 <- lapply(1:length(area1), function(i) {rads[area1[i],]})
# waves1 <- rep(list(waves), length(area1))
# org.list <- OrganizeXi(locs[area1,], locs[aream,])
# pca1 <- SpecFPCAv4(rads1, waves1, optns = list(FVEthreshold = 0.99, ftprint = locs$ftprint[area1],
#                                               X = org.list$X, long = locs$long[area1], lat = locs$lat[area1],
#                                               userMu = NULL, na.action = 'smooth'))
# # area2: land area
# rads2 <- lapply(1:length(area2), function(i) {rads[area2[i],]})
# waves2 <- rep(list(waves), length(area2))
# org.list <- OrganizeXi(locs[area2,], locs[aream,])
# pca2 <- SpecFPCAv4(rads2, waves2, optns = list(FVEthreshold = 0.99, ftprint = locs$ftprint[area2],
#                                              X = org.list$X, long = locs$long[area2], lat = locs$lat[area2],
#                                              userMu = NULL, na.action = 'smooth'))
# saveRDS(list('water' = pca1, 'land' = pca2), file = "./sim_data/water_land_fpca.rds")
# 
# ## variogram model for principle component scores
# # area1: water, with spatial dependence
# coordinates(locs) <- ~long+lat
# proj4string(locs) <- CRS("+proj=longlat +datum=WGS84")
# locs1 <- locs[area1,]; locs1$xi <- pca1$xiEst[,1]
# variog1 <- variogram(xi~1, data = locs1, cutoff = 8, width = 0.4)
# vgm1 <- fit.variogram(variog1, model = vgm("Exp"), fit.method = 1)
# ## range = 9, psill = 480
# # area2: land, with spatial dependence
# locs2 <- locs[area2,]; locs2$xi <- pca2$xiEst[,1]
# variog2 <- variogram(xi~1, data = locs2, cutoff = 20, width = 1)
# vgm2 <- fit.variogram(variog2, model = vgm("Exp"))
# ## range = 5, psill = 940
# 
# ## variance for the second and third princple components
# apply(pca1$xiEst, 2, sd) # 16.75
# apply(pca2$xiEst, 2, sd) # 21.33, 5.05


#############################functions preparation##########################
## function for simulating water, land radiance, and then mixing
# locs: first water then land locations, mixing footprint is marked as NA
# vgm: the list with water and land variogram model
# mu: mean function at waves: the first is water and second is land
# phi: list of eigen functions at waves: first one is water and second is land
# waves: wavelength index
# xi.sd: standard deviation of principle component without spatial dependence
MixingSim <- function(locs, mu, phi, vgm, xi.sd, rsd.e) {
  n <- nrow(locs)
  mu.w <- cbind(1, locs$lat) %*% mu[[1]]; mu.l <- cbind(1, locs$lat) %*% mu[[2]]
  phi.w <- phi[[1]]; phi.l <- phi[[2]]
  vgm.w <- vgm[[1]]; vgm.l <- vgm[[2]]
  xi.sd.w <- xi.sd[[1]]; xi.sd.l <- xi.sd[[2]]
  # generate principle component score based on variogram model
  dist.mat <- distm(locs[, c("long", "lat")], fun = distGeo)/1000
  cov.mat.w <- (vgm.w$psill)*exp(-dist.mat/vgm.w$range)
  xi1.w <- mvrnorm(n = 1, mu = rep(0, n), Sigma = cov.mat.w)
  xi.w <- cbind(xi1.w, mvrnorm(n = n, mu = rep(0, length(xi.sd.w)), Sigma = diag(xi.sd.w^2, length(xi.sd.w))))
  cov.mat.l <- (vgm.l$psill)*exp(-dist.mat/vgm.l$range)
  xi1.l <- mvrnorm(n = 1, mu = rep(0, n), Sigma = cov.mat.l)
  xi.l <- cbind(xi1.l, mvrnorm(n = n, mu = rep(0, length(xi.sd.l)), Sigma = diag(xi.sd.l^2, length(xi.sd.l))))
  # water and land radiance
  #rads.w <- matrix(rep(mu.w, n), byrow = TRUE, nrow = n) + xi.w %*% t(phi.w)
  rads.w <- mu.w + xi.w %*% t(phi.w)
  #rads.l <- matrix(rep(mu.l, n), byrow = TRUE, nrow = n) + xi.l %*% t(phi.l)
  rads.l <- mu.l + xi.l %*% t(phi.l)
  # combine these two processes, mixing proportion is uniform (0,1)
  alpha <- runif(n = 1, min = 0, max = 1)
  mix.idx <- which(is.na(locs$mix))
  mix.rad <- alpha*rads.l[mix.idx,] + (1-alpha)*rads.w[mix.idx,]
  rads <- rbind(rads.w[1:(mix.idx-1), ], mix.rad, rads.l[(mix.idx+1):n, ])
  # standard deviation is proportional to the mean functions
  sde.w <- colMeans(mu.w) * rsd.e
  sde.l <- colMeans(mu.l) * rsd.e
  sde.m <- 0.5*(colMeans(mu.w) + colMeans(mu.l)) * rsd.e
  e.w <- mvrnorm(n = mix.idx-1, mu = rep(0, length(sde.w)), Sigma = diag(sde.w^2))
  e.l <- mvrnorm(n = n-mix.idx, mu = rep(0, length(sde.l)), Sigma = diag(sde.l^2))
  e.m <- mvrnorm(n = 1, mu = rep(0, length(sde.m)), Sigma = diag(sde.m^2))
  res <- rads + rbind(e.w, e.m, e.l)
  return(list(rads = res, alpha = alpha))
}

# test
# vgm1 <- list(psill = 480, range = 9); vgm2 <- list(psill = 940, range = 5)
# xc <- MixingSim(locs, mu = list(pca1$mu, pca2$mu), phi = list(pca1$phi, pca2$phi), vgm = list(vgm1, vgm2), 
#           xi.sd = list(c(16.75), c(21.33, 5.05)), rsd.e = 0.01)


# function that uses nearest two footprints to estimate alpha
# first water and then land footprints
# locs is data frame with coordinates, rads is radiance simulated
IntpEst <- function(locs, rads) {
  mix.idx <- which(is.na(locs$mix))
  radw <- rads[(mix.idx-1), ]
  radl <- rads[(mix.idx+1), ]
  radm <- rads[mix.idx,]
  ahat <- sum((radm-radw)*(radl-radw), na.rm = TRUE)/sum((radw-radl)^2, na.rm = TRUE)
  return(ahat)
}

# function that uses water and land kriging to estimate alpha
# imputatio without smoothing principle component scores
KrigeEst <- function(locs, rads, waves, cutoff = 20, width = 1) {
  n <- nrow(locs)
  coordinates(locs) <- ~long+lat
  proj4string(locs) <- CRS("+proj=longlat +datum=WGS84")
  mix.idx <- which(is.na(locs$mix)); new.loc <- locs[mix.idx,]
  areaw <- 1:(mix.idx-1); locsw <- locs[areaw,]
  areal <- (mix.idx+1):n; locsl <- locs[areal,]
  # imputation on water ftprints
  radsw <- lapply(1:length(areaw), function(i) {rads[areaw[i],]})
  wavesw <- rep(list(waves), length(areaw))
  pcaw <- SpecFPCAv40(radsw, wavesw, optns = list(FVEthreshold = 0.99, userMu = NULL, X = cbind(1, locsw$lat), na.action = 'remove'))
  locsw$xi <- pcaw$xiEst[,1]
  variog <- variogram(xi~1, data = locsw, cutoff = cutoff, width = width)
  vgm <- tryCatch(fit.variogram(variog, model = vgm("Exp")), error = function(c) {NA}, warning = function(c) {
    tryCatch(fit.variogram(variog, model = vgm("Exp"), fit.method = 1), error = function(c) {NA}, warning = function(c) {NA})
  })
  if (is.na(vgm) || vgm$range[2] < 0) {
    vgm <- tryCatch(fit.variogram(variog, 
                                  model = vgm(psill = max(variog$gamma), "Exp", range = cutoff/2, nugget = min(variog$gamma)),
                                  fit.method = 1), error = function(c) {NA}, warning = function(c) {NA})
    if (is.na(vgm) || vgm$range[2] < 0) {
      xi1 <- mean(locsw$xi)
    } else {
      gw <- gstat(id = "xi", formula = xi~1, data = locsw, model = vgm)
      xi1 <- predict(gw, newdata = new.loc)$xi.pred
    }
  } else {
    gw <- gstat(id = "xi", formula = xi~1, data = locsw, model = vgm)
    xi1 <- predict(gw, newdata = new.loc)$xi.pred
  }
  pred.w <- t(c(1, new.loc$lat))%*%pcaw$beta.mat + c(xi1, colMeans(pcaw$xiEst[, -1, drop = FALSE])) %*% t(pcaw$phi)
  # imputation on land footprints
  radsl <- lapply(1:length(areal), function(i) {rads[areal[i],]})
  wavesl <- rep(list(waves), length(areal))
  pcal <- SpecFPCAv40(radsl, wavesl, optns = list(FVEthreshold = 0.99, userMu = NULL, X = cbind(1, locsl$lat), na.action = 'remove'))
  locsl$xi <- pcal$xiEst[,1]
  variog <- variogram(xi~1, data = locsl, cutoff = cutoff, width = width)
  vgm <- tryCatch(fit.variogram(variog, model = vgm("Exp")), error = function(c) {NA}, warning = function(c) {
    tryCatch(fit.variogram(variog, model = vgm("Exp"), fit.method = 1), error = function(c) {NA}, warning = function(c) {NA})
  })
  if (is.na(vgm) || vgm$range[2] < 0) {
    vgm <- tryCatch(fit.variogram(variog, 
                                  model = vgm(psill = max(variog$gamma), "Exp", range = cutoff/2, nugget = min(variog$gamma)),
                                  fit.method = 1), error = function(c) {NA}, warning = function(c) {NA})
    if (is.na(vgm) || vgm$range[2] < 0) {
      xi1 <- mean(locsl$xi)
    } else {
      gl <- gstat(id = "xi", formula = xi~1, data = locsl, model = vgm)
      xi1 <- predict(gl, newdata = new.loc)$xi.pred
    }
  } else {
    gl <- gstat(id = "xi", formula = xi~1, data = locsl, model = vgm)
    xi1 <- predict(gl, newdata = new.loc)$xi.pred
  }
  pred.l <- t(c(1, new.loc$lat))%*%pcal$beta.mat + c(xi1, colMeans(pcal$xiEst[, -1, drop = FALSE])) %*% t(pcal$phi)
  # estimate land fraction
  prad.m <- rads[mix.idx,]
  ahat <- sum((prad.m-pred.w)*(pred.l-pred.w), na.rm = TRUE)/sum((pred.w-pred.l)^2, na.rm = TRUE)
  return(ahat)
}

