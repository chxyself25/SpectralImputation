########################## functions #################################
## function for smoothing radiance function
# file is the directory of the data; new.id is the radiance id we want to smooth
# waves is the wavelength index of radiance we want to smooth (usually is known before smoothing)
# rads is radiance before smoothing, new.id all id for smoothing
rads_smooth <- function(wavedata, ids, new.id, rads, waves, type = FALSE) {
  options(scipen = 999)
  #nc1 <- nc_open(file)
  n <- length(rads)
  #type <- grepl("crossing", file)
  if (type) {
    #wavedata <- t(ncvar_get(nc1, "wavelength"))
    waveidx <- waves # index w.r.t wavelength
    #id <- as.character(ncvar_get(nc1, "sounding_id"))
  } else {
    #wavedata <- t(ncvar_get(nc1, "SpectralParameters/wavelength"))
    waveidx <- lapply(waves, function(w) {1:length(w)}) # index w.r.t wavelength
    #id <- as.character(ncvar_get(nc1, "RetrievalHeader/sounding_id"))
  }
  res <- list()
  for (i in 1:n) {
    xx <- as.vector(wavedata[ids == new.id[i], waveidx[[i]]])
    yy <- as.vector(rads[[i]])
    hcv <- regCVBwSelC(xx,yy,deg=1,kernel=gaussK,interval=c(seq(1e-05,1e-07,length.out=1000)))
    #print(hcv)
    r_rot <- locpol(yy~xx, data = data.frame(xx=xx, yy=yy), deg = 1, bw = hcv, kernel = gaussK, xeval = xx)
    f.s <- fitted(r_rot)
    res[[new.id[i]]] <- f.s
  }
  options(scipen = 0)
  return(res)
}

## function for calculating component score by integration
## xx is the covariate vector of where f.s located
INScores <- function(f.s, pca, xx) {
  beta.mat <- pca$beta.mat
  phi <- pca$phi; waves <- pca$waves
  nnaidx <- which(!is.na(f.s))
  inxi <- sapply(1:pca$kChoosen, function(k) {
    trapzRcpp(waves[nnaidx], (f.s[nnaidx] - t(xx)%*%beta.mat[,nnaidx])*phi[nnaidx,k])
  })
  return(inxi)
}
## function for doing removing gaps
# remove one gap each time, and do imputation for the removed region, compare results with smoothing radiance
# return a list of several removals, each element is a dataframe with id, lat/long, ftprint, mse, prede
# c.id is the center point in the 8by8 grid region for imputation
remove_gaps <- function(locs, rads, waves, c.id, optns = list(), smdata = list(), num.gap=8) {
  tks <- unique(locs$track)
  tkidx <- which(tks == locs$track[locs$id == c.id])
  # radiance smoothing
  if (num.gap%%2 == 0) {
    h <- (num.gap-1)%/%2
    candtk <- tks[(tkidx-h-1):(tkidx+h)]
  }else {
    h <- (num.gap-1)/2
    candtk <- tks[(tkidx-h):(tkidx+h)]
  }
  r.all <- which(locs$track %in% candtk)
  r.id <- locs$id[r.all]
  f.sm <- rads_smooth(wavedata = smdata$wavedata, ids = smdata$ids, r.id, rads[r.all], waves[r.all])
  # remove gaps and imputation
  impute.mat <- list(); res <- NULL
  for (i in 1:num.gap) {
    if (i == 1) {
      candtk <- tks[tkidx]
    }else {
      if (i%%2 == 0) {
        j <- (i-1)%/%2
        candtk <- tks[(tkidx-j-1):(tkidx+j)]
      }else {
        j <- (i-1)/2
        candtk <- tks[(tkidx-j):(tkidx+j)]
      }
    }
    r.lines <- which(locs$track %in% candtk)
    r.locs <- locs[-r.lines,]
    r.rads <- rads[-r.lines]
    r.waves <- waves[-r.lines]
    new.locs <- locs[r.lines,]
    latdiff <- sapply(r.lines, function(x) {min(abs(locs$lat[x]-r.locs$lat))})
    pred <- SpecImpute(r.rads, r.waves, r.locs, new.locs, optns = optns)
    rad.pred <- pred$rad.pred
    wgrid <- pred$pca$waves
    f.smi <- f.sm[new.locs$id]
    # map smoothing function to one universal grid
    f.s <- lapply(1:length(r.lines), function(x) {
      f <- rep(NA, length(wgrid))
      wavex <- waves[r.lines][[x]]
      grid <- intersect(wgrid, wavex)
      f[wgrid %in% grid] <- f.smi[[x]][wavex %in% grid]
      f
    })
    # calculate true xi
    xi <- do.call("rbind", lapply(1:length(r.lines), function(x) {
      INScores(f.s[[x]], pred$pca, pred$X0[x,])
    }))
    mse <- lapply(1:length(r.lines), function(i) {
      mean((rad.pred[i,]-f.s[[i]])^2/(f.s[[i]])^2, na.rm = TRUE)
    })
    emse <- apply((xi-pred$xi.pred) %*% t(pred$pca$phi), 1, function(x) {mean(x^2, na.rm = TRUE)})
    dist <- sapply(as.numeric(as.factor(new.locs$track)), function(x) {min(x-0, i+1-x)})
    impute.mat[[i]] <- rad.pred # predicted radiance function for all removed points
    resi <- cbind(new.locs, data.frame(mse = unlist(mse), emse = emse, K = pred$pca$kChoosen,
                                      gap = i, latdiff = latdiff, dist = dist, dps = sum(pred$dps)))
    res <- rbind(res, resi)
  }
  return(list(res = res, impute.mat = impute.mat))
}


