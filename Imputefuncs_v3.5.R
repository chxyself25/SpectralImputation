## imputation functions for doing FPCA and predicting xi
## mean function is beta_p * mu, footprint based mean function

## function for testing if the principle component has spatial dependence, a permutation based test
SDtest <- function(xi0, locs, cutoff, width) {
  n <- nrow(locs)
  locs$xi <- xi0
  variog <- variogram(xi~1, data = locs, cutoff = cutoff, width = width)
  Fvalue <- variog$gamma[1]/var(xi0)
  Fs <- sapply(1:200, function(i) {
    permute <- sample(1:n, size = n)
    locs$xi <- xi0[permute]
    variog <- variogram(xi~1, data = locs, cutoff = cutoff, width = width)
    variog$gamma[1]/var(xi0)
  })
  Fcdf <- ecdf(Fs)
  if (Fcdf(Fvalue) < 0.25 || Fcdf(Fvalue) > 0.975) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

## function for testing dependence for all principal components
SDtestAll <- function(xi, locs, pca, cutoff, width) {
  dps <- rep(NULL, pca$kChoosen)
  for (k in 1:pca$kChoosen) {
    if (k > 1 && pca$cumFVE[k]-pca$cumFVE[k-1] < 1) {
      dps[k] <- FALSE
    } else {
      dps[k] <- SDtest(xi[,k], locs, cutoff = cutoff, width = width)
    }
  }
  return(dps)
}

## function that reorder data by footprint, and calculate design matrix
OrganizeXi <- function(locs, new.locs) {
  # group by footprint
  idx.p <- lapply(1:8, function(p) {which(locs$ftprint == p)})
  # design matrix X
  X <- as.matrix(bdiag(lapply(1:8, function(p) {
    locsp <- locs[idx.p[[p]],]
    # cbind(1, locsp$lat - lat.mean[p])
    cbind(1, locsp$lat)
  })))
  # predictor matrix X0
  X0 <- do.call('rbind', lapply(1:nrow(new.locs), function(i) {
    fp <- new.locs$ftprint[i]; lat <- new.locs$lat[i]
    temp <- rep(0, 16)
    temp[(2*(fp-1)+1):(2*fp)] <- c(1, lat)
    temp
  }))
  # reorder the locations into footprint groups
  return(list(order.idx = unlist(idx.p), X = X, X0 = X0))
}

## function that use kriging to predict xi
KrigeXi <- function(xi, taus, locs, cutoff = 20, width = 1, new.locs = NULL, vgm.str = list(), ft.sp = FALSE) {
  # add geospatial information
  locs$xi <- xi
  coordinates(new.locs) <- ~long+lat
  proj4string(new.locs) <- CRS("+proj=longlat +datum=WGS84")
  # calculate variogram cloud
  variogcloud <- as.data.frame(variogram(xi~1, data = locs, cutoff = cutoff, width = width, cloud = TRUE))
  gamma2 <- variogcloud$gamma - 0.5*(taus[locs$ftprint[variogcloud$left]] + taus[locs$ftprint[variogcloud$right]])
  # transform to sample variogram
  variog <- list(np = NULL, dist = NULL, gamma = NULL)
  grid <- seq(0, cutoff, by = width)
  for (b in 2:length(grid)) {
    idxb <- which(variogcloud$dist <= grid[b] & variogcloud$dist > grid[b-1])
    if (length(idxb) != 0) {
      variog$np <- append(variog$np, as.numeric(length(idxb)))
      variog$dist <- append(variog$dist, mean(variogcloud$dist[idxb]))
      variog$gamma <- append(variog$gamma, mean(gamma2[idxb]))
    }
  }
  variog <- data.frame(variog, dir.hor = 0, dir.ver = 0)
  class(variog) <- c("gstatVariogram", "data.frame")
  # bias correction
  vgm <- eval(parse(text = vgm.str[[1]]))
  if (is.na(vgm) || any(vgm$range < 0)) {
    cat("Set initial values for fitting variogram", "\n")
    vgm <- eval(parse(text = vgm.str[[2]]))
  }
  if (!any(is.na(vgm))) {
    variog$gamma <- variog$gamma - vgm$psill[1]
  }
  # fit variogram model
  vgm <- eval(parse(text = vgm.str[[1]]))
  if (is.na(vgm) || any(vgm$range < 0)) {
    cat("Set initial values for fitting variogram", "\n")
    vgm <- eval(parse(text = vgm.str[[2]]))
    if (is.na(vgm) || any(vgm$range < 0)) {
      cat("This principal component may not have strong spatial dependence, use linear model to fit.", "\n")
      lmres <- LmXi(xi, locs, new.locs, ft.sp)
      return(list(xi.pred = lmres, variog = NULL))
    }
  }
  # covariance matrix
  vgm$psill[1] <- 0 # only keep the covariance part
  Sigmak <- variogramLine(vgm, dist_vector = distm(locs, fun = distGeo)/1000, covariance = TRUE) + diag(taus[locs$ftprint], nrow = nrow(locs))
  sigmak <- variogramLine(vgm, dist_vector = distm(new.locs, locs, fun = distGeo)/1000, covariance = TRUE)
  # beta
  X <- rep(1, nrow(locs))
  beta <- solve(t(X)%*%solve(Sigmak)%*%X) %*% t(X) %*% solve(Sigmak) %*% xi
  # prediction
  xi.pred <- rep(beta, nrow(new.locs)) + sigmak %*% solve(Sigmak) %*% (xi - X %*% beta)
  return(list(xi.pred = xi.pred, variog = variog, sigmak = sigmak, Sigmak = Sigmak))
}

LmXi <- function(xi, locs, new.locs = NULL, ft.sp = FALSE) {
  if (ft.sp) {
    xi.pred <- sapply(new.locs$ftprint, function(x) {
      mean(xi[locs$ftprint == x])
    })
    return(xi.pred)
  } else {
    xi.pred <- rep(mean(xi), nrow(new.locs))
    return(xi.pred)
  }
}

## function that imputes radiance for a missing location
## options: FVEthreshold, userMu, cutoff, width, sp.dp, ft.sp
SpecImpute <- function(rads, waves, locs, new.locs, optns = list()) {
  # reorder and organize xi and locations
  org.list <- OrganizeXi(locs, new.locs)
  rads <- rads[org.list$order.idx]; waves <- waves[org.list$order.idx]
  locs <- locs[org.list$order.idx,]; num.locs <- length(rads)
  cat("Doing FPCA over spatial locations input.", "\n")
  optns.pca <- list(FVEthreshold = ifelse(is.null(optns$FVEthreshold), 0.99, optns$FVEthreshold),
                    ftprint = locs$ftprint, X = org.list$X, long = locs$long, lat = locs$lat,
                    userMu = optns$userMu, na.action = 'remove')
  pca <- SpecFPCAv4(rads, waves, optns = optns.pca)
  xi <- pca$xiEst
  phi <- pca$phi; beta.mat <- pca$beta.mat
  cat("Test spatial dependence: ", "\n")
  coordinates(locs) <- ~long+lat
  proj4string(locs) <- CRS("+proj=longlat +datum=WGS84")
  cutoff <- ifelse(is.null(optns$cutoff), 20, optns$cutoff)
  width <- ifelse(is.null(optns$width), 1, optns$width)
  dps <- optns$sp.dp
  if (is.null(dps)) {
    dps <- SDtestAll(xi, locs, pca, cutoff, width)
    # dps <- rep(FALSE, pca$kChoosen)
    # no warning is allowed
    vgm.str1 <- "tryCatch(fit.variogram(variog, model = vgm('Exp')), error = function(c) {NA}, warning = function(c) {
    tryCatch(fit.variogram(variog, model = vgm('Exp'), fit.method = 1), error = function(c) {NA}, warning = function(c) {NA})})"
    vgm.str2 <- "tryCatch(fit.variogram(variog, model = vgm(psill = max(variog$gamma), 'Exp', range = cutoff/2, nugget = min(variog$gamma)),
    fit.method = 1), error = function(c) {NA}, warning = function(c) {NA})"
  } else {
    # warning is allowed, but not used in this paper
    vgm.str1 <- "tryCatch(fit.variogram(variog, model = vgm('Exp')), error = function(c) {NA}, warning = function(c) {
    tryCatch(fit.variogram(variog, model = vgm('Exp'), fit.method = 1), error = function(c) {NA})})"
    vgm.str2 <- "tryCatch(fit.variogram(variog, model = vgm(psill = max(variog$gamma), 'Exp', range = cutoff/2, nugget = min(variog$gamma)),
    fit.method = 1), error = function(c) {NA})"
  }
  cat("Spatial dependence: ", dps, '\n')
  cat("Calculate PC scores: ", "\n")
  xi.pred <- matrix(0, ncol = pca$kChoosen, nrow = nrow(new.locs))
  for (i in 1:pca$kChoosen) {
    if (dps[i]) {
      taui <- sapply(pca$taus, '[', i)
      krigres <- KrigeXi(xi[,i], taui, locs, cutoff, width, new.locs, vgm.str = list(vgm.str1, vgm.str2), ft.sp = optns$ft.sp)
      xi.pred[,i] <- krigres$xi.pred
      dps[i] <- !is.null(krigres$variog)
    } else {
      xi.pred[,i] <- LmXi(xi[,i], locs, new.locs, optns$ft.sp)
    }
  }
  # imputation results
  new.mumat <- org.list$X0 %*% beta.mat
  rad.pred <- new.mumat + xi.pred %*% t(phi)
  return(list(xi.pred = xi.pred, rad.pred = rad.pred, pca = pca, #variog = krigres$variog, 
              dps = dps, X0 = org.list$X0))
}

