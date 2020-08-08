## functions for estimating land fractions and mixing process


# function for doing radiance imputation and smoothing the PC scores
# analog to SpecImpute, but with one more step of smoothing PC scores
MixingImpute <- function(rads, waves, locs, new.locs, optns = list(), hxi = "smooth") {
  # reorder and organize xi and locations
  org.list <- OrganizeXi(locs, new.locs)
  rads <- rads[org.list$order.idx]; waves <- waves[org.list$order.idx]
  locs <- locs[org.list$order.idx,]; num.locs <- length(rads)
  cat("Doing FPCA over spatial locations input.", "\n")
  optns.pca <- list(FVEthreshold = ifelse(is.null(optns$FVEthreshold), 0.99, optns$FVEthreshold),
                    ftprint = locs$ftprint, X = org.list$X, long = locs$long, lat = locs$lat,
                    userMu = optns$userMu, na.action = 'smooth')
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
  }
  cat("Spatial dependence: ", dps, '\n')
  # warning is allowed
  vgm.str1 <- "tryCatch(fit.variogram(variog, model = vgm('Exp')), error = function(c) {NA}, warning = function(c) {
    tryCatch(fit.variogram(variog, model = vgm('Exp'), fit.method = 1), error = function(c) {NA})})"
  vgm.str2 <- "tryCatch(fit.variogram(variog, model = vgm(psill = max(variog$gamma), 'Exp', range = cutoff/2, nugget = min(variog$gamma)),
    fit.method = 1), error = function(c) {NA})"
  cat("Calculate PC scores: ", "\n")
  xi.pred <- matrix(0, ncol = pca$kChoosen, nrow = nrow(new.locs))
  for (i in 1:pca$kChoosen) {
    if (dps[i]) {
      # remove outliers
      xibox <- boxplot(xi[,i], plot = FALSE)
      out.idx <- (xi[,i] > xibox$stats[5]) # | (xi[,i] < xibox$stats[1]) # outlier indicator
      out.dist <- min(abs(locs$lat[which.max(xi[,i])] - new.locs$lat)) # minimum distance from outlier to mixing region
      xii <- xi[!out.idx,i]; locsi <- locs[!out.idx,]
      #locsi <- locs; xii <- xi[,i]
      #xii[out.idx] <- xibox$stats[5]
      taui <- sapply(pca$taus, '[', i)
      # doing linear local smoother for PC scores
      if (hxi == "smooth") {
        if (sum(out.idx) > 5 & out.dist > 0.005 & out.dist < 0.25) { # the condition is specified to the data
          #hxii <- regCVBwSelC(locsi$lat, xii, deg=1, kernel=EpaK, interval=c(seq(0.005, 1, length.out = 500)))
          hxii <- 0.1
        } else {
          hxii <- regCVBwSelC(locsi$lat, xii, deg=1, kernel=EpaK, interval=c(seq(0.005, 1, length.out = 500)))
          #hxii <- 0.1
        }
        print(hxii)
        loc.fit <- locpol(y~x, data = data.frame(x = locsi$lat, y = xii), deg = 1, bw = hxii, kernel = EpaK, xeval = locsi$lat)
        xii <- fitted(loc.fit)[rank(locsi$lat)]
      }
      krigres <- KrigeXi(xii, taui, locsi, cutoff, width, new.locs, vgm.str = list(vgm.str1, vgm.str2), ft.sp = optns$ft.sp)
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


# function that warps up the procedure
# mix is the latitude boundaries for mixing region to deal with
# threshold: how many points are valid for doing FPCA, otherwise drop off
# optns: option list for doing radiance imputation 
# hxi: bandwidth for doing local linear smoother for principle component scores: vector for the first area and second area
MixingEst <- function(rads, waves, locs, mix, threshold = 100, optns = list(),
                      hxi = c(NULL, NULL)) {
  cat("Recognize area1 and area2, and check: ", "\n")
  L1 <- mix[1]; L2 <- mix[2]
  area1 <- with(locs, which(lat <= L1 & lat >= (L1 - 0.6)))
  area2 <- with(locs, which(lat >= L2 & lat <= (L2 + 0.6)))
  aream <- with(locs, which(lat < L2 & lat > L1))
  # check if this is valid for processing, use only pure pixels for FPCA
  lf1 <- mean(locs$mix[area1], na.rm = TRUE)
  lf2 <- mean(locs$mix[area2], na.rm = TRUE)
  type1 <- ifelse(lf1 > 70, "land", ifelse(lf1 < 30, "water", "not sure"))
  type2 <- ifelse(lf2 > 70, "land", ifelse(lf2 < 30, "water", "not sure"))
  if (type1 == type2 || any(c(type1, type2) == "not sure")) {
    return(NULL)
  }
  if (type1 == "water") {
    area1 <- intersect(area1, with(locs, which(mix == 0)))
    area2 <- intersect(area2, with(locs, which(mix == 100)))
  } else {
    area1 <- intersect(area1, with(locs, which(mix == 100)))
    area2 <- intersect(area2, with(locs, which(mix == 0)))
  }
  if (length(area1) < threshold | length(area2) < threshold) {
    return(NULL)
  }
  cat("Imputation for the mixing region: ", "\n")
  rads1 <- lapply(1:length(area1), function(i) {rads[area1[i],]})
  waves1 <- rep(list(waves), length(area1))
  impute1 <- MixingImpute(rads1, waves1, locs[area1,], locs[aream,], optns = optns, hxi = hxi[1])
  rad.preds1 <- impute1$rad.pred
  rads2 <- lapply(1:length(area2), function(i) {rads[area2[i],]})
  waves2 <- rep(list(waves), length(area2))
  impute2 <- MixingImpute(rads2, waves2, locs[area2,], locs[aream,], optns = optns, hxi = hxi[2])
  rad.preds2 <- impute2$rad.pred
  cat("Calculate land fraction: ", "\n")
  radsm <- rads[aream,]
  a <- locs$mix[aream]/100
  # sigma_est <- (impute1$pca$sigma_est + impute2$pca$sigma_est)*0.5
  # sigmam <- sigma_est[locs$ftprint[aream],]
  if (type1 == "water") {
    num <- rowSums((radsm - rad.preds1)*(rad.preds2 - rad.preds1))
    den <- rowSums((rad.preds1 - rad.preds2)^2)
    ahat <- sapply(num/den, function(x) {ifelse(x > 1, 1, ifelse(x < 0, 0, x))})
    mse <- rowSums((radsm - diag(a) %*% rad.preds2 - diag(1-a) %*% rad.preds1)^2)
    msehat <- rowSums((radsm - diag(ahat) %*% rad.preds2 - diag(1-ahat) %*% rad.preds1)^2)
  } else {
    num <- rowSums((radsm - rad.preds2)*(rad.preds1 - rad.preds2))
    den <- rowSums((rad.preds1 - rad.preds2)^2)
    ahat <- sapply(num/den, function(x) {ifelse(x > 1, 1, ifelse(x < 0, 0, x))})
    mse <- rowSums((radsm - diag(a) %*% rad.preds1 - diag(1-a) %*% rad.preds2)^2)
    msehat <- rowSums((radsm - diag(ahat) %*% rad.preds1 - diag(1-ahat) %*% rad.preds2)^2)
  }
  res <- data.frame(id = locs$id[aream], alpha = a, ahat = ahat, mse = mse, msehat = msehat, lat = locs$lat[aream])
  return(res)
}

