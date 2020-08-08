## function for doing FPCA on radiance across all locations in the selected region
## mean function: average, covariance function: sample covariance, eigen decomposition
## xi estimation: integral over wavelength index, measurement error: average over 8 footprints
## optns: FVEthreshold, ftprint, userMu

SpecFPCAv3 <- function(Ly, Lt, optns = list(FVEthreshold = 0.99, ftprint = NULL, userMu = NULL)) {
  
  # Check the data validity for further analysis
  CheckData(Ly,Lt)
  
  inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <- inputData$Ly;
  Lt <- inputData$Lt;
  
  # Set the options structure members that are still NULL
  FVEthreshold <- optns$FVEthreshold
  ftprint <- optns$ftprint
  userMu <- optns$userMu
  optns$shrink <- FALSE
  
  # Check the options validity for the PCA function. 
  numOfCurves = length(Ly);
  
  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new data
  # nRegGrid = length(obsGrid) is data is dense
  obsGrid = sort(unique( c(unlist(Lt))));
  ymat <- List2Mat(Ly, Lt) # ymat is on obsGrid
  
  ## common Mean function
  # If the user provided a mean function use it
  if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
    mu <- userMu$mu # must be on the obsGrid
  } else { # cross-sectional mean
    mu <- colMeans(ymat, na.rm = TRUE)
  }
  
  ## footprint-based mean function estimation: multiplicator * mu
  ft.count <- as.vector(table(ftprint)[c('1', '2', '3', '4', '5', '6', '7', '8')])
  ft.idx <- lapply(1:8, function(i) {which(ftprint == i)})
  betas <- sapply(1:8, function(i) {
    sum(apply(ymat[ft.idx[[i]],], 1, "*", mu), na.rm = TRUE)/(ft.count[i]*sum(mu^2, na.rm= TRUE))
  })
  mumat <- matrix(rep(mu, numOfCurves), byrow = TRUE, nrow = numOfCurves)
  for (i in 1:8) {
    mumat[ft.idx[[i]],] <- betas[i] * mumat[ft.idx[[i]],]
  }

  ## Covariance function and sigma2
  sigma_estp <- do.call('rbind', lapply(1:8, function(i) {
    colMeans(diff(ymat[ft.idx[[i]], ], differences = 2)^2, na.rm = TRUE)/6
  }))
  sigma_est <- colSums(sweep(sigma_estp, 1, ft.count, "*"), na.rm = TRUE)/(numOfCurves - 1)
  scsObj = GetCovDense(ymat, mu, optns, sigma_est, mumat)
  sigma2 <- scsObj[['sigma2']]
  
  eigObjObs = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, obsGrid, FVEthreshold, muWork = mu)
  phi <- eigObjObs$phi
  lambda <- eigObjObs$lambda
  
  # Get scores  
  scoresObj <- GetINScores(ymat, obsGrid, optns, mu, lambda, phi, sigma2, mumat) # we don't use sigma2 unless shrink is set
  
  # Make the return object
  ret <- list(waves = obsGrid, mu=mu, smoothCov = scsObj$smoothCov, fittedCov = eigObjObs$fittedCov, sigma_est = sigma_est,
              phi=phi, lambda=lambda, cumFVE = eigObjObs$cumFVE, kChoosen = eigObjObs$kChoosen, xiEst = scoresObj$xiEst,
              betas = betas)
  
  return(ret)
}

