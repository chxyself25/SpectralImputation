## function for doing FPCA on radiance across all locations in the selected region
## mean function: two dimensional, linear in wvelength direction, covariance function: sample covariance, eigen decomposition
## xi estimation: integral over wavelength index, measurement error: average over 8 footprints
## optns: FVEthreshold, ftprint, userMu

SpecFPCAv40 <- function(Ly, Lt, optns = list(FVEthreshold = 0.99, userMu = NULL)) {
  
  # Check the data validity for further analysis
  CheckData(Ly,Lt)
  
  inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <- inputData$Ly;
  Lt <- inputData$Lt;
  
  # Set the options structure members that are still NULL
  FVEthreshold <- optns$FVEthreshold
  userMu <- optns$userMu
  optns$shrink <- FALSE
  X <- optns$X
  
  # Check the options validity for the PCA function. 
  numOfCurves = length(Ly);
  
  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new data
  # nRegGrid = length(obsGrid) is data is dense
  obsGrid = sort(unique( c(unlist(Lt))));
  ymat <- List2Mat(Ly, Lt) # ymat is on obsGrid
  
  # remove wavelength with NA or do interpolation to impute NAs according to options settings
  if (optns$na.action == 'remove') {
    na.col <- which(apply(ymat, 2, function(x) {any(is.na(x))}))
    if (length(na.col) != 0) {
      obsGrid <- obsGrid[-na.col]
      ymat <- ymat[,-na.col] 
    }
  } else if (optns$na.action == 'smooth') {
    ymat <- NAFill(ymat, optns$long, optns$lat, obsGrid)
  }
  
  ## Mean function estimation
  # If the user provided a mean function use it
  if ( is.list(userMu) && (length(userMu$mu) == length(userMu$t))){
    mu <- userMu$mu # must be on the obsGrid
    mumat <- NULL
  } else { # linear regression by fixing wavelength
    beta.mat <- apply(ymat, 2, function(x) {
      solve(t(X) %*% X) %*% t(X) %*% x
    })
    mumat <- X%*%beta.mat
    mu <- colMeans(mumat)
  }
  
  ## Covariance function and sigma2
  sigma_est <- colMeans(diff(ymat, differences = 2)^2, na.rm = TRUE)/6
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
              beta.mat = beta.mat, mumat = mumat)
  
  return(ret)
}

