## function for doing FPCA on radiance across all locations in the selected region
## mean function: given location-dependent mean functions, covariance function: sample covariance, eigen decomposition
## xi estimation: substract common mean and integral over wavelength index, measurement error: average over 8 footprints
## optns: FVEthreshold, Mu, Mumat, Sigma_est
## the order should be by footprint

SpecFPCAv2 <- function(Ly, Lt, optns = list(FVEthreshold = 0.99, Mu = NULL, Mumat = NULL, Sigma_est = NULL)) {
  
  # Check the data validity for further analysis
  CheckData(Ly,Lt)
  
  inputData <- HandleNumericsAndNAN(Ly,Lt);
  Ly <- inputData$Ly;
  Lt <- inputData$Lt;
  
  # Set the options structure members that are still NULL
  FVEthreshold <- optns$FVEthreshold
  mu <- optns$Mu
  mumat <- optns$Mumat
  sigma_est <- optns$sigma_est
  optns$shrink <- FALSE
  optns$userMu <- NULL
  
  # Check the options validity for the PCA function. 
  numOfCurves = length(Ly);
  
  # Generate basic grids:
  # obsGrid:  the unique sorted pooled time points of the sample and the new data
  # nRegGrid = length(obsGrid) is data is dense
  obsGrid = sort(unique( c(unlist(Lt))));
  ymat <- List2Mat(Ly, Lt) # ymat is on obsGrid
  
  ## Covariance function and sigma2
  scsObj <- GetCovDense(ymat, mu, optns, sigma_est, mumat)
  sigma2 <- scsObj[['sigma2']]
  
  ## eigen decomposition
  eigObjObs = GetEigenAnalysisResults(smoothCov = scsObj$smoothCov, obsGrid, FVEthreshold, muWork = mu)
  phi <- eigObjObs$phi
  lambda <- eigObjObs$lambda
  
  # Get scores  
  scoresObj <- GetINScores(ymat, obsGrid, optns, mu, lambda, phi, sigma2) # we don't use sigma2 unless shrink is set
  
  # fitted error of FPCA
  # err = sum(apply(scoresObj$fittedY - ymat, 1, function(x) {trapzRcpp(X = obsGrid[!is.na(x)], Y = (x[!is.na(x)])^2)}))
  err <- sum(apply((scoresObj$fittedY - ymat)/ymat, 1, function(x) {mean(x^2, na.rm = TRUE)}))
  
  # Make the return object
  ret <- list(waves = obsGrid, mu = mu, smoothCov = scsObj$smoothCov, #fittedCov = eigObjObs$fittedCov,
              phi=phi, lambda=lambda, cumFVE = eigObjObs$cumFVE, kChoosen = eigObjObs$kChoosen,
              xiEst = scoresObj$xiEst, fittedY = scoresObj$fittedY, err = err)
  
  return(ret)
}

