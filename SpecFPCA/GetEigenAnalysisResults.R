# phi: a nRegGrid * no_FVE
# The input smoothCov is possibly truncated.

GetEigenAnalysisResults <- function(smoothCov, obsGrid, FVEthreshold, muWork = NULL) {
  
  gridSize <- obsGrid[2] - obsGrid[1]
  numGrids <- nrow(smoothCov)
  
  eig <- eigen(smoothCov)

  positiveInd <- eig[['values']] >= 0
  if (sum(positiveInd) == 0) {
    stop('All eigenvalues are negative. The covariance estimate is incorrect.')
  }
  d <- eig[['values']][positiveInd]
  eigenV <- eig[['vectors']][, positiveInd, drop=FALSE]

  # thresholding for corresponding FVE option 
  #(not before to avoid not being able to reach the FVEthreshold when pos eigenvalues > maxk)
  # i.e. default FVE 0.9999 outputs all components remained here.
  FVE <- cumsum(d) / sum(d) * 100  # cumulative FVE for all available eigenvalues from fitted cov
  no_opt <- min(which(FVE >= FVEthreshold * 100)) # final number of component chosen based on FVE
  
  # normalization
  if (is.null(muWork)) {
    muWork = 1:dim(eigenV)[1]
  }

  phi <- apply(eigenV, 2, function(x) {
                    x <- x / sqrt(trapzRcpp(obsGrid, x^2)) 
                    if ( 0 <= sum(x*muWork) )
                      return(x)
                    else
                      return(-x)
  })
  lambda <- gridSize * d; # there maybe an issue that obsGrid is not equally spaced

  # selection on lambda and phi
  lambda = lambda[1:no_opt]; phi = phi[,1:no_opt, drop=FALSE]
  fittedCov <- phi %*% diag(x=lambda, nrow = length(lambda)) %*% t(phi)

  return(list(lambda = lambda, phi = phi, cumFVE = FVE, kChoosen=no_opt, fittedCov=fittedCov))
}
