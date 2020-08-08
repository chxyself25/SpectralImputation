# This function obtains the sample covariance matrix at observed grid
# for dense regular functional data

######
# Input:
######  
#  ymat: n by p matrix of dense regular functional data
#  mu: p-dim vector, estimated cross-sectional mean
#  optns: options for FPCA function
#  y: list of amplitude information
#  t: list of time information
######
# Output: 
######
#  a SmoothCov object containing: 
#    - p by p matrix of sample cov surface estimation on observed grid
#    - NULL for all other entires
##########################################################################

GetCovDense <- function(ymat, mu, optns, sigma_est, mumat = NULL) {
  n = nrow(ymat)
  m = ncol(ymat)
  if (!is.null(optns$userMu)) {
    ymat = ymat - matrix(rep(times= nrow(ymat), mu), ncol= ncol(ymat), byrow=TRUE)
    K = matrix( rep(0,m^2), m)
    for( i in (1:m)){
      for( j in (1:m)){
        XcNaNindx = which(is.na(ymat[,i]));
        YcNaNindx = which(is.na(ymat[,j]));
        NaNrows = union(XcNaNindx, YcNaNindx);
        # Find inner product of the columns with no NaN values
        indx = setdiff( 1:n, NaNrows)
        K[i,j] =  sum(ymat[indx,i] * ymat[indx,j]) * (1/(n-1-length(NaNrows)));  
      }
    }
  } else if (!is.null(mumat)) {
    # residuals
    ymat <- ymat - mumat
    K = matrix( rep(0,m^2), m)
    for( i in (1:m)){
      for( j in (1:m)){
        prodij <- ymat[,i]*ymat[,j]
        K[i,j] <- sum(prodij, na.rm = TRUE)/(sum(!is.na(prodij))-1)  
      }
    }
  } else {
    K = cov(ymat, use = 'pairwise.complete.obs')
  }
  K = 0.5 * (K + t(K)) # ensure that K is symmetric
  
  # substract sigma2
  ord <- 2
  sigma2 <- mean(diff(t(ymat), differences=ord)^2, na.rm=TRUE) / choose(2 * ord, ord)
  if (is.null(sigma_est)) {
    diag(K) <- diag(K) - colMeans(diff(ymat, differences = 2)^2, na.rm = TRUE)/6
  } else {
    #diag(K) <- diag(K) - sigma_est
    K <- K - sigma_est
  }
  
  return(list('smoothCov' = K, 'sigma2' = sigma2))
}
