## calculate covariance function for measurement error process
## this is for one ftprint
GetCovEps <- function(ymat) {
  w.num <- ncol(ymat)
  res <- matrix(0, ncol = w.num, nrow = w.num)
  deltaw <- apply(ymat, 2, diff, differences = 2)
  for (w1 in 1:(w.num-1)) {
    for (w2 in (w1+1):w.num) {
      res[w1,w2] <- mean(deltaw[,w1]*deltaw[,w2], na.rm = TRUE)/6
    }
  }
  res <- res+t(res)
  diag(res) <- colMeans(deltaw^2, na.rm = TRUE)/6
  return(res)
}

## calculate variance of e_ik
## this is for one footprint
GetTauEps <- function(sigma_estp, phi, grid) {
  K <- ncol(phi); W <- nrow(phi)
  taus <- lapply(1:K, function(k) {
    temp <- sapply(1:W, function(w) {trapzRcpp(grid, sigma_estp[,w]*phi[,k])})
    trapzRcpp(grid, temp*phi[,k])
  })
  return(unlist(taus))
}