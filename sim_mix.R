## simulation on land fraction estimation
## change relative measurement error and compare kriging based method and linear interpolation
library(doParallel)
registerDoParallel(cores = 8)


## prepare for the simulation
source("./sim_mix_prep.R")
locs <- readRDS("./sim_data/simulation_locs.rds")
pca <- readRDS("./sim_data/water_land_fpca.rds")
#emu <- readRDS("./sim_data/error_mean_func.rds")
pcaw <- pca$water; pcal <- pca$land; waves <- pcaw$waves
vgmw <- list(psill = 5, range = 10); vgml <- list(psill = 10, range = 7)

rsd.e <- seq(0.05, 0.15, length.out = 26)
rmse <- list()
for (s in 1:length(rsd.e)) {
  rse <- foreach(i = 1:200, .combine = "rbind") %dopar% {
    sims <- MixingSim(locs, mu = list(pcaw$beta.mat[7:8,], pcal$beta.mat[7:8,]), phi = list(pcaw$phi, pcal$phi), vgm = list(vgmw, vgml),
              xi.sd = list(c(sqrt(2)), c(sqrt(2), 1)), rsd.e = rsd.e[s])
    rads <- sims$rads
    alpha <- sims$alpha
    pred1 <- IntpEst(locs, rads)
    pred2 <- KrigeEst(locs, rads, waves)
    c((pred1-alpha)/alpha, (pred2-alpha)/alpha)
  }
  rmse[[s]] <- rse
}
saveRDS(rmse, file = "./sim_data/sim_5w10l_10w7l_2w2&1l.rds")

res1 <- readRDS("./sim_data/sim_5w10l_10w7l_2w2&1l_eps_1.rds")
res2 <- readRDS("./sim_data/sim_5w10l_10w7l_2w2&1l_eps_2.rds")
rsd.e <- c(res1[[1]], res2[[1]]); rmse <- c(res1[[2]], res2[[2]])
pred1 <- c()
for (s in 1:length(rsd.e)) {
  pred1 <- c(pred1, mean(abs(rmse[[s]][,1]), trim = 0.1))
}
pred2 <- c()
for (s in 1:length(rsd.e)) {
  pred2 <- c(pred2, mean(abs(rmse[[s]][,2]), trim = 0.1))
}
res <- data.frame(r.sd = rep(rsd.e, 2), method = rep(c('Interpolation', 'Unmixing'), each = length(rsd.e)), 
                  r.bias = c(pred1, pred2))
ggplot(res, aes(r.sd, r.bias, color = method)) + geom_point(size = 1) + geom_line() +
  ylab("mean relative absolute error") + xlab("standard deviation of measurement error") + 
  labs(color = "Method") + theme_bw() + theme(text = element_text(size = 13))
ggsave("./sim_data/sim_mix_plot_eps.pdf", width = 6, height = 4)


