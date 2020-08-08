## evaluation of imputation algorithm
## 128 center points: remove one track each time and repeat 8 times, impute for all removed points
## return metrics: id, long, lat, ftprint, gap, mse, prede, miminum latitude difference around
## region selected for FPCA: plus minus 0.25 of center point
library(fdapace)
library(gstat)
library(sp)
library(geosphere)
library(Matrix)
library(ncdf4)
library(locpol)
library(akima)
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("./SpecFPCA/", trace = FALSE)
source("./model_v3.5/Imputefuncs_v3.5.R")
source("./model_v3.5/helpers.R")
#source("./Datafuncs.R")


########################## imputation ###################################
# read in data and center points for imputation
cpts <- readRDS("./pacific_data/nasa_14793_validpoints.rds")
cpts <- subset(cpts, num == 64 & total >= 164)
locs0 <- readRDS("./pacific_data/nasa_14793_locs.rds")
rads0 <- readRDS("./pacific_data/nasa_14793_rads.rds")
waves0 <- readRDS("./pacific_data/nasa_14793_waves.rds")
### visualize dataset
library(ggmap)
register_google(key = "AIzaSyBCnjdc7irPc2bZUUA8o8AuN14G6_bqs00")
locs0s <- subset(locs0, mix == 0 & lat < 70)
all(cpts$id %in% locs0s$id)
sbbox <- make_bbox(lon = locs0s$long, lat = locs0s$lat, f = 0.1)
sbbox[1] <- 100; sbbox[3] <- 179.9
sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google", zoom = 3)
ggmap(sq_map) + geom_point(data = locs0s, mapping = aes(x = long, y = lat)) + 
  theme(text = element_text(size = 13)) + xlab("longitude") + ylab("latitude")
ggsave("../CSSMpresentation/14793_all_map.pdf", width = 5, height = 4)
ggmap(sq_map) + geom_point(data = subset(locs0s, id %in% cpts$idx[seq(1, 255, 2)]), aes(long, lat)) +
  theme(text = element_text(size = 13)) + xlab("longitude") + ylab("latitude")
ggsave("../CSSMpresentation/14793_valid_map.pdf", width = 5, height = 4)
### wavelength data preparation
# file <- "./pacific_data/oco2_L2DiaGL_14793a_170413_B8100r_170910085430.h5"
# nc1 <- nc_open(file)
# wavedata <- t(ncvar_get(nc1, "SpectralParameters/wavelength"))
# ids <- as.character(ncvar_get(nc1, "RetrievalHeader/sounding_id"))
# smdata <- list(wavedata = wavedata, ids = ids)
### wavelength data preparation
smdata <- readRDS("../pacific_data/nasa_14793_wavelength_id.rds")
library(doParallel)
registerDoParallel(cores = 8)
res <- foreach(id = cpts$idx[seq(1, 255, 2)]) %dopar% {
  idlat <- locs0$lat[locs0$id == id]
  # select a region
  s.idx <- with(locs0, lat < idlat+0.25 & lat > idlat-0.25) # 9.25 and 8.75
  locs <- locs0[s.idx,]
  rads <- rads0[s.idx]
  waves <- waves0[s.idx]
  impute <- remove_gaps(locs, rads, waves, c.id=id,
                        optns = list(FVEthreshold = 0.99, userMu = NULL, cutoff = 20, width = 1, sp.dp = NULL, ft.sp = FALSE), 
                        smdata, 8)
  impute$res
}
saveRDS(res, file = "./model_v3.5/nasa_14793_255by2_impute_remove.rds")

######################### summary and visualization #################################
res <- readRDS(file = "./Results/nasa_14793_255by2_impute_remove_r.rds")
## heatmap
library(gridExtra)
evals <- list()
for (i in 1:8) {
  resi <- lapply(res, function(x) {xi <- subset(x, gap == i); xi[order(xi$track, xi$ftprint),]})
  msei <- sapply(resi, function(x) {
    #dat <- subset(x, gap == i)
    sqrt(x$mse)
  })
  emsei <- sapply(resi, function(x) {
    #dat <- subset(x, gap == i)
    sqrt(x$emse)
  })
  mati <- cbind(expand.grid(footprint = 1:8, track = 1:i), rmse = rowMeans(msei, na.rm = TRUE), 
                rpmse = rowMeans(emsei, na.rm = TRUE))
  mati$gap <- i
  evals[[i]] <- mati
}

library(ggplot2)
# root mean square error
mserg <- range(do.call("rbind", evals)$rmse)
ggplot(do.call("rbind", evals[1:4]), aes(footprint, track)) + geom_raster(aes(fill = rmse)) +
  scale_fill_gradientn(colours = c("blue", "red"), limits = mserg, name = "RMSE") + coord_equal() + 
  facet_grid(.~gap) + scale_x_continuous(breaks = 1:8) + ylab("cross-track") + 
  theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank(),text = element_text(size = 15))
ggsave("./Results/remove_NA/14793_rmse_1-4.png", width = 12, height = 3)
ggplot(do.call("rbind", evals[5:8]), aes(footprint, track)) + geom_raster(aes(fill = rmse)) +
  scale_fill_gradientn(colours = c("blue", "red"), limits = mserg, name = "RMSE") + coord_equal() + 
  facet_grid(.~gap) + scale_x_continuous(breaks = 1:8) + ylab("cross-track") + 
  theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank(),text = element_text(size = 15))
ggsave("./Results/remove_NA/14793_rmse_5-8.png", width = 12, height = 3)
# htmaps <- list()
# for (g in 1:8) {
#   htmaps[[g]] <- ggplot(evals[[g]], aes(footprint, track)) + geom_raster(aes(fill = rmse)) + 
#     scale_fill_gradientn(colours = c("blue", "red"), limits = mserg) + coord_equal() + theme_bw()
# }
# root predictive mean square error
rpmserg <- range(do.call("rbind", evals)$rpmse)
ggplot(do.call("rbind", evals[1:4]), aes(footprint, track)) + geom_raster(aes(fill = rpmse)) +
  scale_fill_gradientn(colours = c("blue", "red"), limits = rpmserg, name = "RPMSE") + coord_equal() + 
  facet_grid(.~gap) + scale_x_continuous(breaks = 1:8) + ylab("cross-track") + 
  theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank(),text = element_text(size = 15))
ggsave("./Results/remove_NA/14793_rpmse_1-4.png", width = 12, height = 3)
ggplot(do.call("rbind", evals[5:8]), aes(footprint, track)) + geom_raster(aes(fill = rpmse)) +
  scale_fill_gradientn(colours = c("blue", "red"), limits = rpmserg, name = "RPMSE") + coord_equal() + 
  facet_grid(.~gap) + scale_x_continuous(breaks = 1:8) + ylab("cross-track") + 
  theme_bw() + theme(strip.background = element_blank(), strip.text.x = element_blank(),text = element_text(size = 15))
ggsave("./Results/remove_NA/14793_rpmse_5-8.png", width = 12, height = 3)
# htmaps <- list()
# for (g in 1:8) {
#   htmaps[[g]] <- ggplot(evals[[g]], aes(footprint, track)) + geom_raster(aes(fill = rpmse)) + 
#     scale_fill_gradientn(colours = c("blue", "red"), limits = rpmserg) + coord_equal() + theme_bw()
# }
# ggsave("./model_v3.5/14793_rpmse_1-4.png", grid.arrange(grobs = htmaps[1:4], nrow = 1, ncol = 4), width = 16, height = 4)
# ggsave("./model_v3.5/14793_rpmse_5-8.png", grid.arrange(grobs = htmaps[5:8], nrow = 1, ncol = 4), width = 16, height = 4)

library(dplyr)
## mse against gaps
msegaps <- lapply(res, function(x) {
  ctk <- sort(unique(x$track))[(i%/%2+1)] # i is 8 from above
  subset(x, track == ctk)
}) %>% do.call("rbind",.)
msegaps <- msegaps %>% 
  group_by(gap, ftprint) %>% 
  summarise(rmse = mean(sqrt(mse)),
            rpmse = mean(sqrt(emse))) %>% as.data.frame
ggplot(msegaps, aes(gap, rmse, color = as.factor(ftprint))) + geom_line() + geom_point(size = 1) +
  scale_color_discrete(name="footprint") + xlab("# of cross-tracks removed") + ylab("RMSE") +
  theme_bw() + theme(text = element_text(size = 13))
ggsave("./Results/remove_NA/rmse_gap_ctk_plot.pdf", width = 6, height = 4)
ggplot(msegaps, aes(gap, rpmse, color = as.factor(ftprint))) + geom_line() + geom_point(size = 1) +
  scale_color_discrete(name="footprint") + xlab("# of cross-tracks removed") + ylab("RPMSE") + 
  theme_bw() + theme(text = element_text(size = 13))
ggsave("./Results/remove_NA/rpmse_gap_ctk_plot.pdf", width = 6, height = 4)

