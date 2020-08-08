## look into measurement error in OCO2 data
## it depends on footprint and is a function of wavelength
library(ncdf4)
library(fdapace)
library(ggplot2)
library(reshape2)
library(ggmap)
register_google(key = "AIzaSyDYmdcRD3ItDPhdf16avSm9ITjP9TmWs1w")
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("./SpecFPCA/", trace = FALSE)
source("./Datafuncs.R")

###### global locations map using NASA data#############
# locs_info <- NULL
# for (d in c(90,180,270)) {
#   locsd <- readRDS(paste0("../Downloads/2017_day", d, "_locations.rds"))[, c("id", "orbit", "long", "lat", "date", "ftprint")]
#   locsd <- subset(locsd, date %in% c("20170331", "20170629", "20170927"))
#   locs_info <- rbind(locs_info, locsd)
# }
locs_info <- readRDS("../Downloads/2017_day180_glint_locs.rds")
locs_info <- subset(locs_info, date == "20170629")
sbbox <- make_bbox(lon = locs_info$long, lat = locs_info$lat, f = 0.01)
sq_map <- get_stamenmap(sbbox, maptype = "toner-lite", zoom = 1, color = "bw")
#sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google", zoom = 1)
ggmap(sq_map) + geom_point(data = locs_info, aes(long, lat, color = as.factor(orbit)), size = 1) + theme_bw() +
  theme(text = element_text(size = 13)) + xlab("longitude") + ylab("latitude") + labs(color = "Orbit") +
  guides(colour = guide_legend(override.aes = list(size=2)))
ggsave("../CSSMpresentation/global_locs_map_0629.png", width = 8, height = 4)
## radiance plot example
rads <- readRDS("../Downloads/2017_day180_glint_rads.rds")
waves <- readRDS("../Downloads/2017_day180_glint_waves.rds")
wavelength <- readRDS("../Downloads/2017_day180_glint_wavelength.rds")
idx <- with(locs_info, which(lat < 35.26 & lat > 35.23 & orbit == "15921" & ftprint == 4))[1]
locs_info[idx,]
ggplot(data.frame(waves=waves[[idx]], rads = rads[[idx]]), aes(waves, rads)) + geom_line(color = "black") + 
  xlab("wavelength index") + ylab("radiance") + 
  theme_bw() + theme(text= element_text(size = 16))
ggsave("../CSSMpresentation/radiance_example_15921.pdf", width = 5.2, height = 4)
ggplot(data.frame(wavelength=wavelength[[idx]], rads = rads[[idx]]), aes(wavelength, rads)) + geom_line(color = "black") + 
  xlab("wavelength") + ylab("radiance") + 
  theme_bw() + theme(text= element_text(size = 16))
ggsave("../CSSMpresentation/radiance_example_15921_wavelength.pdf", width = 5.2, height = 4)

## spatial layout map for illustration in section 2
ggplot(subset(locs_info[-idx,], lat < 35.5 & lat > 35 & orbit == "15921"), aes(long, lat)) + 
  geom_point(aes(color = as.factor(ftprint))) + xlab("longitude") + ylab("latitude") +
  scale_color_discrete(name = "footprint") + theme_bw() + theme(text= element_text(size = 16)) +
  guides(colour = guide_legend(override.aes = list(size=1.5))) + 
  geom_point(aes(locs_info[idx,"long"], locs_info[idx,"lat"]), colour = "black")
ggsave("../CSSMpresentation/ftprint_layout_15921.pdf", width = 5.2, height = 4)

############### start using crossing data ################
seg <- Getdata(file = "./coastal_data/crossing_7003545_oco2_radiance.h5")
rads <- seg$rads; waves <- seg$waves
locs <- seg$locs; wavelength <- seg$wavelength
## spatial presentation colored by 1st radiance
library(plotly)
locs$rads <- rads[,1]
fig <- plot_ly(subset(locs, lat > 34.5), x = ~long, y = ~lat, z = ~rads, size = 2, color = I('black'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'longitude', tickfont = list(size = 15), titlefont = list(size = 24)),
                      yaxis = list(title = 'latitude', tickfont = list(size = 15), titlefont = list(size = 24)),
                      zaxis = list(title = 'radiance at index 91', tickfont = list(size = 15), titlefont = list(size = 24))))
fig
## function for calculating measurement error
Getmerror <- function(rads) {
  res <-colMeans(diff(rads, differences = 2)^2, na.rm = TRUE)/6
  return(res)
}

## measurement error estimation for using only water pixels
w.idx <- with(locs, which(mix == 0 & lat >= 34.3 & lat <= 34.8))
# by footprint
e.res <- NULL
for (ft in 1:8) {
  ft.idx <- intersect(w.idx, which(locs$ftprint == ft))
  e.ft <- Getmerror(rads[ft.idx,])
  mu.ft <- colMeans(rads[ft.idx,])
  e.res <- rbind(e.res, data.frame(error = e.ft, waves = waves, mu = mu.ft, ftprint = ft, 
                                   wavelength = colMeans(wavelength[ft.idx,], na.rm = TRUE)))
}
# measurement error against wavelength
ggplot(e.res, aes(waves, sqrt(error))) + geom_line(aes(color = as.factor(ftprint))) + xlab("wavelength index") + 
  ylab("standard error") + scale_color_discrete(name = "footprint") + 
  theme_bw() + theme(text = element_text(size = 13))
ggsave("../CSSMpresentation/sd_vs_waves.pdf", width = 5.2, height = 4)
ggplot(e.res, aes(wavelength, sqrt(error))) + geom_line(aes(color = as.factor(ftprint))) + xlab("wavelength") + 
  ylab("standard error") + scale_color_discrete(name = "footprint") + 
  theme_bw() + theme(text = element_text(size = 13))
ggsave("../CSSMpresentation/sd_vs_wavelength.pdf", width = 5.2, height = 4)
# measurement error against mean function
ggplot(e.res, aes(mu, sqrt(error))) + geom_point(aes(color = as.factor(ftprint))) + xlab("mean radiance") + 
  ylab("standard error") + scale_color_discrete(name = "footprint") + theme_bw()
# get the proportion
prop <- c()
for (ft in 1:8) {
  res.ft <- subset(e.res, ftprint == ft)
  lm.ft <- lm(sqrt(error) ~ mu, data = res.ft)
  prop <- c(prop, coef(lm.ft)[2])
}
prop

## spatial layout map for illustration in section 2
ggplot(subset(locs, lat < 35.5 & lat > 35), aes(long, lat)) + geom_point(aes(color = as.factor(ftprint))) + xlab("longitude") + ylab("latitude") +
  scale_color_discrete(name = "footprint") + theme_bw() + theme(text= element_text(size = 13))
ggsave("../CSSMpresentation/ftprint_layout.pdf", width = 5.2, height = 4)
## radiance plot example
ggplot(data.frame(waves=waves, rads = rads[100,]), aes(waves, rads)) + geom_line(color = "black") + 
  xlab("wavelength index") + ylab("radiance") + 
  theme_bw() + theme(text= element_text(size = 13))
ggsave("../CSSMpresentation/radiance_example.pdf", width = 5.2, height = 4)
#plot(waves, rads[100,], type = "l", xlab = "wavelength index", ylab = "radiance")

## mean function in region latitude in (35, 35.5)
# determine randiance range in plot:
rads_range <- c(0.25, 22.8)
idx1 <- which(locs$lat <= 35.5 & locs$lat >= 35)
ft.list <- lapply(1:8, function(i) {which(locs$ftprint == i)})
mean.ft <- sapply(1:8, function(i) {
  colMeans(rads[intersect(idx1, ft.list[[i]]), ], na.rm = TRUE)
})
wave.ft <- sapply(1:8, function(i) {
  colMeans(wavelength[intersect(idx1, ft.list[[i]]), ], na.rm = TRUE)
})
colnames(mean.ft) <- as.character(1:8)
res1 <- melt(as.data.frame(mean.ft), variable.name = "ftprint", value.name = "radiance")
res1$waves <- rep(waves, 8)
res1$wavelength <- c(wave.ft)
ggplot(res1, aes(waves, radiance)) + geom_line(aes(color = ftprint)) + theme_bw() + ylim(rads_range) + 
  xlab("wavelength index") + ylab("mean radiance") + labs(color = "footprint") +
  theme(text= element_text(size = 13))
ggsave("../CSSMpresentation/ftmean_35.pdf", width = 5.2, height = 4)
ggplot(res1, aes(wavelength, radiance)) + geom_line(aes(color = ftprint)) + theme_bw() + ylim(rads_range) + 
  xlab("wavelength") + ylab("mean radiance") + labs(color = "footprint") +
  theme(text= element_text(size = 13))
ggsave("../CSSMpresentation/ftmean_35_wavelength.pdf", width = 5.2, height = 4)

## mean function in region latitude in (34 34.5)
idx2 <- which(locs$lat <= 34.5 & locs$lat >= 34)
mean.ft <- sapply(1:8, function(i) {
  colMeans(rads[intersect(idx2, ft.list[[i]]), ], na.rm = TRUE)
})
wave.ft <- sapply(1:8, function(i) {
  colMeans(wavelength[intersect(idx2, ft.list[[i]]), ], na.rm = TRUE)
})
colnames(mean.ft) <- as.character(1:8)
res2 <- melt(as.data.frame(mean.ft), variable.name = "ftprint", value.name = "radiance")
res2$waves <- rep(waves, 8)
res2$wavelength <- c(wave.ft)
ggplot(res2, aes(waves, radiance)) + geom_line(aes(color = ftprint)) + theme_bw() + ylim(rads_range) + 
  xlab("wavelength index") + ylab("mean radiance") + labs(color = "footprint") +
  theme(text= element_text(size = 13))
ggsave("../CSSMpresentation/ftmean_34.pdf", width = 5.2, height = 4)
ggplot(res2, aes(wavelength, radiance)) + geom_line(aes(color = ftprint)) + theme_bw() + ylim(rads_range) + 
  xlab("wavelength") + ylab("mean radiance") + labs(color = "footprint") +
  theme(text= element_text(size = 13))
ggsave("../CSSMpresentation/ftmean_34_wavelength.pdf", width = 5.2, height = 4)

# radiance is linear with latitude for fixed wavelength
idx <- which(locs$ftprint == 1)
plot(locs$lat[idx], rads[idx,100], pch = 20)

##########################################################################
############ some plots for oral prelim slides ###########################
##########################################################################
## visualization for illustrating large amount of missing locations
locs0 <- readRDS("./pacific_data/nasa_14793_locs.rds")
locs0s <- subset(locs0, lat < 0 & lat > -5)
ggplot(locs0s, aes(long, lat)) + geom_point(aes(color = as.factor(ftprint))) + 
  scale_color_discrete(name = "footprint") + theme_bw()
## radiance missing plots
files <- list.files(path = "./coastal_data/")
all.rads <- list(NULL, NULL, NULL)
for (f in files) {
  nc1 <- nc_open(paste0("./coastal_data/", f))
  meas_rads <- t(ncvar_get(nc1, "measured_radiance"))
  for (i in 1:3) {
    radsi <- meas_rads[,((i-1)*1016+1):(i*1016)]
    narow <- na_rows(radsi)
    if (!is.null(narow)) {
      radsi <- radsi[-narow, ]
    }
    all.rads[[i]] <- rbind(all.rads[[i]], radsi)
  }
}
miss.res <- NULL; band <- c("O2", "WCO2", "SCO2")
for (i in 1:3) {
  na.per <- apply(all.rads[[i]], 2, function(x) {sum(is.na(x))/length(x)})
  miss.res <- rbind(miss.res, data.frame(na.per = na.per, waves = ((i-1)*1016+1):(i*1016), band = band[i]))
}
ggplot(miss.res, aes(waves, na.per)) + geom_point(aes(color = band)) + xlab("wavelength index") + 
  ylab("Percentage of missing radiance") + theme_bw()
## visualization for illustrating land fraction is bad
seg <- Getdata(file = "./coastal_data/crossing_7001333_oco2_radiance.h5")
locs <- seg$locs
seg_map <- readRDS("../Mixing/Coastal_mixing_maps_nopoints.rds")
ggmap(seg_map[[2]]) + geom_point(data = locs, mapping = aes(x = long, y = lat, color = mix), size = 2) + 
  scale_color_gradient2(low = "white", high = "red", mid = "yellow", midpoint = 50, name = "land fraction")
true_frac <- readRDS("./Boundary_Data/mixing_locs_bds_true_all4.rds")
ggmap(seg_map[[2]]) + geom_point(data = subset(true_frac, orbit == 5216), mapping = aes(x = long, y = lat, color = true_frac), size = 2) + 
  scale_color_gradient2(low = "white", high = "red", mid = "yellow", midpoint = 0.5, name = "land fraction")
## principle component scores plot
library(tikzDevice)
seg <- Getdata(file = "./coastal_data/crossing_7001333_oco2_radiance.h5");
waves <- seg$waves
data <- Organizedata(seg, upper.bd = 30)
locs <- data$locs; rads <- data$rads
mix <- data$bds[[2]]+c(-0.005, 0.005)
L1 <- mix[1]; L2 <- mix[2]
area1 <- with(locs, which(lat <= L1 & lat >= (L1 - 0.5)))
area2 <- with(locs, which(lat >= L2 & lat <= (L2 + 0.5)))
rads2 <- lapply(1:length(area2), function(i) {rads[area2[i],]})
waves2 <- rep(list(waves), length(area2))
locs2 <- locs[area2,]
pca <- SpecFPCA(rads2, waves2, optns = list(FVEthreshold = 0.99, ftprint = locs2$ftprint, userMu = NULL))
locs2$xi <- pca$xiEst[,1]
library(latex2exp)
ggplot(locs2, aes(lat, xi)) + geom_point(aes(color = as.factor(ftprint))) + scale_color_discrete(name = "footprint") +
  ylab(TeX("$\\hat{\\xi}_1$")) + theme_bw()
seg <- Getdata(file = "../Coastal/crossing_7001513_oco2_radiance.h5")



