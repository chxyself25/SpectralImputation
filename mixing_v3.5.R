## test for land fraction estimation, make it work
library(fdapace)
library(gstat)
library(sp)
library(geosphere)
library(Matrix)
library(ncdf4)
library(locpol)
library(akima)
library(dplyr)
library(reshape2)
library(ggplot2)
sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}
sourceDir("./SpecFPCA/", trace = FALSE)
source("./model_v3.5/Imputefuncs_v3.5.R")
source("./model_v3.5//Mixingfuncs_v3.5.R")
source("./Datafuncs.R")

## test on Coastal data
true_frac <- readRDS("./Boundary_Data/mixing_locs_bds_true_all4.rds")
# 5216: second region: use 2nd pc for second area, hcv = c(0.01,0.1), c(0.02,0.02)
seg <- Getdata(file = "./coastal_data/crossing_7001333_oco2_radiance.h5");
# 5449: second region: use 2nd pc for second area, hcv = c(0.1,0.1), c(0.02,0.1)
#seg <- Getdata(file = "./coastal_data/crossing_7001513_oco2_radiance.h5")
waves <- seg$waves
data <- Organizedata(seg, upper.bd = 30)
locs <- data$locs; rads <- data$rads
res <- NULL
for (i in 1:data$reg.num) {
  resi <- MixingEst(rads, waves, locs, mix = data$bds[[i]]+c(-0.005, 0.005), threshold = 90, 
                    optns = list(FVEthreshold = 0.99, userMu = NULL, sp.dp = NULL, 
                                 cutoff = 20, width = 2, ft.sp = FALSE),
                    hxi = c("smooth", "smooth"))
  if (!is.null(resi)) {
    resi$region <- i 
  }
  res <- rbind(res, resi)
}
saveRDS(res, file = "./Results/crossing_5449_landfraction_estimate.rds")

###### visualize results ###############
orbit <- 5449
res <- readRDS(paste0("./Results/crossing_", orbit, "_landfraction_estimate.rds"))
res2 <- left_join(res, true_frac[, c("long", "orbit", "id", "true_frac")], by = "id")
res2 %>% group_by(region) %>% summarize(ssehat = mean((ahat-true_frac)^2),
                                        sse = mean((alpha-true_frac)^2)) %>% as.data.frame
forplot <- melt(res2[,c("id", "lat", "region", "alpha", "ahat", "true_frac")], id.vars = c("id", "lat", "region"))
for (reg in 1:2) {
  ggplot(subset(forplot, region == reg), aes(lat, value, color = variable)) + geom_point(size = 1.5) + geom_line() + 
    theme_bw() + xlab("latitude") + ylab("Land Fraction") +
    theme(text = element_text(size=15), legend.position = "none", plot.margin = margin(5, 18, 5, 5)) +  
    scale_color_manual(breaks=c("alpha", "ahat", "true_frac"), 
                       values = c("red", "black", "blue"),
                       labels = c("from OCO-2 dataset", "by unmixing estimation", "by manual calculation"), 
                       name = "Land Fraction") 
  ggsave(paste0("./Results/", orbit, "_mix", reg, "_nl.pdf"), device = 'pdf', width = 5.5, height = 4)
}
##### compare land fraction in map ###########
seg_map <- readRDS("../Mixing/Coastal_mixing_maps_nopoints.rds")
locs$ahat <- locs$mix/100
locs$ahat[match(res2$id,locs$id)] <- res2$ahat
library(ggmap)
a.map <- ggmap(seg_map[[6]]) + geom_point(data = locs, mapping = aes(x = long, y = lat, color = mix/100), size = 2) + 
  scale_color_gradient2(low = "white", high = "red", mid = "yellow", midpoint = 0.5, name = "land fraction")
ahat.map <- ggmap(seg_map[[6]]) + geom_point(data = locs, mapping = aes(x = long, y = lat, color = ahat), size = 2) + 
  scale_color_gradient2(low = "white", high = "red", mid = "yellow", midpoint = 0.5, name = "land fraction")
ggsave("~/OralPrelim/5449upper_compare.pdf", grid.arrange(a.map, ahat.map, nrow = 1), width = 9, height = 4)

# visualize dataset
library(ggmap)
register_google(key = "AIzaSyBCnjdc7irPc2bZUUA8o8AuN14G6_bqs00")
for (orbit in c(5216, 5449)) {
  if (orbit == 5216) {
    seg <- Getdata(file = "./coastal_data/crossing_7001333_oco2_radiance.h5");
  }
  if (orbit == 5449) {
    seg <- Getdata(file = "../Coastal/crossing_7001513_oco2_radiance.h5")
  }
  data <- Organizedata(seg, upper.bd = 30)
  locs <- data$locs; print(unique(locs$orbit) == orbit)
  sbbox <- make_bbox(lon = locs$long, lat = locs$lat, f = 0.1)
  sbbox[1] <- sbbox[1] - 0.3; sbbox[3] <- sbbox[3] + 0.3
  sq_map <- get_map(location = sbbox, maptype = "satellite", source = "google", zoom = 8)
  ggmap(sq_map) + geom_point(data = locs, mapping = aes(x = long, y = lat, color = mix)) + 
    scale_color_gradient2(low = "blue", high = "red", mid = "yellow", midpoint = 50, name = "Land fraction") + 
    geom_hline(yintercept = c(data$bds[[1]]+c(-0.005, 0.005), data$bds[[2]]+c(-0.005, 0.005))) + 
    ylab("latitude") + xlab("longitude") + theme(text = element_text(size = 13))
  ggsave(paste0("../CSSMpresentation/", orbit, "_all_locs_map.pdf"), width = 5, height = 4)
}

