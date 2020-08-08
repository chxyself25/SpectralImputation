## data selection and testing, explore several datasets to choose one as an approporaite validation set

dat17 <- readRDS("../Downloads/OceanLocations17.rds") # one track all through 2017
dat177 <- readRDS("../Downloads/OceanLocations177_1.rds") # first part of July of 2017 in the region selected
dat <- readRDS("../2017July/17xiOrbits.rds") # first pc score in one track of 2017
dat <- readRDS("../2017July/177xiOrbits.rds") # first PC score in July of 2017 in selected region


## tracks on pacific next to China and Austrilia: select from multiple tracks of 2017 set
rad17 <- readRDS("../Downloads/OceanRadiance17.rds")
wvs17 <- readRDS("../Downloads/OceanWavelength17.rds")
for (obt in unique(dat17$orbit)) {
  idx <- which(dat17$orbit == obt)
  wvs <- wvs17[idx]
  cat(obt, "\n")
  #print(range(sapply(wvs, length))); 
  print(length(wvs))
}
## 14793, 16191, 17123 have three most locations

## tracks on pacific ocean and next to United States: select from 5 days in 2017
# 15985, 15970, 15942 have three most locations



## explore the dataset, if there can be appropriate points for imputation
sets <- c(14793, 15985, 17123)
for (s in sets) {
  locs <- readRDS(paste0("./pacific_data/nasa_", s, "_locs.rds"))
  locs <- subset(locs, mix == 0 & lat < 70) # only water pixels
  # find center points where 8by8 grids around satisfies conditions: no track missing (if time period is more than 0.5)
  # and no more than 2 points missing each track
  # FPCA selected region around center point have at least 200 points
  tks <- unique(locs$track)
  cpts <- which(locs$ftprint==4)
  valids <- NULL
  for (idx in cpts) {
    pcareg <- subset(locs, lat < locs$lat[idx]+0.25 & lat > locs$lat[idx]-0.25)
    if (nrow(pcareg) < 164) { # 164 = 100+8*8
      #cat("not enough points for doing FPCA", "\n")
      next
    }
    # check no track is missing
    tkidx <- which(tks == locs$track[idx])
    candtk <- tks[(tkidx-4):(tkidx+3)]
    candsec <- sapply(candtk, function(x) {
      as.numeric(substr(x,1,2))*3600+as.numeric(substr(x,3,4))*60+as.numeric(substr(x,5,6))+as.numeric(substr(x,7,7))/10
    })
    if (any(diff(candsec) > 0.5)) {
      cat("tracks are not continuous", "\n")
      next
    }
    locsi <- subset(locs, track %in% candtk)
    if (nrow(locsi) > 64) {
      stop("something is wrong, check code!")
    }
    mincheck <- min(table(locsi$track)) < 6
    maxcheck <- max(table(locsi$track)) < 8
    if (mincheck || maxcheck) {
      cat("too much missing in the 8by8 grid", "\n")
      next
    }
    # if all above passed
    valids <- rbind(valids, data.frame(idx = locs$id[idx], diff = max(diff(candsec)), num = nrow(locsi), 
                                       total = nrow(pcareg), stringsAsFactors = FALSE))
  }
  saveRDS(valids, file = paste0("./pacific_data/nasa_", s, "_validpoints.rds"))
}

