# evaluate time-dependent effect of SLP pc1 variability on EBS sea ice
library(tidyverse)
library(sde)
library(rstan)
library(ggpubr)
library(FactoMineR)

# plot settings
theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# begin with North Pacific SLP EOF1

# load slp and cell weights
slp <- read.csv("./data/north.pacific.slp.anom.csv", row.names = 1)

weights <- read.csv("./data/north.pacific.slp.weights.csv", row.names = 1)

pca <- svd.triplet(cov(slp), col.w=weights[,1]) #weighting the columns

pc1 <- as.matrix(slp) %*% pca$U[,1]

# and scale!
pc1 <- as.vector(scale(pc1))

# get winter (NDJFM) means

year <- as.numeric(as.character(years(dates(rownames(slp)))))

# check
change <- year > 2013
year[change] <- year[change]-100
month <- months(dates(rownames(slp)))

# set up winter years
winter.year <- if_else(month %in% c("Nov", "Dec"), year + 1, year)
pc1 <- pc1[month %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
winter.year <- winter.year[month %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]

winter.pc1 <- tapply(pc1, winter.year, mean)

winter.pc1 <- data.frame(year = 1969:2008,
                         slp.pc1 = winter.pc1[names(winter.pc1) %in% 1969:2008])



### load and process ice cover data ------------------------

# note that there are separate time series for 1950-1978 and 1979-present

nc1 <- nc_open("./Data/ERA5_ice_1950-1978.nc")

# process

ncvar_get(nc1, "time")   # hours since 1-1-1900
raw <- ncvar_get(nc1, "time")
h <- raw/24
d1 <- dates(h, origin = c(1,1,1900))
m1 <- months(d1)
yr1 <- years(d1)

x1 <- ncvar_get(nc1, "longitude")
y1 <- ncvar_get(nc1, "latitude")

ice1 <- ncvar_get(nc1, "siconc", verbose = F)
dim(ice1) # 87 long, 37 lat, 203 months

# reverse lat for plotting
ice1 <- ice1[,37:1,]

# reverse y too
y1 <- rev(y1)

ice1 <- aperm(ice1, 3:1)

ice1 <- matrix(ice1, nrow=dim(ice1)[1], ncol=prod(dim(ice1)[2:3]))

# plot to check

ice.mean <- colMeans(ice1, na.rm=T)
z <- t(matrix(ice.mean,length(y1)))
image.plot(x1,y1,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x1, y1, z, add=T)
map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")  

# now the second time series

nc2 <- nc_open("./Data/ERA5_ice_1979-2022.nc")

# process

ncvar_get(nc2, "time")   # hours since 1-1-1900
raw <- ncvar_get(nc2, "time")
h <- raw/24
d2 <- dates(h, origin = c(1,1,1900))
m2 <- months(d2)
yr2 <- years(d2)

x2 <- ncvar_get(nc2, "longitude")
y2 <- ncvar_get(nc2, "latitude")

# expver <-  ncvar_get(nc2, "expver", verbose = F)
# expver # 1 and 5??


ice2 <- ncvar_get(nc2, "siconc", verbose = F)
dim(ice2) # 87 long, 37 lat, 2 expver, 203 months

# expver1 - this is ERA5

# ice2 <- ice2[,,1,]

# reverse lat for plotting
ice2 <- ice2[,37:1,]

# reverse y too
y2 <- rev(y2)

ice2 <- aperm(ice2, 3:1)

ice2 <- matrix(ice2, nrow=dim(ice2)[1], ncol=prod(dim(ice2)[2:3]))


# plot to check

ice.mean <- colMeans(ice2, na.rm=T)
z <- t(matrix(ice.mean,length(y2)))
image.plot(x2,y2,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x2, y2, z, add=T)
map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3")  


# check dimensions
identical(x1, x2)
identical(y1,y2)


# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y1, length(x1))
lon <- rep(x1, each = length(y1))

ice <- rbind(ice1, ice2)

# drop W of 175 and S of 56
drop <- lon < -175 | lat < 56
ice[,drop] <- NA

# plot to check
ice.mean <- colMeans(ice, na.rm=T)
z <- t(matrix(ice.mean,length(y1)))
image.plot(x1,y1,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x1, y1, z, add=T)
map('world2Hires', c('usa', 'USSR'),  fill=T,add=T, lwd=1, col="lightyellow3") # perfecto

dimnames(ice) <- list(as.character(c(d1, d2)), paste("N", lat, "E", lon, sep=""))

f <- function(x) colMeans(x, na.rm = T)

m <- c(as.character(m1), as.character(m2))

yr <- c(as.numeric(as.character(yr1)), as.numeric(as.character(yr2)))

means <- data.frame(month = m,
                    year = as.numeric(as.character(yr)),
                    ice = rowMeans(ice, na.rm = T)) 


ggplot(means, aes(year, ice, color = month)) +
  geom_line()

# drop Oct - Dec
means <- means %>%
  filter(!month %in% c("Oct", "Nov", "Dec"))


# pivot wider
means <- means %>% 
  pivot_wider(values_from = ice, names_from = month) %>%
  filter(year %in% 1953:2022) 

means[,2:5] <- apply(means[,2:5], 2, scale)

plot <- means %>%
  pivot_longer(cols = -year)

ggplot(plot, aes(year, value, color = name)) +
  geom_line()

# generate Jan-Apr means
means$JanApr_ice <- apply(means[,2:5], 1, mean)

# clean up
means <- means %>%
  dplyr::select(year, JanApr_ice)


ggplot(means, aes(year, JanApr_ice)) +
  geom_line() +
  geom_point()

# combine with slp pc1
dat <- left_join(means, winter.pc1) %>%
  dplyr::mutate(era = if_else(year %in% 1969:1988, "1969-1988", "1989-2008")) %>%
  filter(year %in% 1969:2008)

ggplot(dat, aes(slp.pc1, JanApr_ice, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)
