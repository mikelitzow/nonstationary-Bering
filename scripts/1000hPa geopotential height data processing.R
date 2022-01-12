library(tidyverse)
library(ncdf4)
library(maps)
library(maptools)
library(mapdata)
library(fields)
library(oce)
library(chron)
library(zoo)
library(mapproj)


# paste URL into browser to download
URL <- "https://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_941a_c863_cca9.nc?hgt[(1948-01-01):1:(2021-12-01T00:00:00Z)][(1000):1:(1000.0)][(20):1:(90.0)][(0.0):1:(357.5)]"

dat <- nc_open("./data/hawaii_soest_941a_c863_cca9_a2ed_1d59_9885.nc")

x <- ncvar_get(dat, "longitude")
y <- ncvar_get(dat, "latitude")
hgt <- ncvar_get(dat, "hgt", verbose = F)
dim(hgt) # 144 long, 29 lat, 888 months


# first, extract dates
raw <- ncvar_get(dat, "time") # seconds since 1-1-1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- as.numeric(as.character(years(d)))

# make vectors of lat/long and add (with date) as dimnames
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y)) 


# Change data into a matrix with months / cells for rows / columns
hgt <- aperm(hgt, 3:1)  
hgt <- matrix(hgt, nrow=dim(hgt)[1], ncol=prod(dim(hgt)[2:3]))  

dimnames(hgt) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to confirm that everything is ok
z <- colMeans(hgt, na.rm=T) # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "")

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', add=T, lwd=1)


# limit to 1950:2013
hgt <- hgt[yr %in% 1950:2013,]
m <- m[yr %in% 1950:2013]
yr <- yr[yr %in% 1950:2013]


# get monthly anomalies (remove seasonal signal)

f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(hgt, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

hgt.anom <- hgt - mu   # compute matrix of anomalies

# save for SDE analysis!
write.csv(hgt.anom, "./data/northern.hemisphere.hgt.anom.csv")

# and also save weights for EOF!
weight <- sqrt(cos(lat*pi/180))

# save
write.csv(weight, "./data/north.pacific.hgt.weights.csv")

