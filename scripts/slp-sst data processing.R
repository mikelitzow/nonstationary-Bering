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
# library(sinkr)


# to load sinkr, use package devtools and uncomment the following:
 # devtools::install_github("marchtaylor/sinkr")

# load slp
# identify latest year and month needed
year <- 2021
month <- "12"
query <- c("f19d_3925_d70b.nc?slp")

variable <- "slp"

URL <- paste("http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_", query, "[(1948-01-01):1:(", year, "-",
               month, "-01T00:00:00Z)][(19.99970054626):1:(69.52169799805)][(120):1:(249.375)]", sep="")

# paste URL into browser to download

dat <- nc_open("data/hawaii_soest_f19d_3925_d70b_f63c_43ff_14ca.nc")

x <- ncvar_get(dat, "longitude")
y <- ncvar_get(dat, "latitude")
slp <- ncvar_get(dat, "slp", verbose = F)
dim(slp) # 53 long, 21 lat, 888 months

# need to reverse latitude for plotting!
# y <- rev(y)
# slp <- slp[,21:1,]

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
slp <- aperm(slp, 3:1)  
slp <- matrix(slp, nrow=dim(slp)[1], ncol=prod(dim(slp)[2:3]))  

dimnames(slp) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to confirm that everything is ok
z <- colMeans(slp, na.rm=T) # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "")

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', add=T, lwd=1)


# limit to 1950:2013
slp <- slp[yr %in% 1950:2013,]
m <- m[yr %in% 1950:2013]
yr <- yr[yr %in% 1950:2013]


# let's get a regional PCA for the area affecting the GOA and Bering more or less directly
regional.slp <- slp
poly.x <- c(139, 139, 241, 241, 139) 
poly.y <- c(29, 68, 68, 29, 29)

# get x and y for plotting later

write.csv(x[x >=139 & x <= 241], "./output/slp_x.csv", row.names = F)


write.csv(y[y >=29 & y <= 68], "./output/slp_y.csv", row.names = F)


xp <- cbind(poly.x, poly.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

regional.slp[,!check] <- NA

# plot to check

z <- colMeans(regional.slp, na.rm=T) # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", xlim=c(100,260), ylim=c(30,70))

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', add=T, lwd=1)


# regional subset

# identify columns containing NA
out <- is.na(colMeans(regional.slp)) 

X <- regional.slp[,!out] 

# get anomalies

f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(X, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

slp.anom <- X - mu   # compute matrix of anomalies

# save for SDE analysis!
write.csv(slp.anom, "./data/north.pacific.slp.anom.csv")

# and also save weights for EOF!
# get a vector of weights (square root of the cosine of latitude)

temp1 <- str_split(colnames(X), "E", simplify = T)[,1]
X.lats <- as.numeric(str_split(temp1, "N", simplify = T)[,2])

weight <- sqrt(cos(X.lats*pi/180))

# save
write.csv(weight, "./data/north.pacific.slp.weights.csv")

## extract SST
####################
# add ERSST v5

## load and process ERSST ------------------------

# download.file("https://coastwatch.pfeg.noaa.gov/erddap/griddap/nceiErsstv5.nc?sst[(1950-01-01):1:(2013-12-01T00:00:00Z)][(0.0):1:(0.0)][(20):1:(66)][(110):1:(250)]", "~temp")



# load and process SST data
# nc <- nc_open("~temp")

nc <- nc_open("./data/nceiErsstv5_b279_4822_ef14.nc")

# process

ncvar_get(nc, "time")   # seconds since 1-1-1970
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- years(d)

x <- ncvar_get(nc, "longitude")
y <- ncvar_get(nc, "latitude")

SST <- ncvar_get(nc, "sst", verbose = F)

SST <- aperm(SST, 3:1)

SST <- matrix(SST, nrow=dim(SST)[1], ncol=prod(dim(SST)[2:3]))

# Keep track of corresponding latitudes and longitudes of each column:
lat <- rep(y, length(x))
lon <- rep(x, each = length(y))
dimnames(SST) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to check

temp.mean <- colMeans(SST, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "")
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR', 'Japan', 'Mexico', 'South Korea', 'North Korea', 'China', 'Mongolia'), fill=T,add=T, lwd=1, col="lightyellow3")

# extract regional areas
# EBS polygon
ebs.x <- c(183, 183, 203, 203, 191) 
ebs.y <- c(53, 65, 65, 57.5, 53)

polygon(ebs.x, ebs.y, border = "red", lwd = 2)

# GOA polygon
goa.x <- c(201, 201, 205, 208, 225, 231, 201)
goa.y <- c(55, 56.5, 59, 61, 61, 55, 55)

polygon(goa.x, goa.y, border = "red", lwd = 2)

# BC polygon
bc.x <- c(231, 238, 228, 225, 225)
bc.y <- c(55, 49, 49, 53, 55)

polygon(bc.x, bc.y, border = "red", lwd = 2)

# northern CCE polygon
ncc.x <- c(238, 238, 233, 233, 238)
ncc.y <- c(49, 41, 41, 49, 49)

polygon(ncc.x, ncc.y, border = "red", lwd = 2)

# southern CCE polygon
scc.x <- c(238, 243, 237, 233, 233, 238)
scc.y <- c(41, 33, 33, 39, 41, 41)

polygon(scc.x, scc.y, border = "red", lwd = 2)

# those  areas look fine

# define cells within each polygon and plot to check

#first, ebs
ebs.sst <- as.data.frame(SST)

xp <- cbind(ebs.x, ebs.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ebs.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# now, goa

goa.sst <- as.data.frame(SST)

xp <- cbind(goa.x, goa.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

goa.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(goa.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# BC coast

bc.sst <- as.data.frame(SST)

xp <- cbind(bc.x, bc.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

bc.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(bc.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(170, 240), ylim = c(40, 66))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")


# NCC 

ncc.sst <- as.data.frame(SST)

xp <- cbind(ncc.x, ncc.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ncc.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(ncc.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(210, 245), ylim = c(30, 52))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")

# SCC
scc.sst <- as.data.frame(SST)

xp <- cbind(scc.x, scc.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)
scc.sst[,!check] <- NA

# plot to check
temp.mean <- colMeans(scc.sst, na.rm=T)
z <- t(matrix(temp.mean,length(y)))
image.plot(x,y,z, col=oceColorsPalette(64), xlab = "", ylab = "", xlim = c(210, 245), ylim = c(30, 52))
contour(x, y, z, add=T)
map('world2Hires',c('Canada', 'usa', 'USSR'), fill=T,add=T, lwd=1, col="lightyellow3")


## summarize time series for each region ------------------------------

# first, make a list of data frames

# npac is the full grid

npac.sst <- as.data.frame(SST)

sst.data.frames <- list()
sst.data.frames[[1]] <- npac.sst
sst.data.frames[[2]] <- ebs.sst
sst.data.frames[[3]] <- goa.sst
sst.data.frames[[4]] <- bc.sst
sst.data.frames[[5]] <- ncc.sst
sst.data.frames[[6]] <- scc.sst

# save names of each df for processing

sst.data.names <- c("npac_sst",
                    "ebs_sst",
                    "goa_sst",
                    "bc_sst",
                    "ncc_sst",
                    "scc_sst")

# and a vector of clean names
sst.clean.names <- c("North_Pacific",
                     "Eastern_Bering_Sea",
                     "Gulf_of_Alaska",
                     "British_Columbia_Coast",
                     "Northern_California_Current",
                     "Southern_California_Current")

# loop through each df, process, summarize, combine, save

# create vector of latitudes to weight mean sst by cell area

cell.weight <- sqrt(cos(lat*pi/180))

# plot to check
hist(cell.weight, breaks = 50)
unique(cell.weight) # looks right

# create a weighted mean function
weighted.cell.mean <- function(x) weighted.mean(x, w = cell.weight, na.rm = T)

# create a function to compute monthly anomalies
monthly.anomalies <- function(x) tapply(x, m, mean) 


# create blank data frame for catching results
temp.anomaly.time.series <- data.frame()

# processing loop

for(i in 1: length(sst.data.names)){
  
  # i <- 1
  # pull out sst dataframe of interest
  temp.dat <- sst.data.frames[[i]]
  
  # first, calculate monthly anomalies
  mu <- apply(temp.dat, 2, monthly.anomalies)	# compute monthly means for each time series (cell)
  mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location
  
  temp.anom <- temp.dat - mu   # compute matrix of anomalies
  
  # calculate weighted monthly means
  temp.monthly.anom <- apply(temp.anom, 1, weighted.cell.mean)
  
  # add to df
  temp.anomaly.time.series <- rbind(temp.anomaly.time.series,
                                    data.frame(region = sst.clean.names[i],
                                    monthly.anom = temp.monthly.anom,
                                    date = lubridate::parse_date_time(x = paste(yr,as.numeric(m),"01"),orders="ymd",tz="America/Anchorage")))
  
}

ggplot(temp.anomaly.time.series, aes(date, monthly.anom, color = region)) +
  geom_line(lwd = 0.02) +
  ylab("Monthly anomaly") +
  theme(axis.title.x = element_blank())
  
# and save
write.csv(temp.anomaly.time.series, "./data/regional_monthly_sst.csv", row.names = F)
  
