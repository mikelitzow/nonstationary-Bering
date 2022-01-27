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
library(FactoMineR)

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


# get monthly anomalies (remove seasonal signal)

f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(hgt, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

hgt.anom <- hgt - mu   # compute matrix of anomalies

# define weights for EOF!
weights <- sqrt(cos(lat*pi/180))

# AO.jfm vs. EOF1.jfm correlations

# load canonical AO
AO.canon <- read.csv("./data/canonical.ao.csv")

AO.canon.jfm <- AO.canon %>%
  dplyr::filter(month %in% 1:3) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(AO.jfm = mean(ao))

AO.jfm.out <- data.frame()
AO.jfm.rolling.window <- data.frame()

for(i in 1970:2021){
  # i <- 1970
  
  temp.hgt <- hgt.anom[yr %in% c((i-19):i),] # year before to year of 20-yr window
  
  temp.pca <- svd.triplet(cov(temp.hgt), col.w=weights) #weighting the columns
  
  temp.pc1 <- scale(as.matrix(temp.hgt) %*% temp.pca$U[,1])
  
  # separate out Jan - Mar values
  temp.m <- months(rownames(temp.hgt))
  temp.yr <- as.numeric(as.character(years(rownames(temp.hgt))))
  
  # correct
  change <- temp.yr > 2030
  temp.yr[change] <- temp.yr[change] - 100
  
  # JFM only
  temp.pc1.jfm <- temp.pc1[temp.m %in% c("Jan", "Feb", "Mar")]
  temp.yr <- temp.yr[temp.m %in% c("Jan", "Feb", "Mar")]
  temp.m <- temp.m[temp.m %in% c("Jan", "Feb", "Mar")]
  
  # winter means
  jfm.pc1 <- tapply(temp.pc1.jfm, temp.yr, mean)
  
  # get corresponding AO.jfm values
  AO.temp <- as.vector(AO.canon.jfm$AO.jfm[AO.canon.jfm$year %in% (i-19):i])
  
  # and save correlation
  AO.jfm.out <- rbind(AO.jfm.out,
                      data.frame(end.year = i,
                                 r = cor(AO.temp, jfm.pc1, use ="p")))
}  

########################

ggplot(AO.jfm.out, aes(end.year, -r)) + 
  geom_line() +
  geom_point()


ggsave("./figs/AO-potential_hgt_EOF1_rolling_cor_1950-2021.png", width = 6, height = 4, units = 'in')