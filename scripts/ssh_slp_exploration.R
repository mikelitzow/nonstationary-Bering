library(ncdf4)
library(zoo)
library(gplots)
library(dplyr)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(tidyr)
library(nlme)
library(ggplot2)
library(oce)
library(pracma)

# plot settings
theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_f29d_0f00_8b34.nc?ssh[(1950-01-15):1:(2010-12-15T00:00:00Z)][(52.25):1:(62.25)][(180.25):1:(230.25)]

nc <- nc_open("./data/hawaii_soest_f29d_0f00_8b34_2138_dccb_f22e.nc")

# view dates (middle of month):
raw <- ncvar_get(nc, "time")
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))

ssh <- ncvar_get(nc, "ssh") # get all the data!
x <- ncvar_get(nc, "longitude")     # view longitudes (degrees East)
y <- ncvar_get(nc, "latitude")     # view latitudes
lat <- rep(y, length(x))   # Vector of latitudes
lon <- rep(x, each = length(y))   # Vector of longitudes

# process!
ssh <- aperm(ssh, 3:1)  # First, reverse order of dimensions ("transpose" array)

ssh <- matrix(ssh, nrow=dim(ssh)[1], ncol=prod(dim(ssh)[2:3]))  # Change to matrix

dimnames(ssh) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))


# check
SSH.mean <- colMeans(ssh)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)


# limit to GOA
ssh_goa <- ssh

ssh_goa[,lon < 195] <- NA

ssh_goa[,lat > 57 & lon < 201] <- NA

# and dropping the southeast section which seems to dominate EOF1 for godas!
# ssh_goa[,lat < 57 & lon > 220] <- NA

yy <- c(54, 58)
xx <- c(195, 204.5)

mod <- lm(yy ~ xx)

newdata <- data.frame(xx=lon)
predlat <- predict(mod, newdata = newdata)

drop <- NA
for(i in 1:ncol(ssh_goa)){

   drop[i] <- ifelse(lon[i] > 204.5, FALSE,
                     ifelse(lon[i] < 204.5 & lat[i] > predlat[i], TRUE, FALSE))
 }

ssh_goa[,drop] <- NA

# check
SSH.mean <- colMeans(ssh_goa)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# limit to EBS
ssh_ebs <- ssh

ssh_ebs[,lon > 204] <- NA

# ssh_ebs[,lat > 57 & lon < 201] <- NA


yy <- c(52, 58)
xx <- c(188, 204.5)

mod <- lm(yy ~ xx)

newdata <- data.frame(xx=lon)
predlat <- predict(mod, newdata = newdata)

drop <- NA
for(i in 1:ncol(ssh_ebs)){
  # i <- 1
  drop[i] <- ifelse(lon[i] < 204.5 & lat[i] > predlat[i], FALSE, TRUE)
}

ssh_ebs[,drop] <- NA

# check
SSH.mean <- colMeans(ssh_ebs)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# now limit to FMA and fit EOF!
m <- months(d)
yr <- years(d)

ssh_goa_fma <- ssh_goa[m %in% c("Feb", "Mar", "Apr"),]
yr.fma <- yr[m %in% c("Feb", "Mar", "Apr")]

# use SVD on covariance
# need to drop NAs!
land <- is.na(colMeans(ssh_goa_fma))
ssh_goa_fma <- ssh_goa_fma[,!land]

eof <- svd(cov(ssh_goa_fma))

# plot to check
z <- rep(NA, ncol(ssh))
z[!land] <- eof$u[,1]

z <- t(matrix(z,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# and the time series of PC1
goa_ssh_pc1 <- ssh_goa_fma %*% eof$u[,1]

# and annual FMA means for this value
goa_ssh_pc1 <- tapply(goa_ssh_pc1, yr.fma, mean)

# plot to check
plot(1950:2010, goa_ssh_pc1, type="l") # looks right re. 1976/77

## fit eof to ebs data
ssh_ebs_fma <- ssh_ebs[m %in% c("Feb", "Mar", "Apr"),]

# check
SSH.mean <- colMeans(ssh_ebs_fma)
z <- t(matrix(SSH.mean,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)


# use SVD on covariance
# need to drop NAs!
land <- is.na(colMeans(ssh_ebs_fma))
ssh_ebs_fma <- ssh_ebs_fma[,!land]

eof <- svd(cov(ssh_ebs_fma))

# plot to check
z <- rep(NA, ncol(ssh))
z[!land] <- eof$u[,1]

z <- t(matrix(z,length(y)))  
image(x,y,z, col=tim.colors(64), xlim=c(180,250), ylim=c(50,64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,xlim=c(130,250), ylim=c(20,66),add=T, lwd=2)

# and the time series of PC1
ebs_ssh_pc1 <- ssh_ebs_fma %*% eof$u[,1]

# and annual FMA means for this value
ebs_ssh_pc1 <- tapply(ebs_ssh_pc1, yr.fma, mean)

# plot to check
plot(1950:2010, ebs_ssh_pc1, type="l") # looks right re. 1976/77

## slp pc1 ---------------------------------------

# load slp and cell weights
slp <- read.csv("./data/north.pacific.slp.anom.csv", row.names = 1)

# limit to FMA
d <- dates(rownames(slp))
year <- as.numeric(as.character(years(d)))
change <-  year > 2013
year[change] <- year[change]-100

month <- months(d)

slp_fma <- slp[m %in% c("Feb", "Mar", "Apr"),]
yr.fma <- year[m %in% c("Feb", "Mar", "Apr")]

weights <- read.csv("./data/north.pacific.slp.weights.csv", row.names = 1)

pca <- svd.triplet(cov(slp_fma), col.w=weights[,1]) #weighting the columns

pc1_slp <- as.matrix(slp_fma) %*% pca$U[,1]

# and scale!
pc1_slp <- as.vector(scale(pc1_slp))

# and annual FMA means for this value
pc1_slp <- tapply(pc1_slp, yr.fma, mean)

# plot to check
plot(1950:2013, pc1_slp, type="l") # looks right re. 1976/77

# plot
plot_dat <- data.frame(year = 1950:2013,
                       slp_pc1 = pc1_slp) %>%
  dplyr::filter(year %in% 1969:2008)

ssh_dat <- data.frame(year = 1950:2010,
                      ebs_ssh_pc1 = ebs_ssh_pc1,
                      goa_ssh_pc1 = goa_ssh_pc1)

plot_dat <- left_join(plot_dat, ssh_dat) %>%
  pivot_longer(cols = c(-year, -slp_pc1)) %>%
  dplyr::mutate(era = if_else(year <= 1988, "1969-1988", "1989-2008"))

# ggplot(plot_dat, aes(goa_ssh_pc1, ebs_ssh_pc1)) +
#   geom_point() # not unexpected!


ggplot(plot_dat, aes(slp_pc1, value, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~name, scales = "free_y")

## load sst ------------------------

# load regional sst anomaly
sst <- read.csv("./data/regional_monthly_sst.csv")
sst$year <- str_sub(sst$date, 1, -7)
sst$month <- str_sub(sst$date, 6, -4)

sst_fma <- sst %>%
  dplyr::filter(region %in% c("Eastern_Bering_Sea", "Gulf_of_Alaska"),
                month %in% c("02", "03", "04")) %>%
  group_by(region, year) %>%
  summarise(sst_fma = mean(monthly.anom))

# change name
change <- sst_fma$region == "Eastern_Bering_Sea"
sst_fma$region[change] <- "ebs_ssh_pc1"

change <- sst_fma$region == "Gulf_of_Alaska"
sst_fma$region[change] <- "goa_ssh_pc1"

sst_fma <- sst_fma %>%
  dplyr::rename(name = region) %>%
  mutate(year = as.integer(year))

sst_ssh_dat <- left_join(plot_dat, sst_fma)

ggplot(sst_ssh_dat, aes(sst_fma, value, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~name)

## compare sst and wind
trend <- read.csv("./output/wind.ice.recruit.dfa.trends.csv")

# get winter (Nov-Mar) SST
sst_ndjfm <- sst %>%
  dplyr::filter(region == "Eastern_Bering_Sea",
                month %in% c("11", "12", "01", "02", "03")) %>%
  dplyr::mutate(year = as.numeric(year),
                winter.year = if_else(month %in% c("11", "12"), year+1, year)) %>%
  dplyr::group_by(winter.year) %>%
  dplyr::summarise(sst_ndjfm = mean(monthly.anom)) %>%
  filter(winter.year %in% 1951:2013) %>%
  rename(year = winter.year)

# limit to wind and join
trend <- trend %>%
  dplyr::select(-recruit.trend)

dat <- left_join(sst_ndjfm, trend) %>%
  dplyr::filter(year %in% 1951:2008) %>%
  dplyr::mutate(era = case_when(year <= 1968 ~ "1951-1968",
                                year %in% 1969:1988 ~ "1969-1988",
                                year %in% 1989:2008 ~ "1989-2008"))

ggplot(dat, aes(wind.trend, sst_ndjfm, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

mod1 <- nlme::gls(sst_ndjfm ~ wind.trend*era, correlation = corAR1(),
                 data = dat)

summary(mod1)

mod2 <- lm(sst_ndjfm ~ wind.trend*era, data = dat)

summary(mod2)

# rolling windows
dat <- left_join(sst_ndjfm, trend) 

plot_out <- data.frame()

for(i in 1970:2013){
  # i <- 1970
temp <- dat %>%
  filter(year %in% i:(i-19))
  
sst.wind.cor <- cor(temp$sst_ndjfm, temp$wind.trend)
ice.wind.cor <- cor(temp$ice.trend, temp$wind.trend)

plot_out <- rbind(plot_out,
                  data.frame(window.end = i,
                             `Wind-sst correlation` = sst.wind.cor,
                             `Wind-ice correlation` = ice.wind.cor))
}


plot_out <- plot_out %>%
  pivot_longer(cols = -window.end) 

ggplot(plot_out, aes(window.end, value, color = name)) +
  geom_hline(yintercept = 0, color = "dark grey") +
  geom_line() +
  geom_point() + 
  labs(x = "Window end",
       y = "Correlation coefficient") +
  scale_color_manual(values = c("black", cb[7])) +
  scale_y_continuous(breaks = seq(-0.75, 0.75, 0.25)) +
  theme(legend.title = element_blank())

ggsave("./figs/wind_trend_v_ice_sst_rolling_windows.png", width = 6, height = 3, units = 'in')
