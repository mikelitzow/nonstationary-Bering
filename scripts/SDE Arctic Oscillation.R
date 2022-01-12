# evaluate time-dependent AO-sst relationships affecting Bering Sea ecosystem

library(tidyverse)
library(sde)
library(rstan)
library(ggpubr)
library(FactoMineR)

# plot settings
theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# Northern Hemisphere 1000 hPa height EOF1 - Arctic Oscillation

# load hgt and cell weights
hgt.anom <- read.csv("./data/northern.hemisphere.hgt.anom.csv", row.names = 1)
yr <- as.numeric(as.character(chron::years(rownames(hgt.anom))))
fix <- yr > 2030
yr[fix] <- yr[fix]-100

weights <- read.csv("./data/northern.hemisphere.hgt.weights.csv", row.names = 1)

pca <- svd.triplet(cov(hgt.anom), col.w=weights[,1]) #weighting the columns

pc1 <- as.matrix(hgt.anom) %*% pca$U[,1]

# and scale!
pc1 <- as.vector(scale(pc1))

# plot to check
z <- temp.pca$U[,1]

z <- t(matrix(z,length(y)))  
image(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)

# load mean EBS sst anomaly
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)

# examine strongest correlative relationships -
# at various degrees of hgt smoothing and lags

# smooth hgt pc1
pc1.sm.1 <- pc1
pc1.sm.2 <- rollmean(pc1, 2, fill = NA, align = "right")
pc1.sm.3 <- rollmean(pc1, 3, fill = NA, align = "right")
pc1.sm.4 <- rollmean(pc1, 4, fill = NA, align = "right")
pc1.sm.5 <- rollmean(pc1, 5, fill = NA, align = "right")
pc1.sm.6 <- rollmean(pc1, 6, fill = NA, align = "right")
pc1.sm.7 <- rollmean(pc1, 7, fill = NA, align = "right")
pc1.sm.8 <- rollmean(pc1, 8, fill = NA, align = "right")
pc1.sm.9 <- rollmean(pc1, 9, fill = NA, align = "right")

# examine cross-correlations
ccf(sst[,2], pc1.sm.1, lag.max = 10)
ccf(sst[2:768,2], pc1.sm.2[2:768], lag.max = 10)
ccf(sst[3:768,2], pc1.sm.3[3:768], lag.max = 10)
ccf(sst[4:768,2], pc1.sm.4[4:768], lag.max = 10)
ccf(sst[5:768,2], pc1.sm.5[5:768], lag.max = 10)
ccf(sst[6:768,2], pc1.sm.6[6:768], lag.max = 10)
ccf(sst[7:768,2], pc1.sm.7[7:768], lag.max = 10)
ccf(sst[8:768,2], pc1.sm.8[8:768], lag.max = 10)
ccf(sst[9:768,2], pc1.sm.9[9:768], lag.max = 10)

print(ccf(sst[5:768,2], pc1.sm.5[5:768], lag.max = 10))
print(ccf(sst[6:768,2], pc1.sm.6[6:768], lag.max = 10))
print(ccf(sst[7:768,2], pc1.sm.7[7:768], lag.max = 10))
print(ccf(sst[8:768,2], pc1.sm.8[8:768], lag.max = 10))
print(ccf(sst[9:768,2], pc1.sm.9[9:768], lag.max = 10))

# so - very weak relationships - much weaker than N. Pacific SLP!

# compare with ice and wind time series

# load ice and wind dfa trends
wind.ice <- read.csv("./data/wind.ice.dfa.trends.csv")

# also load Oct-Apr SE wind
# this loads heavily on wind DFA trend1 and matches the winter season of peak AO variance

# also load canonical AO to compare with our rolling window EOF
# and scale wind
clim.dat <- read.csv("./data/Bering climate data.csv")
clim.dat <- clim.dat %>%
  dplyr::select(year, AO.jfm, SE.wind.Oct.Apr) %>%
  mutate(SE.wind.Oct.Apr = scale(SE.wind.Oct.Apr))


ggplot(clim.dat, aes(AO.jfm, SE.wind.Oct.Apr, color = year)) + 
  geom_point() 

# loop through 20-year windows and regress ice/wind on EOF1

ice.out <- data.frame()
eof1.load.out <- data.frame()
ice.rolling.window <- data.frame()

se.wind.oct.apr.out <- data.frame()
se.wind.oct.apr.rolling.window <- data.frame()

AO.jfm.out <- data.frame()
AO.jfm.rolling.window <- data.frame()

for(i in 1970:2013){
  # i <- 1970

  temp.hgt <- hgt.anom[yr %in% c((i-20):i),] # year before to year of 20-yr window
  
  temp.pca <- svd.triplet(cov(temp.hgt), col.w=weights[,1]) #weighting the columns
  
  temp.pc1 <- scale(as.matrix(temp.hgt) %*% temp.pca$U[,1])
  
  # save loadings to plot
  eof1.load.out <- rbind(eof1.load.out,
                         data.frame(end.year = i,
                                    loadings = temp.pca$U[,1]))
  # separate out Nov - Mar values
  temp.m <- months(rownames(temp.hgt))
  temp.yr <- as.numeric(as.character(years(rownames(temp.hgt))))

  # correct
  change <- temp.yr > 2030
  temp.yr[change] <- temp.yr[change] - 100
  
  # and advance to winter
  temp.winter <- temp.yr
  temp.winter[temp.m %in% c("Nov", "Dec")] <- temp.winter[temp.m %in% c("Nov", "Dec")] + 1
  
  # winter only
  temp.pc1.winter <- temp.pc1[temp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
  temp.winter <- temp.winter[temp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
  temp.m <- temp.m[temp.m %in% c("Nov", "Dec", "Jan", "Feb", "Mar")]
  
  # winter means
  win.pc1 <- tapply(temp.pc1.winter, temp.winter, mean)
  
  # and remove first and last (incomplete)
  win.pc1 <- win.pc1[2:21]

  # get corresponding ice time series
  ice.temp <- as.vector(wind.ice[wind.ice$year %in% (i-19):i,2])
  
  # save ice.temp and win.pc1 to save 
  ice.rolling.window <- rbind(ice.rolling.window,
                              data.frame(end.year = i, 
                                year = (i-19):i, 
                                ice.trend = ice.temp,
                                height.pc1 = win.pc1))
  

  # and regress
  mod <- lm(ice.temp ~ win.pc1)
  
  ice.out <- rbind(ice.out,
                   data.frame(end.year = i, 
                              coef = coef(mod)[2],
                              r = cor(ice.temp, win.pc1, use ="p"),
                              r.sq = summary(mod)$r.squared))
  
  #####################
  # correlation / regression for SE wind Oct-Apr
  
  # separate out Nov - Mar values
  temp.m <- months(rownames(temp.hgt))
  temp.yr <- as.numeric(as.character(years(rownames(temp.hgt))))
  
  # correct
  change <- temp.yr > 2030
  temp.yr[change] <- temp.yr[change] - 100
  
  # and advance to winter
  temp.winter.oct.apr <- temp.yr
  temp.winter.oct.apr[temp.m %in% c("Oct", "Nov", "Dec")] <- temp.winter.oct.apr[temp.m %in% c("Oct", "Nov", "Dec")] + 1
  
  # winter only
  temp.pc1.winter.oct.apr <- temp.pc1[temp.m %in% c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")]
  temp.winter.oct.apr <- temp.winter.oct.apr[temp.m %in% c("Oct","Nov", "Dec", "Jan", "Feb", "Mar", "Apr")]
  temp.m.oct.apr <- temp.m[temp.m %in% c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr")]
  
  # winter means
  win.pc1.oct.apr <- tapply(temp.pc1.winter.oct.apr, temp.winter.oct.apr, mean)
  
  # and remove first and last (incomplete)
  win.pc1.oct.apr <- win.pc1.oct.apr[2:21]
  
  # get corresponding SE wind oct apr time series
  se.wind.temp <- as.vector(clim.dat[clim.dat$year %in% (i-19):i,3])
  
  # save se.wind.temp and win.pc1.oct.apr 
  se.wind.oct.apr.rolling.window <- rbind(se.wind.oct.apr.rolling.window,
                              data.frame(end.year = i, 
                                         year = (i-19):i, 
                                         se.wind.trend = ice.temp,
                                         height.pc1 = win.pc1.oct.apr))
  
  
  # and regress
  mod <- lm(se.wind.temp ~ win.pc1.oct.apr)
  
  se.wind.oct.apr.out <- rbind(se.wind.oct.apr.out,
                   data.frame(end.year = i, 
                              coef = coef(mod)[2],
                              r = cor(se.wind.temp, win.pc1.oct.apr, use ="p"),
                              r.sq = summary(mod)$r.squared))
  #####################
  # AO.jfm vs. EOF1.jfm correlations
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
  
  # and remove first (don't need as we aren't including winter months from prior calendar year, i.e. Nov and Dec)
  # note that if we want to report this result we may want to refit to years (i-19):i rather than
  # (i - 20):i as is currently done above to accommodate NDJFM at the top of the loop
  
  jfm.pc1 <- jfm.pc1[2:21]
  
  # get corresponding AO.jfm values
  AO.temp <- as.vector(clim.dat[clim.dat$year %in% (i-19):i,2])
  
  # save ice.temp and win.pc1 to save 
  AO.jfm.rolling.window <- rbind(AO.jfm.rolling.window,
                              data.frame(end.year = i, 
                                         year = (i-19):i, 
                                         AO.jfm = AO.temp,
                                         height.pc1 = jfm.pc1))
  # and save correlation
  AO.jfm.out <- rbind(AO.jfm.out,
                   data.frame(end.year = i,
                              r = cor(AO.temp, jfm.pc1, use ="p")))
}

ice.out <- ice.out %>%
  pivot_longer(cols = -end.year)

ggplot(ice.out, aes(end.year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y", ncol = 1)

ggsave("./figs/20_yr_rolling_association_ice_trend_hgt_eof1.png", width = 4, height = 6, units = 'in')

ggplot(AO.jfm.out, aes(end.year, r)) + 
  geom_line() +
  geom_point()


ggsave("./figs/AO-potential_hgt_EOF1_rolling_cor.png", width = 6, height = 4, units = 'in')

se.wind.oct.apr.out  <- se.wind.oct.apr.out %>%
  pivot_longer(cols = -end.year)
  
ggplot(se.wind.oct.apr.out, aes(end.year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y", ncol = 1)

ggsave("./figs/20_yr_rolling_association_SE_wind_oct_apr_hgt_EOF1_rolling_cor.png", width = 4, height = 6, units = 'in')

# save output!
write.csv(ice.out, "./output/geopotential_height_VS.ice_rolling-association.csv", row.names = F)
write.csv(ice.rolling.window, "./output/rolling_window_ice_trend_and_geopotential_hgt_pc1.csv", row.names = F)
write.csv(eof1.load.out, "./output/rolling_window_geopotential_hgt_eof1.csv", row.names = F)

write.csv(se.wind.oct.apr.out, "./output/geopotential_height_vs_se_wind_oct_apr_rolling_association.csv", row.names = F)
write.csv(se.wind.window, "./output/rolling_window_se_wind_and_geopotential_hgt_pc1.csv", row.names = F)

write.csv(AO.jfm.out, "./output/geopotential_height_vs_AO_jfm.csv")

# plot to check
z <- temp.pca$U[,1]

z <- t(matrix(z,length(y)))  
image(x,y,z, col=tim.colors(64))
contour(x, y, z, add=T)  
map('world2Hires',fill=F,add=T, lwd=2)
