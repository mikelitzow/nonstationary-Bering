# evaluate time-dependent slp-sst relationships affecting Bering Sea ecosystem

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

# load mean EBS sst anomaly
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)

# examine cross-correlations
ccf(sst[,2], pc1, lag.max = 10)
print(ccf(sst[,2], pc1, lag.max = 10))
# lags 1-4 are broadly similar


# fit AR model to the entire time series

# make a data frame
yr <- as.numeric(as.character(chron::years(sst$date)))

# and fix incorrect years!
fix <- yr > 2030
yr[fix] <- yr[fix] - 100

m <- as.numeric(months(sst$date))


dat <- data.frame(date = lubridate::parse_date_time(x = paste(yr, m, "01"), orders="ymd", tz="America/Anchorage"),
                sst = sst[,2],
                slp = c(NA, pc1[1:767])) # lagging slp by one month

# and drop NAs
dat <- na.omit(dat)


# ar_ls calculates the process deviations after
# accounting for forcing variables and autocorrelation,
# (1-gamma)
ar_ls = function(time,forcing,gamma) {
  #S(t+1) = (1-GAMMA*DT)*S(t) + F(t)*DT
  forcing = c(forcing - mean(forcing))
  T=length(forcing)
  sig = 0
  
  for(t in 1:(T-1)) {
    #sig[t+1] = -theta*sig[t] + forcing[t]
    sig[t+1] = (1-gamma)*sig[t] + forcing[t]
  }
  
  # next estimates are linearly de-trended
  #s.sig = sig
  sig = sig - lm(sig ~ time)$fitted.values
  # interpolate output on the original time grid
  s.sig=(sig[-1]+sig[-T])/2 # midpoint
  # final step is normalize
  s.sig=s.sig/sd(s.sig)
  return(s.sig)
}

calc_ss = function(theta) {
  pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, theta)
  ss = sum((pred_ts - dat$sst)^2) # return SS for optim (minimizes by default)
}

# optimize by default is minimizing (with maximum = FALSE)
o = optimize(f=calc_ss, interval = c(0,1), maximum=FALSE)

# pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp,
#                 gamma = o$minimum)

pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp,
                gamma = 1/6)

pred.sst = data.frame(t = dat$date,
                      sst = dat$sst,
                      integrated.slp = c(0,-as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP


cor(pred.sst$sst, pred.sst$integrated.slp)


pred.sst <- pred.sst %>%
  mutate(integrated.slp = integrated.slp) %>%
  pivot_longer(cols = -t)

sst.reconstruct <- ggplot(pred.sst, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP", "SST")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.95),
        axis.title.x = element_blank()) +
  ggtitle("EBS SST (r = 0.32)") +
  ylab("Anomaly")

sst.reconstruct

# questions!

# this suggests that tau = 1...seems implausible based on the physics
# can we fix tau?
# just tried fixing tau at 6 months and got r = 0.32!
o$minimum
1/o$minimum
