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

# examine strongest correlative relationships -
# at various degrees of slp smoothing and lags

# smooth slp pc1
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

# so - for a first cut, 6 month smooths, lag 1 (i.e., slp averaged over lags -1 : -7)
# appears reasonable

# now we can go back and fit EOF1 to each 20-year rolling window, then fit
# sst-slp SDE model to those data, and save resulting correlation


# first, fit to the entire time series

# make a data frame
yr <- as.numeric(as.character(chron::years(sst$date)))

# and fix incorrect years!
fix <- yr > 2030
yr[fix] <- yr[fix] - 100

m <- as.numeric(months(sst$date))


dat <- data.frame(date = lubridate::parse_date_time(x = paste(yr, m, "01"), orders="ymd", tz="America/Anchorage"),
                sst = sst[,2],
                slp = c(NA, pc1.sm.6[1:863]))

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

pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp,
                gamma = o$minimum)

# save gamma estimate for comparison to Bayes version below
pdo.gamma.ls <- o$minimum

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
  ggtitle("EBS SST (r = 0.27)") +
  ylab("Anomaly")

sst.reconstruct

ggsave("./figs/N.Pacific_slp_pc1_EBS_sst_SDE_1950-2013.png", width = 6, height = 4, units = 'in')


# now plot / compare the low-frequency correspondence

pred.sst.sm = data.frame(t = dat$date,
                      sst = rollmean(dat$sst, 13, fill = NA),
                      integrated.slp = rollmean(c(0,-as.numeric(pred_ts)), 13, fill = NA)) ## NB - reversing the sign of integrated SLP

cor(pred.sst.sm$sst, pred.sst.sm$integrated.slp, use = "p")

pred.sst.sm <- pred.sst.sm %>%
  mutate(integrated.slp = integrated.slp) %>%
  pivot_longer(cols = -t)

sst.reconstruct.sm <- ggplot(pred.sst.sm, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP", "SST")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.95),
        axis.title.x = element_blank()) +
  ggtitle("EBS SST - 13 month smooths (r = 0.46)") +
  ylab("Anomaly")

sst.reconstruct.sm

ggsave("./figs/N.Pacific_slp_pc1_EBS_sst_SDE_1950-2013_13-month_smooths.png", width = 6, height = 4, units = 'in')

# now loop through on 240-month (20 year) rolling windows

# make object to capture correlation
cor.out <- data.frame()

for(i in 246:nrow(sst)){
  # i <- 247
  # fit eof to subset of slp data 
  
  # parcel out the 240 months ending at time == i
  temp.slp <- slp[(i-245):(i-1),] 
  
  # and fit EOF
  pca.temp <- svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
  pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]
  
  # scale
  pc1.temp <- as.vector(scale(pc1.temp))
  # names(pc1.temp) <- row.names(temp.slp) # just used this to check the setup!
  # and smooth
  pc1.temp.sm <- rollmean(pc1.temp, 6, fill = NA, align = "right")
  
  # now select sst, lagged by 1 month wrt pc1.temp.sm
  temp.sst <- sst[(i-239):i,]
  
  temp.yr <- as.numeric(as.character(chron::years(temp.sst[nrow(temp.sst),1])))
  temp.m <- as.numeric(months(temp.sst[nrow(temp.sst),1]))
  
  cor.out <- rbind(cor.out,
                   data.frame(end.date = lubridate::parse_date_time(x = paste(temp.yr, temp.m, "01"), orders="ymd", tz="America/Anchorage"),
                              cor = cor(temp.sst[,2], pc1.temp.sm[6:245])))
  
}

# add dates in new format to facilitate plotting
cor.out$end.year <- as.numeric(as.character(chron::years(cor.out$end.date)))
cor.out$end.month <- as.numeric(as.factor(base::months(cor.out$end.date)))

ggplot(cor.out, aes(end.date, -cor)) +
  geom_line() +
scale_x_continuous(labels = c("1951-1970", 
                              "1961-1980", 
                              "1971-1990",
                              "1981-2000",
                              "1991-2010",
                              "2001-2020"))

pdo.cors <- data.frame()

for(i in 121:760){
  # i <- 121
  temp <- pred.pdo[((i-120):(i+120)),]
  
  pdo.cors <- rbind(pdo.cors,
                    data.frame(cor = cor(temp$sst.pc1, temp$integrated.slp)))
}

pdo.cors$date <- pred.pdo$t[121:760]

pdo.cor.plot <- ggplot(pdo.cors, aes(date, cor)) +
  geom_line() +
  ylab("Correlation (r)") +
  theme(axis.title.x = element_blank()) +
  coord_cartesian(xlim = range(pred.pdo$t)) +
  ggtitle("20-year moving window correlations")

pdo.cor.plot

ggplot(pred.pdo, aes(integrated.slp, sst.pc1)) + 
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = F)

pred.pdo <- pred.pdo %>%
  mutate(integrated.slp = integrated.slp) %>%
  pivot_longer(cols = -t)

pdo.reconstruct <- ggplot(pred.pdo, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP PC2", "SST PC2")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.95),
        axis.title.x = element_blank()) +
  ggtitle("PDO (r = 0.32)") +
  ylab("Anomaly")

pdo.reconstruct



####################################################################
####################################################################
# now loop through the 21-year windows and calculate relevant values
# note that I am using 253-month (21 year + 1 month) windows to ease plotting at window centers

for(i in 127:(nrow(ind.weighted.21)-126)){
  
  temp.sst <- X1.anom.detr[(i-126):(i+126),] # using anomalies here, NOT scaled, in order to capture importance of variance
  temp.pca <- svd.triplet(cov(temp.sst), col.w=weight) #weighting the columns
  pc1 <- temp.sst %*% temp.pca$U[,1]

# load data
d = read.csv("data/slp_sst_PCs_1948-2021.csv",
             stringsAsFactors = FALSE)

# these data are PC1/PC2 of NE Pacific SLPa/SSTa for 1/1948 - 5/2021
# each PC time series has already been scaled

# plot to check
plot.dat <- d %>%
  pivot_longer(cols = c(-m, -month, -year, -variable)) %>%
  mutate(decimal.year = year + (m-0.5)/12)

ggplot(plot.dat, aes(decimal.year, value)) +
  geom_line() +
  facet_grid(name~variable)

# looks as we would expect - white noise for the SLP PCs, red noise for the SST PCs (especially PC1)

# find the decorrelation scale for each SST PC to estimate theta

print(acf(filter(d, variable == "sst")$pc1)) # above 0.5 for lags 1-7
print(acf(filter(d, variable == "sst")$pc2)) # above 0.5 for lags 1-3

# and find peak cross correlation for PC1s / PC2s
print(ccf(filter(d, variable == "slp")$pc1, filter(d, variable == "sst")$pc1)) # peak at lag1
print(ccf(filter(d, variable == "slp")$pc2, filter(d, variable == "sst")$pc2)) # also peak at lag1; much weaker


# attempt to reconstruct SST PC1
d$date = lubridate::parse_date_time(x = paste(d$year,d$month,"01"),orders="ymd",tz="Pacific")

# make a bespoke data set with lagged slp
sst <- filter(d, variable == "sst")
slp <- filter(d, variable == "slp")

dat <- data.frame(date = sst$date[2:nrow(sst)],
                  sst = sst$pc1[2:nrow(sst)],
                  slp = slp$pc1[1:nrow(slp)-1])

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
  ss = sum((pred_ts - dat$sst[-1])^2) # return SS for optim (minimizes by default)
}

# optimize by default is minimizing (with maximum = FALSE)
o = optimize(f=calc_ss, interval = c(0,1), maximum=FALSE)

pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp,
                gamma = o$minimum)

# save gamma estimate for comparison to Bayes version below
pdo.gamma.ls <- o$minimum

pred.pdo = data.frame(t = dat$date,
                      sst.pc1 = dat$sst,
                      integrated.slp = c(0,-as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP


cor(pred.pdo$sst.pc1, pred.pdo$integrated.slp)

pdo.cors <- data.frame()

for(i in 121:760){
  # i <- 121
  temp <- pred.pdo[((i-120):(i+120)),]
  
  pdo.cors <- rbind(pdo.cors,
                     data.frame(cor = cor(temp$sst.pc1, temp$integrated.slp)))
}

pdo.cors$date <- pred.pdo$t[121:760]

pdo.cor.plot <- ggplot(pdo.cors, aes(date, cor)) +
  geom_line() +
  ylab("Correlation (r)") +
  theme(axis.title.x = element_blank()) +
  coord_cartesian(xlim = range(pred.pdo$t)) +
  ggtitle("20-year moving window correlations")

pdo.cor.plot

ggplot(pred.pdo, aes(integrated.slp, sst.pc1)) + 
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = F)

pred.pdo <- pred.pdo %>%
  mutate(integrated.slp = integrated.slp) %>%
  pivot_longer(cols = -t)

pdo.reconstruct <- ggplot(pred.pdo, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP PC2", "SST PC2")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.95),
        axis.title.x = element_blank()) +
  ggtitle("PDO (r = 0.32)") +
  ylab("Anomaly")

pdo.reconstruct

png("./figs/PDO_reconstruction_and_moving_windows.png", width=6, height=8, units='in', res=300)
ggarrange(pdo.reconstruct, pdo.cor.plot, ncol = 1, labels = "auto")
dev.off()


### compare with NPGO ---------------------
dat <- data.frame(date = sst$date[2:nrow(sst)],
                  sst = sst$pc2[2:nrow(sst)],
                  slp = slp$pc2[1:nrow(slp)-1])

calc_ss = function(theta) {
  pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, theta)
  ss = sum((pred_ts - dat$sst[-1])^2) # return -SS for optim
}

o = optimize(f=calc_ss, interval = c(0,1))

pred_ts = ar_ls(1:nrow(dat), forcing=dat$slp, gamma = o$minimum)

# save gamma estimate for comparison to Bayes version below
npgo.gamma.ls <- o$minimum

pred.npgo = data.frame(t = dat$date,
                      sst.pc2 = dat$sst,
                      integrated.slp = c(0, -as.numeric(pred_ts))) ## NB - reversing the sign of integrated SLP


cor(pred.npgo$sst.pc2, -pred.npgo$integrated.slp)

nrow(pred.npgo)

npgo.cors <- data.frame()

for(i in 121:760){
  # i <- 121
  temp <- pred.npgo[((i-120):(i+120)),]
  
  npgo.cors <- rbind(npgo.cors,
                     data.frame(cor = cor(temp$sst.pc2, -temp$integrated.slp)))
}

npgo.cors$date <- pred.npgo$t[121:760]

npgo.cor.plot <- ggplot(npgo.cors, aes(date, cor)) +
  geom_line() +
  ylab("Correlation (r)") +
  theme(axis.title.x = element_blank()) +
  coord_cartesian(xlim = range(pred.npgo$t)) +
  ggtitle("20-year moving window correlations")

npgo.cor.plot

ggplot(pred.npgo, aes(-integrated.slp, sst.pc2)) + # also not linear
  geom_point() +
  geom_smooth(method = "gam", formula = y ~ s(x, k=4), se = F)

pred.npgo <- pred.npgo %>%
  mutate(integrated.slp = -integrated.slp) %>%
  pivot_longer(cols = -t)

npgo.reconstruct <- ggplot(pred.npgo, aes(t, value, color = name)) +
  geom_line() +
  scale_color_manual(values = cb[c(2,6)], labels = c("Integrated SLP PC2", "SST PC2")) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.15),
        axis.title.x = element_blank()) +
  ggtitle("NPGO (r = 0.35)") +
  ylab("Anomaly")

npgo.reconstruct

png("./figs/NPGO_reconstruction_and_moving_windows.png", width=6, height=8, units='in', res=300)
ggarrange(npgo.reconstruct, npgo.cor.plot, ncol = 1, labels = "auto")
dev.off()


png("./figs/PDO_NPGO_reconstructions.png", width=6, height=8, units='in', res=300)
ggarrange(pdo.reconstruct, npgo.reconstruct, ncol=1)
dev.off()


### Bayesian version --------------------
# M = 1 here, because the forcing variable obs_x
# is already discretized and there's no missing values

# redefine data_list to make it clear what we're working on
dat <- data.frame(date = sst$date[2:nrow(sst)],
                  sst = sst$pc2[2:nrow(sst)],
                  slp = slp$pc2[1:nrow(slp)-1])

data_list = list(
  M = 1,
  N = nrow(dat),
  obs_y = c(dat$sst),
  obs_x = c(dat$slp)
)

# first Bayesian model is estimating both gamma, obs_sigma, and sigma,
# which scales the forcing.
library(bayesplot)
fit = stan("code/ar1_forcing_ss.stan",
           data = data_list,
           pars = c("sigma","pred_y","gamma","obs_sigma"),
           iter=5000,
           chains=3)
pars = rstan::extract(fit)

mcmc_areas(fit, pars = c("sigma","gamma","obs_sigma"))

## Second, we can fix the process sd, "sigma" at 1 and re-fit
fit2 = stan("code/ar1_forcing_ss_fixsig.stan",
           data = data_list,
           pars = c("pred_y","gamma","obs_sigma"),
           iter=5000,
           chains=3)
pars = rstan::extract(fit2)
mcmc_areas(fit2, pars = c("gamma","obs_sigma"))

## Third, if we knew the observation error standard deviation of these SST
# measurements, we could pass that in as a known quantity and only estimate
# 'gamma' -- this is most directly similar to the SS approach in the optim
# function above. Note: the mean estimates of the trajectory won't change if
# you vary obs_sigma, but the precision of the estimates will
data_list$obs_sigma = 0.1
fit3 = stan("code/ar1_forcing_ss_fixbothsig.stan",
            data = data_list,
            pars = c("pred_y","gamma"),
            iter=5000,
            chains=3)
pars = rstan::extract(fit3)

# save gamma estimate for comparison to ls approach
npgo.gamma.bayes <- median(pars$gamma)

mcmc_areas(fit3, pars = c("gamma"))
ggsave("./figs/gamma_npgo_model3.png", width = 6, height = 4, units = 'in')


## Now we can also compare estimates across the three models
library(broom.mixed)
y1 = tidy(fit, pars=c("pred_y"))
y1$time = seq(1,nrow(y1))
y1$model = "ss"
y2 = tidy(fit2, pars=c("pred_y"))
y2$time = seq(1,nrow(y2))
y2$model = "fixsig"
y3 = tidy(fit3, pars=c("pred_y"))
y3$time = seq(1,nrow(y3))
y3$model = "fixboth"
png("figs/bayes_estimates.png", width=6, height=8, units='in', res=300)
rbind(y1,y2,y3) %>%
  ggplot(aes(time, estimate,fill=model,col=model)) +
  geom_ribbon(aes(ymin=estimate-std.error,ymax=estimate+std.error),alpha=0.1) +
  geom_line(alpha=0.1)
dev.off()

# look at the std error rather than the estimate. this shows that
png("figs/bayes_sds_of_estimates.png", width=6, height=8, units='in', res=300)
rbind(y1,y2,y3) %>%
  ggplot(aes(time, log(std.error),fill=model,col=model)) +
  #geom_ribbon(aes(ymin=estimate-std.error,ymax=estimate+std.error),alpha=0.1) +
  geom_line(alpha=0.5) +
  ylab("Ln(std.error of estimate")
dev.off()

# recalculate correlations with Bayesian estimates
pred.npgo = data.frame(t = dat$date,
                      sst.pc2 = dat$sst,
                      integrated.slp = c(0,-as.numeric(pred_ts)))
cor(pred.npgo$sst.pc2, -y1$estimate)
cor(pred.npgo$sst.pc2, -y2$estimate)
cor(pred.npgo$sst.pc2, -y3$estimate)

## Compare with Bayesian model of PDO -------------------------------

dat <- data.frame(date = sst$date[2:nrow(sst)],
                  sst = sst$pc1[2:nrow(sst)],
                  slp = slp$pc1[1:nrow(slp)-1])

data_list = list(
  M = 1,
  N = nrow(dat),
  obs_y = c(dat$sst),
  obs_x = c(dat$slp)
)

# fit to model 3 - closest to SS approach

data_list$obs_sigma = 0.1

pdofit3 = stan("code/ar1_forcing_ss_fixbothsig.stan",
            data = data_list,
            pars = c("pred_y","gamma"),
            iter=5000,
            chains=3)
pars = rstan::extract(pdofit3)

# save gamma estimate for comparison to ls approach
pdo.gamma.bayes <- median(pars$gamma)

# and plot this comparison
plot.gamma <- data.frame(pdo_least_squares = pdo.gamma.ls,
                         pdo_bayes = pdo.gamma.bayes,
                         npgo_least_squares = npgo.gamma.ls,
                         npgo_bayes = npgo.gamma.bayes)

plot.gamma <- pivot_longer(plot.gamma, cols = 1:4)


ggplot(plot.gamma, aes(name, value)) +
  geom_col()

ggsave("./figs/gamma_estimates_ls_bayes.png", width = 5, height = 4, units = 'in')



mcmc_areas(pdofit3, pars = c("gamma"))
ggsave("./figs/gamma_pdo_model3.png", width = 6, height = 4, units = 'in')

y3 = tidy(pdofit3, pars=c("pred_y"))
y3$time = seq(1,nrow(y3))


# recalculate correlations with Bayesian estimates
pred.pdo = data.frame(t = dat$date,
                       sst.pc1 = dat$sst,
                       integrated.slp = c(0,-as.numeric(pred_ts)))

cor(pred.pdo$sst.pc1, y3$estimate)

