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

# so - for a first cut, 6 month smooths, lag 1 (i.e., slp averaged over lags -1 : -6)
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
                slp = c(NA, pc1.sm.6[1:767])) # lagging slp - i.e., slp[6] corresponds to slp[7]

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

# make object to capture correlation and regression output
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
  pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!
  # names(pc1.temp) <- row.names(temp.slp) # just used this to check the setup!
  # and smooth
  pc1.temp.sm <- rollmean(pc1.temp, 6, fill = NA, align = "right")
  
  # now select sst, lagged by 1 month wrt pc1.temp.sm
  temp.sst <- sst[(i-239):i,]
  
  temp.yr <- as.numeric(as.character(chron::years(temp.sst[nrow(temp.sst),1])))
  temp.m <- as.numeric(months(temp.sst[nrow(temp.sst),1]))
  
  mod <- lm(temp.sst[,2] ~ pc1.temp.sm[6:245])
  
  cor.out <- rbind(cor.out,
                   data.frame(end.date = lubridate::parse_date_time(x = paste(temp.yr, temp.m, "01"), orders="ymd", tz="America/Anchorage"),
                              coef = coef(mod)[2],
                              r = cor(temp.sst[,2], pc1.temp.sm[6:245]),
                              r.sq = summary(mod)$r.squared))
  

}


cor.out <- cor.out %>%
  pivot_longer(cols = -end.date)

ggplot(cor.out, aes(end.date, value)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)
  

ggsave("./figs/20_year_SLP_SDE_vs_SST.png", width = 4, height = 6, units = 'in')


# and we really just want to plot r values
r.plot <- cor.out %>%
  filter(name == "r")

ggplot(r.plot, aes(end.date, value)) +
  geom_line() +
  labs(x = "End year", 
       y = "r")

ggsave("./figs/20_year_SLP_SDE_vs_SST_correlation.png", width = 6, height = 4, units = 'in')

## fit Bayesian regression to correlation time series -------------------------
# (can we distinguish the trend from 0?)
library(plyr)
library(rstan)
library(brms)
library(bayesplot)
library(bayesdfa)
source("./scripts/stan_utils.R")

# subset data
dat <- cor.out %>%
  dplyr::filter(name == "r") %>%
  dplyr::select(-name) %>%
  dplyr::rename(r = value) %>%
  mutate(year = lubridate::year(end.date),
         month = lubridate::month(end.date)) 

# now identify a winter year corresponding to Oct - Sept
dat$winter.year <- if_else(dat$month %in% 10:12, dat$year + 1, dat$year)

dat <- dat %>%
  dplyr::group_by(winter.year) %>%
  dplyr::summarise(r = mean(r))

# check to plot
# dat$dec.yr <- dat$year + (dat$month-0.5)/12
# 
# ggplot(dat, aes(dec.yr, r)) +
#   geom_line() +
#   geom_point() 

ggplot(dat, aes(winter.year, r)) +
  geom_line() +
  geom_point()
  
# # define model formula
# dat$month_fac <- as.factor(dat$month)

# slp_sst_cor_formula <-  bf(r ~ s(year, k = 5) + (1 | month_fac/year) + ar(p = 1)) # limiting k to ~ number of decades in time series
# 
# # show default priors
# get_prior(slp_sst_cor_formula, dat)
# 
# # set priors
# priors <- c(set_prior("normal(0, 3)", class = "ar"),
#             set_prior("normal(0, 3)", class = "b"),
#             set_prior("normal(0, 3)", class = "Intercept"),
#             set_prior("student_t(3, 0, 3)", class = "sd"),
#             set_prior("student_t(3, 0, 3)", class = "sds"),
#             set_prior("student_t(3, 0, 3)", class = "sigma"))

slp_sst_cor_formula <-  bf(r ~ s(winter.year, k = 4) + ar(p = 1)) # limiting k to ~ number of decades in time series

# show default priors
get_prior(slp_sst_cor_formula, dat)

# set priors
priors <- c(set_prior("normal(0, 3)", class = "ar"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"))


# fit model 
sst_slp_cor_brm <- brm(slp_sst_cor_formula,
                    prior = priors,
                    data = dat,
                    cores = 4, chains = 4, iter = 5000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99999, max_treedepth = 10))

codR_dfa_brm  <- add_criterion(codR_dfa_brm, c("loo", "bayes_R2"), moment_match = TRUE)
saveRDS(codR_dfa_brm, file = "output/codR_dfa_brm.rds")

codR_dfa_brm <- readRDS("./output/codR_dfa_brm.rds")
check_hmc_diagnostics(sst_slp_cor_brm$fit)
neff_lowest(sst_slp_cor_brm$fit)
rhat_highest(sst_slp_cor_brm$fit)
summary(sst_slp_cor_brm)
bayes_R2(sst_slp_cor_brm)
plot(sst_slp_cor_brm$criteria$loo, "k")
plot(conditional_effects(sst_slp_cor_brm), ask = FALSE)

 ## fit model without ar()
slp_sst_cor_formula <-  bf(r ~ s(winter.year, k = 4)) # limiting k to ~ number of decades in time series

# show default priors
get_prior(slp_sst_cor_formula, dat)

# set priors
priors <- c(set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 3)", class = "Intercept"),
            set_prior("student_t(3, 0, 3)", class = "sds"),
            set_prior("student_t(3, 0, 3)", class = "sigma"))


# fit model 
sst_slp_cor_brm <- brm(slp_sst_cor_formula,
                       prior = priors,
                       data = dat,
                       cores = 4, chains = 4, iter = 3000,
                       save_pars = save_pars(all = TRUE),
                       control = list(adapt_delta = 0.9999, max_treedepth = 12))


check_hmc_diagnostics(sst_slp_cor_brm$fit)
neff_lowest(sst_slp_cor_brm$fit)
rhat_highest(sst_slp_cor_brm$fit)
summary(sst_slp_cor_brm)
bayes_R2(sst_slp_cor_brm)
plot(sst_slp_cor_brm$criteria$loo, "k")
plot(conditional_effects(sst_slp_cor_brm), ask = FALSE)

# so these models only show a trend if ar(1) term isn't used!
# either that means a trend can't be demonstrated given the autocorrelation in the data
# or I'm asking too much of the model to effectively use the same 
# term as a covariate and an AR term
# might be a bad idea to fit any statistical models to these highly overlapping windows!



## plot sst-slp regression maps for contrasting eras -----------------------------
# 1950-1968, 1969-1988, 1989-2008

library(maps)
library(maptools)
library(mapdata)
library(fields)

# load slp and cell weights
slp <- read.csv("./data/north.pacific.slp.anom.csv", row.names = 1)

# get a dates vector to deal with!
d <- chron::dates(rownames(slp))
yr <- as.numeric(as.character(chron::years(d)))

# fix bad years
change <- yr > 2030
yr[change] <- yr[change]-100

# get y and z
x <- seq(152.5, 230, 2.5)
y <- seq(40, 62.5, 2.5)

# load mean EBS sst anomaly
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)
sst$year <- as.numeric(as.character(chron::years(sst$date)))

# fix bad year values
fix <- sst$year > 2030
sst$year[fix] <- sst$year[fix] - 100

# select the three sst windows (response)
sst1 <- sst %>%
  filter(year %in% 1950:1968)

sst2 <- sst %>%
  filter(year %in% 1969:1988)

sst3 <- sst %>%
  filter(year %in% 1989:2008)

# select the corresponding slp windows (1 year behind to allow smoothing!)
slp1 <- slp[yr %in% 1950:1968,]

slp2 <- slp[yr %in% 1968:1988,]

slp3 <- slp[yr %in% 1988:2008,]

# smooth and scale slp
ff <- function(x) as.vector(scale(rollmean(x, 6, fill = NA, align = "right")))
slp.sm.1 <- apply(slp1, 2, ff)
d1 <- chron::dates(rownames(slp))[yr %in% 1950:1968]
sst1$date

# remove first sst value to accommodate lag
sst1 <- sst1[2:nrow(sst1),]

# and last row of slp1
slp.sm.1 <- slp.sm.1[1:(nrow(slp.sm.1)-1),]

# now loop through slp columns, fit lm, and keep coef
coef1 <- NA

for(i in 1:ncol(slp.sm.1)){
  # i <- 1
  mod <- lm(sst1$anom ~ slp.sm.1[,i])
  coef1[i] <- mod$coefficients[2]
}

## second window
# smooth and scale slp
slp.sm.2 <- apply(slp2, 2, ff)
d2 <- chron::dates(rownames(slp))[yr %in% 1968:1988]
sst2$date

# trim slp.sm.2
slp.sm.2 <- slp.sm.2[12:(nrow(slp.sm.2)-1),]

# now loop through slp columns, fit lm, and keep coef
coef2 <- NA

for(i in 1:ncol(slp.sm.2)){
  # i <- 1
  mod <- lm(sst2$anom ~ slp.sm.2[,i])
  coef2[i] <- mod$coefficients[2]
}

## third window
# smooth and scale slp
slp.sm.3 <- apply(slp3, 2, ff)
d3 <- chron::dates(rownames(slp))[yr %in% 1988:2008]
sst3$date
d3
# trim slp.sm.2
slp.sm.3 <- slp.sm.3[12:(nrow(slp.sm.3)-1),]

# now loop through slp columns, fit lm, and keep coef
coef3 <- NA

for(i in 1:ncol(slp.sm.3)){
  # i <- 1
  mod <- lm(sst3$anom ~ slp.sm.3[,i])
  coef3[i] <- mod$coefficients[2]
}

#plot

# get range
range <- range(coef1, coef2, coef3)
zlim <- c(range[1], -range[1])

# draw sebs box for sst area
box.x <- c(187, 199, 199, 203, 203, 199, 199, 193, 193, 187)
box.y <- c(61, 61, 59, 59, 57, 57, 55, 55, 53, 53)

png("./figs/era sst-slp regressions.png", 3, 6, units="in", res=300)

par(mfcol=c(3,1), mar=c(0,0.5,2,0.5), oma=c(2.5,1.5,2,1.7), mgp=c(3, 0.2, 0))

z <- coef1 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("d) 1950-1968",  cex=1, side=3, adj=0)

z <- coef2 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("e) 1969-1988",  cex=1, side=3, adj=0)

z <- coef3 # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col= oce::oce.colorsPalette(64), xlab = "", ylab = "", yaxt="n", xaxt="n", zlim=zlim)
contour(x,y,z, add=T, col="grey",vfont=c("sans serif", "bold"))
polygon(box.x, box.y, border = "red")
map('world2Hires', add=T, lwd=1)
mtext("f) 1989-2008",  cex=1, side=3, adj=0)

dev.off()

