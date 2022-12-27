library(tidyverse)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")

theme_set(theme_bw())
# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

## Read in data --------------------------------------------
dat <- read.csv("./data/Bering climate data.csv")

dat <- dat %>%
  dplyr::select(year, SE.wind.Oct.Apr) %>%
  dplyr::rename(winter.year = year)

# temperature data
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)

str(sst)
sst$year <- as.numeric(as.character(chron::years(sst$date)))
sst$month <- as.numeric(months(sst$date))

change <- sst$year > 2030
sst$year[change] <- sst$year[change] - 100

# get Oct-Apr means (matching wind)

sst.oct.apr <- sst %>%
  dplyr::filter(month %in% c(10:12, 1:4)) %>%
  dplyr::mutate(winter.year = dplyr::if_else(month %in% 10:12, year + 1, year)) %>%
  dplyr::group_by(winter.year) %>%
  dplyr::summarise(sst.oct.apr = mean(anom)) 

dat <- left_join(dat, sst.oct.apr)

# define three eras for analysis
dat$era <- case_when(
  dat$winter.year %in% 1950:1968 ~ "1950-1968",
  dat$winter.year %in% 1969:1988 ~ "1969-1988",
  dat$winter.year %in% 1989:2008 ~ "1989-2008",
  dat$winter.year > 2008 ~ "drop"
)

# drop na for analysis
this.dat <- dat %>%
  dplyr::filter(era != "drop")


ggplot(this.dat, aes(sst.oct.apr, SE.wind.Oct.Apr, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggsave("./figs/oct-apr_sst_wind.png", width = 5, height = 4, units = 'in')

## Fit models -----------------------------------------------


## Define model formulas
wind_sst1_formula <-  bf(wind.trend ~ s(sst.ann, k = 4))

wind_sst2_formula <-  bf(wind.trend ~ sst.ann*era)

wind_sst3_formula <-  bf(wind.trend ~ s(sst.ann, k = 4, by = era))

# wind_sst2_formula <-  bf(wind.trend ~ sst.ann*era + ar(gr = era))


## Show default priors
get_prior(wind_sst1_formula, this.dat)
# get_prior(cod0_formula_nb, cod.data, family = nb)

## fit simple mode (no ar(1) term)
wind_sst1 <- brm(wind_sst1_formula,
                    data = this.dat,
                    cores = 4, chains = 4, iter = 2000,
                    save_pars = save_pars(all = TRUE),
                    control = list(adapt_delta = 0.99, max_treedepth = 10))

saveRDS(wind_sst1, file = "./output/wind_sst1.rds")


wind_sst1 <- readRDS("./output/wind_sst1.rds")

check_hmc_diagnostics(wind_sst1$fit)
neff_lowest(wind_sst1$fit)
rhat_highest(wind_sst1$fit)
summary(wind_sst1)
# bayes_R2(wind_sst1)

# y <- this.dat$wind.trend
# yrep_wind_sst1  <- fitted(wind_sst1, scale = "response", summary = FALSE)
# 
# ppc_dens_overlay(y = y, yrep = yrep_wind_sst1[sample(nrow(yrep_wind_sst1), 25), ]) +
#   ggtitle("wind_sst1")
# 
# pdf("./figs/trace_cod0_zinb_k3.pdf", width = 6, height = 4)
# trace_plot(wind_sst1)
# dev.off()

################
wind_sst2 <- brm(wind_sst2_formula,
                 data = this.dat,
                 cores = 4, chains = 4, iter = 2000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))

saveRDS(wind_sst2, file = "./output/wind_sst2.rds")


wind_sst2 <- readRDS("./output/wind_sst2.rds")

check_hmc_diagnostics(wind_sst2$fit)
neff_lowest(wind_sst2$fit)
rhat_highest(wind_sst2$fit)
summary(wind_sst2)
bayes_R2(wind_sst2)

# y <- this.dat$wind.trend
# yrep_wind_sst2  <- fitted(wind_sst2, scale = "response", summary = FALSE)
# 
# ppc_dens_overlay(y = y, yrep = yrep_wind_sst2[sample(nrow(yrep_wind_sst2), 25), ]) +
#   ggtitle("wind_sst2")
# 
# pdf("./figs/trace_cod0_zinb_k3.pdf", width = 6, height = 4)
# trace_plot(wins_sst1)
# dev.off()

####################
wind_sst3 <- brm(wind_sst3_formula,
                 data = this.dat,
                 cores = 4, chains = 4, iter = 2000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))

saveRDS(wind_sst3, file = "./output/wind_sst3.rds")


wind_sst3 <- readRDS("./output/wind_sst3.rds")

check_hmc_diagnostics(wind_sst3$fit)
neff_lowest(wind_sst3$fit)
rhat_highest(wind_sst3$fit)
summary(wind_sst3)
# bayes_R2(wind_sst3)

# y <- this.dat$wind.trend
# yrep_wind_sst3  <- fitted(wind_sst3, scale = "response", summary = FALSE)
# 
# ppc_dens_overlay(y = y, yrep = yrep_wind_sst3[sample(nrow(yrep_wind_sst3), 25), ]) +
#   ggtitle("wind_sst3")
# 
# pdf("./figs/trace_cod0_zinb_k3.pdf", width = 6, height = 4)
# trace_plot(wins_sst1)
# dev.off()


loo(wind_sst1,wind_sst2, wind_sst3)

# one observation with pareto_k > 0.7
# moment_match crashing R

################
# evaluate model with ar(1)

wind_sst4_formula <-  bf(wind.trend ~ sst.ann*era + ar(gr = era))

wind_sst4 <- brm(wind_sst4_formula,
                 data = this.dat,
                 cores = 4, chains = 4, iter = 2000,
                 save_pars = save_pars(all = TRUE),
                 control = list(adapt_delta = 0.99, max_treedepth = 10))

saveRDS(wind_sst4, file = "./output/wind_sst4.rds")

wind_sst4 <- readRDS("./output/wind_sst4.rds")

check_hmc_diagnostics(wind_sst4$fit)
neff_lowest(wind_sst4$fit)
rhat_highest(wind_sst4$fit)
summary(wind_sst4)
bayes_R2(wind_sst4)

# y <- this.dat$wind.trend
# yrep_wind_sst4  <- fitted(wind_sst4, scale = "response", summary = FALSE)
# 
# ppc_dens_overlay(y = y, yrep = yrep_wind_sst4[sample(nrow(yrep_wind_sst4), 25), ]) +
#   ggtitle("wind_sst4")
# 
# pdf("./figs/trace_cod0_zinb_k3.pdf", width = 6, height = 4)
# trace_plot(wins_sst1)
# dev.off()

loo(wind_sst2, wind_sst4)

# wind_sst4 much better - need to fix pareto_k

ce1s_1 <- conditional_effects(wind_sst4, effect = "sst.ann:era", re_formula = NA,
                              probs = c(0.025, 0.975), points = T)  

################################################
# now sst-biology trend

# temperature data
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)

str(sst)
sst$year <- as.numeric(as.character(chron::years(sst$date)))
sst$month <- as.numeric(months(sst$date))

change <- sst$year > 2030
sst$year[change] <- sst$year[change] - 100

# get Oct-Apr means (matching wind)

sst.oct.apr <- sst %>%
  dplyr::filter(month %in% c(10:12, 1:4)) %>%
  dplyr::mutate(winter.year = dplyr::if_else(month %in% 10:12, year + 1, year)) %>%
  dplyr::group_by(winter.year) %>%
  dplyr::summarise(sst.oct.apr = mean(anom)) 


sst.oct.apr$sst.oct.apr <- zoo::rollmean(sst.oct.apr$sst.oct.apr, 3, fill=NA)# three-year rolling mean

trends <- read.csv("./output/wind.ice.recruit.dfa.trends.csv")

trends <- trends %>%
  dplyr::select(year, recruit.trend) %>%
  dplyr::rename(winter.year = year)

this.dat <- dplyr::left_join(trends, sst.oct.apr) %>%
  dplyr::filter(winter.year %in% 1969:2008) %>%
  dplyr::mutate(era = dplyr::if_else(winter.year %in% 1969:1988, "1969-1988", "1989-2008"))

ggplot(this.dat, aes(winter.year, recruit.trend)) +
  geom_line() + 
  geom_point()

ggplot(this.dat, aes(sst.oct.apr, recruit.trend, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se = F)

ggsave("./figs/biology-sst by era.png", width = 6, height = 4, units = 'in')

# fit linear model 
linear_mod_ebs <- nlme::gls(recruit.trend ~ sst.oct.apr*era, correlation = corAR1(),
                            data = this.dat)


summary(linear_mod_ebs)

## add GOA analysis -------------------------------------

g_dat <- read.csv("./data/GOA community data.csv")
names(g_dat)[1] <- "year"

# limit to 1969-2008, transpose, convert to matrix for DFA
# also removing opilio and cod b/c they do not capture a wide 
# range of sst conditions pre-1988/89

dfa.dat <- g_dat %>%
  dplyr::filter(year %in% 1969:2008) %>%
  dplyr::select(-year) %>% 
  t()

colnames(dfa.dat) <- 1969:2008
dfa.dat <- as.matrix(dfa.dat)

# find best error structure for 1-trend model

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# set up forms of R matrices

levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")
model.data = data.frame()

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS::MARSS(dfa.dat, model=dfa.model, control=cntl.list,
                        form="dfa", z.score=TRUE)
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare

model.data$dAICc <- model.data$AICc-min(model.data$AICc)

model.data <- model.data %>%
  arrange(dAICc)

model.data

# diagonal and unequal best, equalvarcov 2nd best

# these best two models return super-small CIs
# fitting 3nd best model - diagonal and equal

cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=1, R="diagonal and unequal")
recruit.mod = MARSS::MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

recruit.CI <- MARSSparamCIs(recruit.mod)

recruit.plot.CI <- data.frame(names=rownames(dfa.dat),
                              mean=recruit.CI$par$Z,
                              upCI=recruit.CI$par.upCI$Z,
                              lowCI=recruit.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)


recruit.plot.CI$names <- reorder(recruit.plot.CI$names, recruit.CI$par$Z)

recruit.loadings.plot <- ggplot(recruit.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

recruit.loadings.plot

ggsave("./figs/GOA recruit dfa loadings 1969-2008.png", width=2.5, height=4, units='in')

# plot trend
recruit.trend <- data.frame(t=1969:2008,
                            estimate=as.vector(recruit.mod$states),
                            conf.low=as.vector(recruit.mod$states)-1.96*as.vector(recruit.mod$states.se),
                            conf.high=as.vector(recruit.mod$states)+1.96*as.vector(recruit.mod$states.se))


recruit.trend.plot <- ggplot(recruit.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cb[2]) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cb[2]) + 
  theme(axis.title.x = element_blank()) +
  ylab("Trend")

recruit.trend.plot

ggsave("./figs/GOA recruit dfa trend 1969-2008.png", width=4, height=2.5, units='in')

# combine plots and save

png("./figs/GOA_recr_dfa_loadings_plot.png", width = 9, height = 3, units = 'in', res = 300)

ggpubr::ggarrange(recruit.loadings.plot, recruit.trend.plot, labels = "auto", nrow = 1)

dev.off()

# load GOA SST data
goa_clim_dat <- read.csv("./data/GOA environmental data.csv")

# combine GOA recruitment trend and SST with EBS data to plot

ebs_dat <- this.dat %>%
  dplyr::rename(sst = sst.oct.apr) %>%
  dplyr::mutate(system = "Eastern Bering Sea")

goa_dat <- recruit.trend %>%
  dplyr::rename(winter.year = t,
         recr.trend = estimate) %>%
  dplyr::select(-conf.high, -conf.low) %>%
  dplyr::mutate(system = "Gulf of Alaska")


goa_clim_dat <- goa_clim_dat %>%
  dplyr::rename(winter.year = X,
                sst = SST) %>%
  dplyr::select(winter.year, sst)

goa_dat <- left_join(goa_dat, goa_clim_dat)

goa_dat <- goa_dat %>%
  dplyr::mutate(era = if_else(winter.year <= 1988, "1969-1988", "1989-2008")) %>%
  dplyr::rename(recruit.trend = recr.trend) %>%
  dplyr::select(winter.year, recruit.trend, sst, era, system)

# reverse trend to match EBS!
goa_dat$recruit.trend = -goa_dat$recruit.trend

ebs_goa_dat <- rbind(ebs_dat, goa_dat)

ggplot(ebs_goa_dat, aes(sst, recruit.trend, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se=F) +
  facet_wrap(~system, scales = "free", ncol = 1) +
  scale_color_manual(values = cb[c(2,6)]) +
  labs(x = "Winter mean SST (°C)",
       y = "Recruitment trend")


ggsave("./figs/EBS_GOA_recruitment_vs_sst.png", width = 6.5, height = 6, units = 'in')

# fit linear model 
linear_mod_goa <- nlme::gls(recruit.trend ~ sst*era, correlation = corAR1(),
                            data = dplyr::filter(ebs_goa_dat, system == "Gulf of Alaska"))


summary(linear_mod_goa)


################################################
# examine SLP EOF1 as the driver


library(tidyverse)
library(sde)
library(FactoMineR)

# load slp and cell weights
slp <- read.csv("./data/north.pacific.slp.anom.csv", row.names = 1)

weights <- read.csv("./data/north.pacific.slp.weights.csv", row.names = 1)

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

#############################

## changing slp pc1 - wind response

#############################
  
  ###########################
  # first window - response = 1951:1968
  temp.slp <- slp[which(rownames(slp) == "10/01/50") : which(rownames(slp) == "04/01/68"),] 
  
  # and fit EOF
  pca.temp <- FactoMineR::svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
  pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]
  
  # scale
  pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!
 
  
  # make a data frame for selecting / correlating
  temp.dat <- data.frame(year = as.numeric(as.character(chron::years(row.names(temp.slp)))),
                         month = as.numeric(months((row.names(temp.slp)))),
                         pc1.temp = pc1.temp)
  
  change <- temp.dat$year > 2030
  temp.dat$year[change] <- temp.dat$year[change] - 100
  
  
  temp.dat$winter.year <- case_when(
    temp.dat$month %in% 9:12 ~ temp.dat$year+1, 
    temp.dat$month %in% 1:3 ~ temp.dat$year,
    temp.dat$month %in% 4:8 ~ 999)
  
  
  temp.dat <- temp.dat %>%
    dplyr::filter(winter.year != 999) %>%
    dplyr::group_by(winter.year) %>%
    dplyr::summarise(pc1 = mean(pc1.temp))
  
  
  # and add wind variable
  clim.dat <- read.csv("./data/Bering climate data.csv")
  
  # subset - select wind 
  wind <- clim.dat %>%
    dplyr::select(year, SE.wind.Oct.Apr) %>%
    dplyr::rename(winter.year = year)
  
  
  slp.wind.dat1 <- left_join(temp.dat, wind)
  slp.wind.dat1$era <- "1951-1968"
  
  ggplot(slp.wind.dat1, aes(pc1, SE.wind.Oct.Apr)) +
    geom_point()
  
  
  # second window - response = 1969:1988
  temp.slp <- slp[which(rownames(slp) == "10/01/68") : which(rownames(slp) == "04/01/88"),] 
  
  # and fit EOF
  pca.temp <- FactoMineR::svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
  pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]
  
  # scale
  pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!
  
  
  # make a data frame for selecting / correlating
  temp.dat <- data.frame(year = as.numeric(as.character(chron::years(row.names(temp.slp)))),
                         month = as.numeric(months((row.names(temp.slp)))),
                         pc1.temp = pc1.temp)
  
  change <- temp.dat$year > 2030
  temp.dat$year[change] <- temp.dat$year[change] - 100
  
  # we want smoothed values one month earlier than wind seasonal window
  # so this would be Sept-March pc1 for Oct-Apr wind
  
  temp.dat$winter.year <- case_when(
    temp.dat$month %in% 9:12 ~ temp.dat$year+1, 
    temp.dat$month %in% 1:3 ~ temp.dat$year,
    temp.dat$month %in% 4:8 ~ 999)
  
  
  temp.dat <- temp.dat %>%
    dplyr::filter(winter.year != 999) %>%
    dplyr::group_by(winter.year) %>%
    dplyr::summarise(pc1 = mean(pc1.temp))
  

  slp.wind.dat2 <- left_join(temp.dat, wind)
  slp.wind.dat2$era <- "1969-1988"
  
  ggplot(slp.wind.dat2, aes(pc1, SE.wind.Oct.Apr)) +
    geom_point()
  
  ## third era
  temp.slp <- slp[which(rownames(slp) == "10/01/88") : which(rownames(slp) == "04/01/08"),] 
  
  # and fit EOF
  pca.temp <- FactoMineR::svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
  pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]
  
  # scale
  pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!
  
  
  # make a data frame for selecting / correlating
  temp.dat <- data.frame(year = as.numeric(as.character(chron::years(row.names(temp.slp)))),
                         month = as.numeric(months((row.names(temp.slp)))),
                         pc1.temp = pc1.temp)
  
  change <- temp.dat$year > 2030
  temp.dat$year[change] <- temp.dat$year[change] - 100
  
  # we want smoothed values one month earlier than wind seasonal window
  # so this would be Sept-March pc1 for Oct-Apr wind
  
  temp.dat$winter.year <- case_when(
    temp.dat$month %in% 9:12 ~ temp.dat$year+1, 
    temp.dat$month %in% 1:3 ~ temp.dat$year,
    temp.dat$month %in% 4:8 ~ 999)
  
  
  temp.dat <- temp.dat %>%
    dplyr::filter(winter.year != 999) %>%
    dplyr::group_by(winter.year) %>%
    dplyr::summarise(pc1 = mean(pc1.temp))
  
  
  slp.wind.dat3 <- left_join(temp.dat, wind)
  slp.wind.dat3$era <- "1989-2008"
  
  ggplot(slp.wind.dat3, aes(pc1, SE.wind.Oct.Apr)) +
    geom_point()
 
  plot.dat <- rbind(slp.wind.dat1, slp.wind.dat2, slp.wind.dat3)
  
  ggplot(plot.dat, aes(pc1, SE.wind.Oct.Apr, color = era)) +
    geom_point() +
    geom_smooth(method = "lm", se = F)
  
  
  ggsave("./figs/slp pc1-wind by era.png", width = 6, height = 4, units = 'in')
  
  ##########################################################
  
  
  
  
##########################################################
  
  # repeat above with ice as response
  
##########################################################
##########################################################
  # first window - response = 1951:1968
  
  # for ice, we look at slp for Nov-Mar
  temp.slp <- slp[which(rownames(slp) == "11/01/50") : which(rownames(slp) == "03/01/68"),] 
  
  # and fit EOF
  pca.temp <- FactoMineR::svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
  pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]
  
  # scale
  pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!

  # make a data frame for selecting / correlating
  temp.dat <- data.frame(year = as.numeric(as.character(chron::years(row.names(temp.slp)))),
                         month = as.numeric(months((row.names(temp.slp)))),
                         pc1.temp = pc1.temp)
  
  change <- temp.dat$year > 2030
  temp.dat$year[change] <- temp.dat$year[change] - 100
  
  # we want smoothed values one month earlier than ice seasonal iceow
  # so this would be Sept-March pc1 for Oct-Apr ice
  
  temp.dat$winter.year <- case_when(
    temp.dat$month %in% 11:12 ~ temp.dat$year+1, 
    temp.dat$month %in% 1:3 ~ temp.dat$year,
    temp.dat$month %in% 4:10 ~ 999)
  
  
  temp.dat <- temp.dat %>%
    dplyr::filter(winter.year != 999) %>%
    dplyr::group_by(winter.year) %>%
    dplyr::summarise(pc1 = mean(pc1.temp))
  
  
  # and add ice variable
  dfa.dat <- read.csv("./output/wind.ice.recruit.dfa.trends.csv")
  
  # subset - select ice 
  ice <- dfa.dat %>%
    dplyr::select(year, ice.trend) %>%
    dplyr::rename(winter.year = year)
  
  
  slp.ice.dat1 <- left_join(temp.dat, ice)
  slp.ice.dat1$era <- "1951-1968"
  
  ggplot(slp.ice.dat1, aes(pc1, ice.trend)) +
    geom_point()
  
  ## second era
  
  temp.slp <- slp[which(rownames(slp) == "11/01/68") : which(rownames(slp) == "03/01/88"),] 
  
  # and fit EOF
  pca.temp <- FactoMineR::svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
  pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]
  
  # scale
  pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!
  
  # make a data frame for selecting / correlating
  temp.dat <- data.frame(year = as.numeric(as.character(chron::years(row.names(temp.slp)))),
                         month = as.numeric(months((row.names(temp.slp)))),
                         pc1.temp = pc1.temp)
  
  change <- temp.dat$year > 2030
  temp.dat$year[change] <- temp.dat$year[change] - 100
  

  temp.dat$winter.year <- case_when(
    temp.dat$month %in% 11:12 ~ temp.dat$year+1, 
    temp.dat$month %in% 1:3 ~ temp.dat$year,
    temp.dat$month %in% 4:10 ~ 999)
  
  
  temp.dat <- temp.dat %>%
    dplyr::filter(winter.year != 999) %>%
    dplyr::group_by(winter.year) %>%
    dplyr::summarise(pc1 = mean(pc1.temp))
  
  slp.ice.dat2 <- left_join(temp.dat, ice)
  slp.ice.dat2$era <- "1969-1988"
  
  ggplot(slp.ice.dat2, aes(pc1, ice.trend)) +
    geom_point()
  
  ## third era
  temp.slp <- slp[which(rownames(slp) == "11/01/88") : which(rownames(slp) == "03/01/08"),] 
  
  # and fit EOF
  pca.temp <- FactoMineR::svd.triplet(cov(temp.slp), col.w=weights[,1]) #weighting the columns
  pc1.temp <- as.matrix(temp.slp) %*% pca.temp$U[,1]
  
  # scale
  pc1.temp <- as.vector(scale(-pc1.temp)) # reversing sign!
  
  # make a data frame for selecting / correlating
  temp.dat <- data.frame(year = as.numeric(as.character(chron::years(row.names(temp.slp)))),
                         month = as.numeric(months((row.names(temp.slp)))),
                         pc1.temp = pc1.temp)
  
  change <- temp.dat$year > 2030
  temp.dat$year[change] <- temp.dat$year[change] - 100
  
  # we want smoothed values one month earlier than ice seasonal iceow
  # so this would be Sept-March pc1 for Oct-Apr ice
  
  temp.dat$winter.year <- case_when(
    temp.dat$month %in% 11:12 ~ temp.dat$year+1, 
    temp.dat$month %in% 1:3 ~ temp.dat$year,
    temp.dat$month %in% 4:10 ~ 999)
  
  
  temp.dat <- temp.dat %>%
    dplyr::filter(winter.year != 999) %>%
    dplyr::group_by(winter.year) %>%
    dplyr::summarise(pc1 = mean(pc1.temp))
  
  slp.ice.dat3 <- left_join(temp.dat, ice)
  slp.ice.dat3$era <- "1989-2008"
  
  ggplot(slp.ice.dat3, aes(pc1, ice.trend)) +
    geom_point()
  
  plot <- rbind(slp.ice.dat1, slp.ice.dat2, slp.ice.dat3)
  
  ggplot(plot, aes(pc1, ice.trend, color = era)) +
    geom_point()
  
  ggsave("./figs/ice vs slp eof1.png", width = 6, height = 4, units = 'in')
  
  #####################
  
  # fix bad year values
  change <- temp.dat$year > 2030
  temp.dat$year[change] <-temp.dat$year[change] - 100

  temp.dat


cor.out <- cor.out %>%
  pivot_longer(cols = -end.date)

ggplot(cor.out, aes(end.date, value)) +
  geom_line() +
  facet_wrap(~name, scales = "free_y", ncol = 1)