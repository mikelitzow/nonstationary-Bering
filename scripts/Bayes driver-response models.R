library(ggplot2)
library(dplyr)
library(plyr)
library(mgcv)
library(rstan)
library(brms)
library(bayesplot)
source("./scripts/stan_utils.R")


## Read in data --------------------------------------------
dat <- read.csv("./output/wind.ice.recruit.dfa.trends.csv")

# temperature data
sst <- read.csv("./data/ebs.sst.anom.csv", row.names = 1)

str(sst)
sst$year <- as.numeric(as.character(chron::years(sst$date)))

# get annual means (for wind)

sst.ann <- data.frame(year = 1950:2013,
                    sst.ann = tapply(sst$anom, sst$year, mean))

dat <- left_join(dat, sst.ann)

# define three eras for analysis
dat$era <- case_when(
  dat$year %in% 1950:1968 ~ "1958-1968",
  dat$year %in% 1968:1988 ~ "1969-1988",
  dat$year %in% 1989:2008 ~ "1989-2008",
  dat$year > 2008 ~ "drop"
)


## Fit models -----------------------------------------------


## Define model formulas
wind_sst1_formula <-  bf(wind.trend ~ sst.ann*era)

wind_sst2_formula <-  bf(wind.trend ~ sst.ann*era + ar(gr = era))

# drop na for analysis
this.dat <- dat %>%
  dplyr::select(wind.trend, sst.ann, era) %>%
  dplyr::filter(era != "drop")

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


check_hmc_diagnostics(wind_sst1)
neff_lowest(wind_sst1)
rhat_highest(wind_sst1)
summary(wind_sst1)
bayes_R2(wind_sst1)
# plot(cod0_zinb_k3$criteria$loo, "k")
plot(conditional_smooths(cod0_zinb_k3), ask = FALSE)
y <- cod.data$cod
yrep_cod0_zinb_k3  <- fitted(cod0_zinb_k3, scale = "response", summary = FALSE)
ppc_dens_overlay(y = y, yrep = yrep_cod0_zinb_k3[sample(nrow(yrep_cod0_zinb_k3), 25), ]) +
  xlim(0, 500) +
  ggtitle("cod0_zinb_k3")
pdf("./figs/trace_cod0_zinb_k3.pdf", width = 6, height = 4)
trace_plot(wins_sst1)
dev.off()
