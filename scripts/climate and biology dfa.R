library(tidyverse)
library(MARSS)

# fit DFA to long-term southeast Bering climate time series 

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ok - that's a good example to compare with...now how about the Bering TS?
# load climate data 
dat <- read.csv("./data/Bering climate data.csv")

# subset - select wind and ice data - including only 1951:2013
wind <- c("year", "NW.wind.May.Sep", "NW.wind.Oct.Apr", "SE.wind.May.Sep", "SE.wind.Oct.Apr", "south.wind.stress.amj")

ice <- c("year", "ice.area.jfma","m4.march.ice", "m5.march.ice")

ice.dat <- dat %>%
  dplyr::select(all_of(ice)) %>%
  filter(year %in% 1950:2013) %>%
  dplyr::select(-year)


ice.dat <- as.matrix(t(ice.dat))

wind.dat <- dat %>%
  dplyr::select(all_of(wind)) %>%
  filter(year %in% 1950:2013) %>%
  dplyr::select(-year)

wind.dat <- as.matrix(t(wind.dat))

colnames(ice.dat) <- colnames(wind.dat) <- 1951:2013


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
  for(m in 1:2) {  # allowing up to 2 trends for wind
    dfa.model = list(A="zero", R=R, m=m)
    kemz = MARSS(wind.dat[,colnames(ice.dat) %in% 1951:2013], model=dfa.model, control=cntl.list,
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

ice.mod.data <- model.data # unconstrained produces no loadings!
wind.mod.data <- model.data # also equal var covar - 1 trend is best!

# fit best models and plot.... 
# fit model
model.list = list(A="zero", m=1, R="diagonal and unequal")
ice.mod = MARSS(ice.dat[,colnames(ice.dat) %in% 1951:2013], model=model.list, z.score=TRUE, form="dfa", control=cntl.list)
model.list = list(A="zero", m=1, R="equalvarcov")
wind.mod = MARSS(wind.dat[,colnames(wind.dat) %in% 1951:2013], model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# get CI and plot loadings...

ice.CI <- MARSSparamCIs(ice.mod)

ice.plot.CI <- data.frame(names=rownames(ice.dat),
                          mean=ice.CI$par$Z,
                          upCI=ice.CI$par.upCI$Z,
                          lowCI=ice.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)

                      
ice.plot.CI$names <- reorder(ice.plot.CI$names, ice.CI$par$Z)

ice.loadings.plot <- ggplot(ice.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

ggsave("./figs/EBS ice dfa loadings.png", width=2.5, height=4, units='in')

wind.CI <- MARSSparamCIs(wind.mod)

wind.plot.CI <- data.frame(names=rownames(wind.dat),
                           mean=wind.CI$par$Z,
                           upCI=wind.CI$par.upCI$Z,
                           lowCI=wind.CI$par.lowCI$Z)

dodge <- position_dodge(width=0.9)


wind.plot.CI$names <- reorder(wind.plot.CI$names, wind.CI$par$Z)

wind.loadings.plot <- ggplot(wind.plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[2]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle=45, hjust=1,  size=12), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

ggsave("./figs/EBS wind dfa loadings.png", width=2.5, height=4, units='in')

# plot trend
ice.trend <- data.frame(t=1951:2013,
                        estimate=as.vector(ice.mod$states),
                        conf.low=as.vector(ice.mod$states)-1.96*as.vector(ice.mod$states.se),
                        conf.high=as.vector(ice.mod$states)+1.96*as.vector(ice.mod$states.se))


ice.trend.plot <- ggplot(ice.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cb[2]) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cb[2]) + xlab("") + ylab("Trend")

ggsave("./figs/EBS ice dfa trend.png", width=4, height=2.5, units='in')

wind.trend <- data.frame(t=1951:2013,
                         estimate=as.vector(wind.mod$states),
                         conf.low=as.vector(wind.mod$states)-1.96*as.vector(wind.mod$states.se),
                         conf.high=as.vector(wind.mod$states)+1.96*as.vector(wind.mod$states.se))


wind.trend.plot <- ggplot(wind.trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cb[2]) +
  geom_hline(yintercept = 0) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cb[2]) + xlab("") + ylab("Trend")

ggsave("./figs/EBS wind dfa trend.png", width=4, height=2.5, units='in')

# and make a combined plot
png("./figs/EBS ice wind dfa plots.png", 8, 6, units='in', res=300)
ggpubr::ggarrange(ice.loadings.plot, ice.trend.plot,
                  wind.loadings.plot, wind.trend.plot,
                  ncol=2, nrow=2, labels=c("a)", "b)", "c)", "d)"),
                  widths = c(0.7, 1))
dev.off()

########################################
# now dfa of recruitment time series
# load data
dat <- read.csv("./data/EBS.recruit.time.series.csv")

# clean up
names(dat)[2:7] <- c("ylfn.age1", "trbt.age0", "fhsl.age0", "opilio.age0", "cod.age0", "poll.age0")

# lag and rename yellowfin
dat$ylfn.age0 <- c(dat$ylfn.age1[1:(nrow(dat)-1)], NA)

dat <- dat %>%
  dplyr::select(-ylfn.age1)

# plot raw ts
plot <- dat %>%
  pivot_longer(cols = -year) 

ggplot(plot, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y")
# looks about right!

# now log-transform, scale, and plot
ff <- function(x) as.vector(scale(log(x)))

dat[,2:7] <- apply(dat[,2:7], 2, ff)

plot <- dat %>%
  pivot_longer(cols = -year) 

ggplot(plot, aes(year, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~name, scales = "free_y") +
  geom_vline(xintercept = c(1968.5, 2008.5), color = "red", lty = 2)

# dashed verticals separate out the window that we're looking
# to use for first DFA: 1969-2008
# note some bogus initial values for fhsl!

# plot fhsl separately

ggplot(filter(plot, name == "fhsl.age0"),
       aes(year, value)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = c(1968.5, 2008.5), color = "red", lty = 2)

# exclude values before 1973
dat$fhsl.age0[dat$year < 1973] <- NA


# limit to 1969-2008, transpose, convert to matrix for DFA
# also removing opilio and cod b/c they do not capture a wide 
# range of sst conditions pre-1988/89

dfa.dat <- dat %>%
  dplyr::filter(year %in% 1969:2008) %>%
  dplyr::select(-year, -opilio.age0, -cod.age0) %>% 
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

# diagonal and unequal best, unconstrained 2nd best

# these best two models return super-small CIs
# fitting 3nd best model - diagonal and equal

cntl.list = list(minit=200, maxit=30000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

model.list = list(A="zero", m=1, R="diagonal and equal")
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

ggsave("./figs/EBS recruit dfa loadings 1969-2008.png", width=2.5, height=4, units='in')


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

ggsave("./figs/EBS recruit dfa trend 1969-2008.png", width=4, height=2.5, units='in')

# now fit 2nd recruitment DFA model to ____


# now save trends for other analyses

wind.save <- wind.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, wind.trend = estimate)

ice.save <- ice.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, ice.trend = estimate)

recruit.save <- recruit.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, recruit.trend = estimate)

wind.ice <- left_join(wind.save, ice.save)

wind.ice.recruit <- left_join(wind.ice, recruit.save)

write.csv(wind.ice.recruit, "./output/wind.ice.recruit.dfa.trends.csv", row.names = F)


