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

# now save trends for other analyses

wind.save <- wind.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, wind.trend = estimate)

ice.save <- ice.trend %>%
  dplyr::select(t, estimate) %>%
  rename(year = t, ice.trend = estimate)

wind.ice <- left_join(wind.save, ice.save)

write.csv(wind.ice, "./data/wind.ice.dfa.trends.csv", row.names = F)