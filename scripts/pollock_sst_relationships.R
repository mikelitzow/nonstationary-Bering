library(tidyverse)

# set colors
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# load EBS data 
ebs_dat <- read.csv("./data/EBS.recruit.time.series.csv") %>%
  dplyr::select(year, pollock.age0.R) %>%
  dplyr::rename(EBS_poll = pollock.age0.R) %>%
  dplyr::filter(year %in% 1968:2008)

ebs_sst <- read.csv("./data/Bering climate data.csv") %>%
  dplyr::select(year, south.sst.ndjfm) %>%
  dplyr::mutate(sst.3 = zoo::rollmean(south.sst.ndjfm, 3, align = "right", fill = NA))

ebs_dat <- left_join(ebs_dat, ebs_sst) %>%
  dplyr::rename(poll_R0 = EBS_poll,
                sst = south.sst.ndjfm) %>%
  dplyr::mutate(system = "Eastern Bering Sea")

# load GOA data
goa_clim_dat <- read.csv("./data/GOA environmental data.csv") %>%
  dplyr::rename(year = X,
                sst = SST) %>%
  dplyr::select(year, sst) %>%
  dplyr::mutate(sst.3 = zoo::rollmean(sst, 3, align = "right", fill = NA))

goa_dat <- read.csv("./data/GOA community data.csv") %>%
  dplyr::rename(year = X,
                poll_R0 = Walleye.pollock) %>% 
  dplyr::select(year, poll_R0) %>%
  dplyr::mutate(system = "Gulf of Alaska")
  
goa_dat <- left_join(goa_dat, goa_clim_dat) %>%
  dplyr::select(year, poll_R0, sst, sst.3, system)

ebs_goa_poll <- rbind(ebs_dat, goa_dat) %>%
  dplyr::filter(year %in% 1969:2008) %>%
  dplyr::mutate(era = if_else(year %in% 1969:1988, "1969-1988", "1989-2008"))
  
ggplot(ebs_goa_poll, aes(sst, poll_R0, color = era)) +
  geom_point() +
  geom_smooth(method = "lm", se=F) +
  facet_wrap(~system, scales = "free", ncol = 1) +
  scale_color_manual(values = cb[c(2,6)]) +
  labs(x = "Winter mean SST (°C)",
       y = "Recruitment trend")


ggsave("./figs/EBS_GOA_recruitment_vs_sst.png", width = 6.5, height = 6, units = 'in')

linear_mod_poll_ebs <- nlme::gls(poll_R0 ~ sst*era, correlation = corAR1(),
                            data = dplyr::filter(ebs_goa_poll, system == "Eastern Bering Sea"))


summary(linear_mod_poll_ebs)

linear_mod_poll_goa <- nlme::gls(poll_R0 ~ sst*era, correlation = corAR1(),
                                 data = dplyr::filter(ebs_goa_poll, system == "Gulf of Alaska"))


summary(linear_mod_poll_goa)
