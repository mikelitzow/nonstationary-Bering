library(dplyr)
library(ncdf4)
library(maps)
library(maptools)
library(mapdata)
library(fields)
library(oce)
library(chron)
library(zoo)
library(mapproj)
library(sinkr)


# to load sinkr, use package devtools and uncomment the following:
 # devtools::install_github("marchtaylor/sinkr")

# load slp
dat <- nc_open("data/NCEP.NCAR.slp.nc")

x <- ncvar_get(dat, "longitude")
y <- ncvar_get(dat, "latitude")
slp <- ncvar_get(dat, "slp", verbose = F)
dim(slp) # 144 long, 29 lat, 864 months

# need to reverse latitude for plotting!
y <- rev(y)
slp <- slp[,29:1,]

# first, extract dates
raw <- ncvar_get(dat, "time") # seconds since 1-1-1970
h <- raw/(24*60*60)
d <- dates(h, origin = c(1,1,1970))
m <- months(d)
yr <- as.numeric(as.character(years(d)))

# make vectors of lat/long and add (with date) as dimnames
lat <- rep(y, length(x))   
lon <- rep(x, each = length(y)) 


# Change data into a matrix with months / cells for rows / columns
slp <- aperm(slp, 3:1)  
slp <- matrix(slp, nrow=dim(slp)[1], ncol=prod(dim(slp)[2:3]))  

dimnames(slp) <- list(as.character(d), paste("N", lat, "E", lon, sep=""))

# plot to confirm that everything is ok
z <- colMeans(slp, na.rm=T) # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "")

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', add=T, lwd=1)

# limit to 1951:2013
slp <- slp[yr %in% 1950:2013,]
m <- m[yr %in% 1950:2013]
yr <- yr[yr %in% 1950:2013]


# let's get a regional PCA for the area affecting the Bering more or less directly

poly.x <- c(150, 150, 230, 230, 150) 
poly.y <- c(38, 64, 64, 38, 38)

ebs.slp <- slp

xp <- cbind(poly.x, poly.y)
loc=cbind(lon, lat)
check <- in.poly(loc, xp=xp)

ebs.slp[,!check] <- NA

# plot to check

z <- colMeans(ebs.slp, na.rm=T) # mean value for each cell
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image(x,y,z, col=tim.colors(64), xlab = "", ylab = "", xlim=c(100,260), ylim=c(30,70))

contour(x,y,z, add=T, col="white",vfont=c("sans serif", "bold"))
map('world2Hires', add=T, lwd=1)

# identify columns containing NA
out <- is.na(colMeans(ebs.slp)) 

X <- ebs.slp[,!out] 

# get anomalies

f <- function(x) tapply(x, m, mean)  # function to compute monthly means for a single time series
mu <- apply(X, 2, f)	# compute monthly means for each time series (cell)
mu <- mu[rep(1:12, length(m)/12),]  # replicate means matrix for each year at each location

slp.anom <- X - mu   # compute matrix of anomalies

# will use winter (NDJFM) SLP 

# define the winter period
win <- c("Nov", "Dec", "Jan", "Feb", "Mar")

# set Nov and Dec equal to the year corresponding to Jan!
win.yr <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr) 

# and restrict both winter year and slp to the winter months
win.yr <- win.yr[m %in% win]
win.anom <- slp.anom[m %in% win,]

# and limit win.anom to the complete winters (1951-2013)
win.anom <- win.anom[win.yr %in% 1951:2013,]
win.yr <- win.yr[win.yr %in% 1951:2013]

# get a vector of weights (square root of the cosine of latitude)
lat.weights <- lat
weight <- sqrt(cos(lat.weights*pi/180))

# EOF by era
# (doing this without weights for now as I don't have FactoMineR loaded)
EOF.early <- svd(cov(win.anom[win.yr <= 1988,])) #weighting the columns
EOF.late <- svd(cov(win.anom[win.yr > 1988,]))

# get loadings for PC1/2 by era
eig.1.early <- EOF.early$u[,1]
eig.2.early <- EOF.early$u[,2]

eig.1.late <- EOF.late$u[,1]
eig.2.late <- EOF.late$u[,2]


# proportion of variance in each era
ev1 <- EOF.early$d
ev2 <- EOF.late$d

ev1[1:5]/sum(ev1)*100
ev2[1:5]/sum(ev2)*100
# so PC1 is the only important axis, and variance explained is very similar between eras!

# and plot
my.col <- oce.colorsPalette(64)
png("figs/exploratory era EOF SLP.png", 1.2*11.4/2.54, 11.4/2.54, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, mfrow=c(2,2), cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.2))

# set the limit for plotting 
lim <- range(eig.1.early, eig.1.late, eig.2.early, eig.2.late, na.rm=T)

z <- rep(NA, ncol(slp))
z[!out] <- eig.1.early
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", xlim=c(130,250), ylim=c(30,68),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', add=T, lwd=1)

mtext("a", adj=0.05, line=-1.3, cex=mt.cex)
mtext("EOF1 1951-1988", cex=0.8)

z <- rep(NA, ncol(slp))
z[!out] <- eig.1.late
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", xlim=c(130,250), ylim=c(30,68),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', add=T, lwd=1)

mtext("b", adj=0.05, line=-1.3, cex=mt.cex)
mtext("EOF1 1989-2013", cex=0.8)

z <- rep(NA, ncol(slp))
z[!out] <- eig.2.early
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", xlim=c(130,250), ylim=c(30,68),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', add=T, lwd=1)

mtext("c", adj=0.05, line=-1.3, cex=mt.cex)
mtext("EOF2 1951-1988", cex=0.8)

z <- rep(NA, ncol(slp))
z[!out] <- eig.2.late
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=my.col, zlim=c(-lim[2], lim[2]), xlab = "", ylab = "", yaxt="n", xaxt="n", xlim=c(130,250), ylim=c(30,68),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', add=T, lwd=1)

mtext("d", adj=0.05, line=-1.3, cex=mt.cex)
mtext("EOF2 1989-2013", cex=0.8)

dev.off()

# so! best to fit a single EOF to the whole time series as there's no meaningful difference in loadings!

EOF.all <- svd(cov(win.anom))

# now! get PCs
pc1 <- win.anom %*% EOF.all$u[,1]


# ad annual winter means
pc1 <- tapply(pc1, win.yr, mean)

plot <- data.frame(year=1951:2013,
                   pc1=scale(pc1))


library(ggplot2)
ggplot(plot, aes(year, pc1)) +
  theme_bw() +
  geom_line() +
  geom_vline(xintercept = 1988.5, lty=2)

# and rolling SD for PC1 SD
# scale pc1!
pc1 <- scale(pc1)
names(pc1) <- 1951:2013

pc1.sd <- data.frame()

for(i in 1958:2006){

  temp <- pc1[names(pc1) %in% (i-7):(i+7)]
  pc1.sd <- rbind(pc1.sd,
                  data.frame(year=i,
                             pc1.sd=sd(temp)))
  
}



cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

ggplot(pc1.sd, aes(year, pc1.sd)) +
  theme_bw() +
  geom_line() +
  geom_vline(xintercept = 1988.5, lty=2)

# save a joint plot of PC1 and rolling SD!
temp <- data.frame(year=1951:2013,
                   pc1=pc1)

joint.plot <- left_join(temp, pc1.sd) %>%
  tidyr::pivot_longer(cols=-year)


ggplot(joint.plot, aes(year, value)) +
  theme_bw() +
  geom_line() +
  facet_wrap(~name, scales="free_y", nrow=2) +
  geom_vline(xintercept = 1988.5, lty=2)

ggsave("figs/N. Pac slp PC1 time series and SD.png", width=3, height=5, units='in')

# and one-panel plot of EOF.all loadings
eig.1.all <- EOF.all$u[,1]
my.col <- oce.colorsPalette(64)
png("figs/EOF1 lodaings SLP 1951-2013.png", 5, 5, units="in", res=300)

# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 0.2
tc.l <- -0.2

par(mar=c(0.5,0.5,1.5,1),  tcl=tc.l, mgp=c(1.5,0.3,0), las=1, cex.axis=0.8, cex.lab=0.8, oma=c(0,0,0,0.5))

# set the limit for plotting 
lim <- range(eig.1.all, na.rm=T)

z <- rep(NA, ncol(slp))
z[!out] <- eig.1.early
z <- t(matrix(z, length(y))) 

image.plot(x,y,z, col=my.col, zlim=c(-0.1, 0.1), xlab = "", ylab = "", yaxt="n", xaxt="n", xlim=c(130,250), ylim=c(30,68),
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x, y, z, add=T, drawlabels = T, lwd=0.7, col="grey") 
map('world2Hires', add=T, lwd=1)

dev.off()

#  quickly, load sst!
dat2 <- read.csv("data/climate data.csv")

dat2 <- dat2 %>%
  select(year, south.sst.ndjfm)
names(dat2)[2] <- "EBS.sst"

plot <- left_join(temp, dat2) 

plot$era <- as.factor(ifelse(plot$year <= 1988, 1, 2))

ggplot(plot, aes(EBS.sst, pc1, color=era)) + 
  theme_bw() +
  geom_point() +
  geom_smooth(method="lm", se=F)

# so sst associated with higher pc in earlier era?? 
# could check with the Bayesian regression....
library(rstan)
library(plyr)
library(rstanarm)
library(bayesplot)

stan.new <- data.frame(year=1951:2013,
                       EBS.sst=plot$EBS.sst,
                       name="slp.pc1",
                       value=plot$pc1)

stan.new$era <- ifelse(stan.new$year <= 1988, "1951-1988", "1989-2013")


# and run regression
slp.sst.stan <- stan_glm(value ~ era + EBS.sst + EBS.sst:era,
                     data = stan.new,
                     chains = 4, cores = 4, thin = 1,
                     warmup = 1000, iter = 4000, refresh = 0,
                     prior = normal(location = 0, scale = 5, autoscale = FALSE),
                     prior_intercept = normal(location = 0, scale = 5, autoscale = FALSE),
                     prior_aux = student_t(df = 3, location = 0, scale = 5, autoscale = FALSE))

##############
# and plot!

lst <- list(slp.sst.stan)


# extract intercepts
lst.int <- lapply(lst, function(x) {
  beta <- as.matrix(x, pars = c("(Intercept)", "era1989-2013"))
  data.frame(key = unique(x$data$name),
             era1 = beta[ , 1],
             era2 = beta[ , 1] + beta[ , 2])
})

coef_indv_arm <- plyr::rbind.fill(lst.int)

mdf_indv_arm <- reshape2::melt(coef_indv_arm, id.vars = "key")

## extract slopes
lst.slope <- lapply(lst, function(x) {
  beta <- as.matrix(x, pars = c("EBS.sst", "era1989-2013:EBS.sst"))
  data.frame(key = unique(x$data$name),
             era1 = beta[ , 1],
             era2 = beta[ , 1] + beta[ , 2])
})
coef_slope <- plyr::rbind.fill(lst.slope)
mdf_slope <- reshape2::melt(coef_slope, id.vars = "key")


# plot intercepts
int <- ggplot(mdf_indv_arm, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3]), labels=c("1950-1988", "1989-2013")) +
  theme(legend.title = element_blank(), legend.position = 'top') +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Intercept (scaled anomaly)",
       y = "Posterior density") +
  facet_wrap( ~ key, scales="free")
print(int)


# plot slopes
slope <- ggplot(mdf_slope, aes(x = value, fill = variable)) +
  theme_bw() +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c(cb[2], cb[3], cb[4]), labels=c("1951-1988", "1989-2013")) +
  theme(legend.title = element_blank(), legend.position = 'top') +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Slope (scaled anomaly)",
       y = "Posterior density")
print(slope)

# not much difference at all!

ggsave("figs/era slopes - SLP PC1 vs EBS sst.png", width=3, height=3, units='in')

############

# will use winter (NDJFM) SLP 

# define the winter period
win <- c("Nov", "Dec", "Jan", "Feb", "Mar")

# set Nov and Dec equal to the year corresponding to Jan!
win.yr <- ifelse(m %in% c("Nov", "Dec"), yr+1, yr) 
# and restrict both winter year and slp to the winter months
win.yr <- win.yr[m %in% win]
win.slp <- ebs.slp[m %in% win,]


# now get annual winter means for each cell
ff <- function(x) tapply(x, win.yr, mean)

win.slp <- apply(win.slp, 2, ff)

# limit to 
win.slp <- win.slp[rownames(win.slp) < 2020,]
ebs.win.slp <- rowMeans(win.slp, na.rm = T)


# quick plot!
plot(ebs.win.slp[names(ebs.win.slp) %in% 1951:2019], sst$EBS.sst[sst$year %in% 1951:2019], type="n")
text(ebs.win.slp[names(ebs.win.slp) %in% 1951:1988], sst$EBS.sst[sst$year %in% 1951:1988], labels=1951:1988, col="blue")
text(ebs.win.slp[names(ebs.win.slp) %in% 1989:2012], sst$EBS.sst[sst$year %in% 1989:2012], labels=1989:2012, col="red")


plot.dat <- data.frame(year=1951:2019,
                       slp=ebs.win.slp[names(ebs.win.slp) %in% 1951:2019],
                       sst=sst$EBS.sst[sst$year %in% 1951:2019])

plot.dat$era <- as.factor(ifelse(plot.dat$year < 1989, 1, 
                       ifelse(plot.dat$year %in% 1989:2013, 2, 3)))


ggplot(filter(plot.dat, era != 3), aes(slp, sst, color=era)) +
  theme_bw() +
  geom_point() +
  geom_smooth(method="lm")


plot.dat$slp <- scale(plot.dat$slp)
new.dat <- plot.dat %>%
  pivot_longer(cols=c(-year, -era))

ggplot(new.dat, aes(year, value, color=name)) +
  theme_bw() +
  geom_line() + 
  geom_vline(xintercept = 1988.5, lty=2)

ggplot(new.dat, aes(value, fill=era)) +
  theme_bw() +
  geom_density(alpha=0.2) +
  facet_wrap(~name) +
  xlim(-3,3)

# separate SLP, PDO, and NPGO into era-specific chunks for regression maps
win.slp1 <- win.slp[rownames(win.slp) %in% 1951:1988,]
win.slp2 <- win.slp[rownames(win.slp) %in% 1989:2013,]

GOA1 <- sst$GOA.sst[sst$year %in% 1951:1988]
GOA2 <- sst$GOA.sst[sst$year %in% 1989:2013]

EBS1 <- sst$EBS.sst[sst$year %in% 1951:1988]
EBS2 <- sst$EBS.sst[sst$year %in% 1989:2013]

# calculate separate regressions in each era!
# make objects to catch results
GOA.regr1 <- GOA.regr2 <- EBS.regr1 <- EBS.regr2 <- NA

# now loop through each cell
for(i in 1:ncol(win.slp1)){
  #  i <- 1
  mod <- lm(win.slp1[,i] ~ GOA1)
  GOA.regr1[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(win.slp2[,i] ~ GOA2)
  GOA.regr2[i] <- summary(mod)$coef[2,1]
  
  mod <- lm(win.slp1[,i] ~ EBS1)
  EBS.regr1[i] <- summary(mod)$coef[2,1] 
  
  mod <- lm(win.slp2[,i] ~ EBS2)
  EBS.regr2[i] <- summary(mod)$coef[2,1] 
}


# calculate differences for each era
GOA.diff <- GOA.regr2 - GOA.regr1
EBS.diff <- EBS.regr2 - EBS.regr1

diff.lim <- range(GOA.diff, EBS.diff) # limit for plotting


######################################


map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)

plot.points <- mapproject(pink.runs$long, pink.runs$lat,projection=my.proj, orientation=my.orien) 
points(plot.points$x, plot.points$y, pch=21, bg ="red", cex=1)
box()  

mtext("a) Pink", adj=0)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(195,240), ylim=c(49,60),
    fill=FALSE, lforce="e")

map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)

plot.points <- mapproject(sock.runs$long, sock.runs$lat,projection=my.proj, orientation=my.orien) 
points(plot.points$x, plot.points$y, pch=21, bg ="red", cex=1)
box()  
mtext("b) Sockeye", adj=0)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(195,240), ylim=c(49,60),
    fill=FALSE, lforce="e")

map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)

plot.points <- mapproject(chum.runs$long, chum.runs$lat,projection=my.proj, orientation=my.orien) 
points(plot.points$x, plot.points$y, pch=21, bg ="red", cex=1)
box()  
mtext("c) Chum", adj=0)

# Re-shape mean SST data  to a matrix with latitudes in columns, longitudes in rows
# get mean value for each cell
coast.mean <- colMeans(coast.sst)
# turn into matrix for plotting
z <- t(matrix(coast.mean,length(y)))  
# and change to a set of polygons for projection
polys <- matrixPoly(x, y, z)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(193,235.5), ylim=c(35,65),
    fill=FALSE, lforce="e")

COLS <- val2col(z, col = tim.colors(64))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}

# make grid for cells used
plot.lat <- lat[keep]
plot.long <- lon[keep]

for(i in 1:length(plot.lat)){
  xgrid <- c(plot.long[i]-1, plot.long[i]+1, plot.long[i]+1, plot.long[i]-1, plot.long[i]-1)
  ygrid <- c(plot.lat[i]+1, plot.lat[i]+1, plot.lat[i]-1, plot.lat[i]-1, plot.lat[i]+1)
  proj.lines <- mapproject(xgrid, ygrid, projection=my.proj, orientation=my.orien) 
  lines(proj.lines$x, proj.lines$y, lwd=0.5)
}

map('world2Hires', 'usa', fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
map('world2Hires', 'Canada',fill=T, add=T, 
    lwd=0.5, col="lightyellow3", proj=my.proj, parameters = NULL, orient=my.orien)
box()

# add legend strip
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.8
l.l <- 1.2
tc.l <- -0.2
image.plot(z, legend.only=TRUE, horizontal =TRUE,  legend.lab = "Mean annual SST (ºC)", 
           smallplot = c(0.25,0.78,0.2,0.23), 
           legend.cex=0.8,
           legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0))) 
mtext("d) SST data", adj=0)

dev.off()

######################################
# and a series of plots

# define projection details
my.proj <- "orthographic"
my.orien <- c(45,180,0)

png("figs/goa era1 slp-sst.png", 3, 3, units="in", res=300)

par(mar=c(1,0,0,0),  tcl=tc.l, mgp=c(1.5,0.3,0), 
    las=1, cex.axis=0.8, cex.lab=0.8) # , oma=c(0,1,1,0)

ylim <- c(40,90)

new.col <- oce.colorsPalette(64)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(0,350), ylim=c(20,90),
    fill=FALSE, lforce="e")

z <- GOA.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))
polys <- matrixPoly(x, y, z)
COLS <- val2col(z, col = new.col, zlim=c(lim[1], -lim[1]))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan', 'Norway', 'Finland', 'Sweden', 'Denmark',
                     'India', 'Nepal', 'Vietnam', 'Iceland', 'Taiwan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3", proj=my.proj, parameters = NULL, orient=my.orien)

mtext("SLP vs. GOA SST 1950-1988", adj=0.5, cex=0.8, side=1)

dev.off()

####
png("figs/goa era2 slp-sst.png", 3, 3, units="in", res=300)

par(mar=c(1,0,0,0),  tcl=tc.l, mgp=c(1.5,0.3,0), 
    las=1, cex.axis=0.8, cex.lab=0.8) # , oma=c(0,1,1,0)

ylim <- c(40,90)

new.col <- oce.colorsPalette(64)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(0,350), ylim=c(20,90),
    fill=FALSE, lforce="e")

z <- GOA.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))
polys <- matrixPoly(x, y, z)
COLS <- val2col(z, col = new.col, zlim=c(lim[1], -lim[1]))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan', 'Norway', 'Finland', 'Sweden', 'Denmark',
                     'India', 'Nepal', 'Vietnam', 'Iceland', 'Taiwan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3", proj=my.proj, parameters = NULL, orient=my.orien)

mtext("SLP vs. GOA SST 1989-2013", adj=0.5, cex=0.8, side=1)

dev.off()

###########

png("figs/ebs era1 slp-sst.png", 3, 3, units="in", res=300)

par(mar=c(1,0,0,0),  tcl=tc.l, mgp=c(1.5,0.3,0), 
    las=1, cex.axis=0.8, cex.lab=0.8) # , oma=c(0,1,1,0)

ylim <- c(40,90)

new.col <- oce.colorsPalette(64)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(0,350), ylim=c(20,90),
    fill=FALSE, lforce="e")

z <- EBS.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))
polys <- matrixPoly(x, y, z)
COLS <- val2col(z, col = new.col, zlim=c(lim[1], -lim[1]))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan', 'Norway', 'Finland', 'Sweden', 'Denmark',
                     'India', 'Nepal', 'Vietnam', 'Iceland', 'Taiwan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3", proj=my.proj, parameters = NULL, orient=my.orien)

mtext("SLP vs. EBS SST 1950-1988", adj=0.5, cex=0.8, side=1)

dev.off()

#######
png("figs/ebs era2 slp-sst.png", 3, 3, units="in", res=300)

par(mar=c(1,0,0,0),  tcl=tc.l, mgp=c(1.5,0.3,0), 
    las=1, cex.axis=0.8, cex.lab=0.8) # , oma=c(0,1,1,0)

ylim <- c(40,90)

new.col <- oce.colorsPalette(64)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(0,350), ylim=c(20,90),
    fill=FALSE, lforce="e")

z <- EBS.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))
polys <- matrixPoly(x, y, z)
COLS <- val2col(z, col = new.col, zlim=c(lim[1], -lim[1]))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan', 'Norway', 'Finland', 'Sweden', 'Denmark',
                     'India', 'Nepal', 'Vietnam', 'Iceland', 'Taiwan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3", proj=my.proj, parameters = NULL, orient=my.orien)

mtext("SLP vs. EBS SST 1989-2013", adj=0.5, cex=0.8, side=1)

dev.off()

##
# GOA diff
png("figs/GOA era diff slp-sst.png", 3, 3, units="in", res=300)

par(mar=c(1,0,0,0),  tcl=tc.l, mgp=c(1.5,0.3,0), 
    las=1, cex.axis=0.8, cex.lab=0.8) # , oma=c(0,1,1,0)

ylim <- c(40,90)

new.col <- oce.colorsPalette(64)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(0,350), ylim=c(20,90),
    fill=FALSE, lforce="e")

z <- GOA.diff  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))
polys <- matrixPoly(x, y, z)
COLS <- val2col(z, col = new.col, zlim=c(-diff.lim[2], diff.lim[2]))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}


map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan', 'Norway', 'Finland', 'Sweden', 'Denmark',
                     'India', 'Nepal', 'Vietnam', 'Iceland', 'Taiwan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3", proj=my.proj, parameters = NULL, orient=my.orien)

mtext("GOA: (1989-2013)-(1951-1988)", adj=0.5, cex=0.8, side=1)

dev.off()
##
# EBS diff
png("figs/EBS era diff slp-sst.png", 3, 3, units="in", res=300)

par(mar=c(1,0,0,0),  tcl=tc.l, mgp=c(1.5,0.3,0), 
    las=1, cex.axis=0.8, cex.lab=0.8) # , oma=c(0,1,1,0)

ylim <- c(40,90)

new.col <- oce.colorsPalette(64)

map("world2Hires", 
    proj=my.proj, parameters = NULL, orient=my.orien,
    xlim=c(0,350), ylim=c(20,90),
    fill=FALSE, lforce="e")

z <- EBS.diff  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))
polys <- matrixPoly(x, y, z)
COLS <- val2col(z, col = new.col, zlim=c(-diff.lim[2], diff.lim[2]))
for(i in seq(polys)){
  tmp <- mapproject(polys[[i]],
                    proj=my.proj, parameters = NULL, orient=my.orien)
  polygon(tmp$x, tmp$y, col=COLS[i], border=COLS[i], lwd=0.1)
}

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan', 'Norway', 'Finland', 'Sweden', 'Denmark',
                     'India', 'Nepal', 'Vietnam', 'Iceland', 'Taiwan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3", proj=my.proj, parameters = NULL, orient=my.orien)

mtext("EBS: (1989-2013)-(1951-1988)", adj=0.5, cex=0.8, side=1)

dev.off()

#################
# old 'square' version
png("figs/goa-ebs slp-sst.png", 8, 6, units="in", res=300)
# setup the layout
mt.cex <- 1.1
l.mar <- 3
l.cex <- 0.6
l.l <- 0.2
tc.l <- -0.2


par(mar=c(0.5,0.5,1,2),  tcl=tc.l, mgp=c(1.5,0.3,0), 
    las=1, mfrow=c(3,2), cex.axis=0.8, cex.lab=0.8) # , oma=c(0,1,1,0)

ylim <- c(40,90)


new.col <- oce.colorsPalette(64)
# GOA first era
z <- GOA.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=ylim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey", drawlabels=F, lwd=0.5)

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

mtext("SLP vs. GOA SST 1950-1988", adj=0.5, cex=0.8)

# GOA second era
z <- GOA.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=ylim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey", drawlabels=F, lwd=0.5)

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

mtext("SLP vs. GOA SST 1989-2013", adj=0.5, cex=0.8)

# EBS first era
z <- EBS.regr1  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=ylim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey", drawlabels=F, lwd=0.5)

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

mtext("SLP vs. EBS SST 1950-1988", adj=0.5, cex=0.8)

# EBS second era
z <- EBS.regr2  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=new.col, zlim=c(lim[1], -lim[1]), ylim=ylim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey", drawlabels=F, lwd=0.5)

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

mtext("SLP vs. EBS SST 1989-2013", adj=0.5, cex=0.8)

# difference first era
z <- era1.diff  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=new.col, zlim=c(-1.7,1.7), ylim=ylim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey", drawlabels=F, lwd=0.5)

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

mtext("GOA-EBS 1950-1988", adj=0.5, cex=0.8)

# difference second era
z <- era2.diff  # replace elements NOT corresponding to land with loadings!
z <- t(matrix(z, length(y)))  # Convert vector to matrix and transpose for plotting
image.plot(x,y,z, col=new.col, zlim=c(-1.7,1.7), ylim=ylim,
           xlab = "", ylab = "",yaxt="n", xaxt="n", legend.mar=l.mar, legend.line=l.l, axis.args=list(cex.axis=l.cex, tcl=tc.l, mgp=c(3,0.3,0)))

contour(x,y,z, add=T, col="grey", drawlabels=F, lwd=0.5)

map('world2Lores', c('Canada', 'usa', 'USSR', 'Mexico', 'China', 'Mongolia', 'Greenland',
                     'South Korea', 'North Korea', 'Japan', 'Finland', 'Norway', 'Denmark', 'Sweden', 'Ukraine'), 
    fill=T,add=T, lwd=0.5, col="darkgoldenrod3")

mtext("GOA-EBS 1989-2013", adj=0.5, cex=0.8)

dev.off()

