
source('~/Documents/GitTest/mixture_distributions_for_forecasting/simulate_functions.R')


set.seed(168)
# lnorm_params <- getLnormParams(par=c(3,1))
lnorm_params <- c(1,.4)

x <- seq(0,8,length.out = 1001)
# y <- dlnorm(x,lnorm_params[1],lnorm_params[2])

y <- EnvStats::dlnormTrunc(x,meanlog = lnorm_params[1],
                           sdlog = lnorm_params[2],min = x[1],
                           max=x[length(x)])

yp <- EnvStats::plnormTrunc(x,meanlog = lnorm_params[1],
                                 sdlog = lnorm_params[2],min = x[1],
                                 max=x[length(x)])

par(mfrow=c(2,2))
plot(y~x,col='blue',type='l',main='(a)',
     ylab = 'p(x)')

wid <- .2
lnbins <- seq(0,8,by=wid)
disc_lnorm <- discretizeDist(bins=lnbins,par=lnorm_params,
                             family='trlnorm',nu=NA)

lnorm_samp <- EnvStats::rlnormTrunc(1000,meanlog = lnorm_params[1],
                                    sdlog = lnorm_params[2],min = x[1],
                                    max=x[length(x)])
hist(lnorm_samp,main='(b)',xlab='x',border='purple',axes=TRUE)
box()

secdf <- ecdf(lnorm_samp)


outl <- disc_lnorm$bins
plot((disc_lnorm$probs*5)~disc_lnorm$bins,type='s',col='red',xlab='x',lwd=1.5,
     ylab='Prob',main='(c)')
lines((disc_lnorm$probs*5)~outl,type='h',col='black')
lines((disc_lnorm$probs*5)~outl,type='h',col='red')



quant_23 <- c(0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 
              0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700,
              0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990)

quant_7 <- c(0.025, 0.100, 0.250, 0.500, 0.750, 0.900, 0.975)
quants <- quant_23

qvals <- EnvStats::qlnormTrunc(quants,meanlog = lnorm_params[1],
                               sdlog = lnorm_params[2],min = x[1],
                               max=x[length(x)])
dqvals <- EnvStats::dlnormTrunc(qvals,meanlog = lnorm_params[1],
                                sdlog = lnorm_params[2],min = x[1],
                                max=x[length(x)])
tspot <- qvals-.5
plot(y~x,col='blue',type='l',main='(d)',
     ylab = 'p(x)')
lines(dqvals~qvals,type='h',col='green')
points(dqvals~qvals,col='forestgreen',pch=19,cex=.8)
# text(dqvals~tspot,labels=quants)


par(mfrow=c(2,2))
plot(yp~x,col='blue',type='l',main='(a)',
     ylab = 'F(x)')


secdf
plot(secdf,col='purple',main='(b)',col.01line = NULL)


disc <- cbind(disc_lnorm,cmf=cumsum(disc_lnorm$probs))

plot(disc$cmf~disc$bins,type='s',col='red',main='(c)',ylab='P(X ≤ x)',
     xlab='x')

plot(quants~qvals,pch=19,cex=.6,col='forestgreen',main='(d)',
     xlab='x',ylab='Quantile')







plot(x, freq = equidist, density = NULL, angle = 45,
     col = NULL, border = par("fg"), lty = NULL,
     main = paste("Histogram of",
                  paste(x$xname, collapse = "\n")),
     sub = NULL, xlab = x$xname, ylab,
     xlim = range(x$breaks), ylim = NULL,
     axes = TRUE, labels = FALSE, add = FALSE,
     ann = TRUE)




















outl <- disc_lnorm$bins
plot(disc_norm$probs~disc_norm$bins,type='h',col='red',xlab='x',lwd=4,
     ylab='Probabilities',main='Discretized and Standardized Normal')
lines(disc_norm$probs~outl,type='h')
# points(disc_norm$probs~disc_norm$bins,pch=20,col='pink')
library(tidyverse)

quant_23 <- c(0.010, 0.025, 0.050, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 
              0.400, 0.450, 0.500, 0.550, 0.600, 0.650, 0.700,
              0.750, 0.800, 0.850, 0.900, 0.950, 0.975, 0.990)

quant_7 <- c(0.025, 0.100, 0.250, 0.500, 0.750, 0.900, 0.975)
quants <- quant_23


qvals <- qnorm(quants,0,1)
dqvals <- dnorm(qvals,0,1)
tspot <- qvals-.5
plot(y~x,col='blue',type='l',main='Standard Normal Quantiles',
     ylab = 'p(x)')
lines(dqvals~qvals,type='h',col='green')
points(dqvals~qvals,col='forestgreen',pch=19,cex=.6)
text(dqvals~tspot,labels=quants)
