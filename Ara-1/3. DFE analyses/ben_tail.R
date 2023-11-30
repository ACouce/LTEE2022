# .. This script tests fit to exponential and plots a zoomed-in version of the beneficial tail of the distributions

# .. load a useful library (to use fitdistr())
require(MASS)

# .. load beneficial tails using thben as threshold for approximate neutrality
thben<-0.015

R_beneficial<-Rtable[which(Rtable$fitted1>thben),]
R_beneficial<-R_beneficial[!duplicated(R_beneficial$site),] # to impede overlapping regions being counted twice

T_beneficial<-Ttable[which(Ttable$fitted1>thben),]
T_beneficial<-T_beneficial[!duplicated(T_beneficial$site),] # to impede overlapping regions being counted twice

F_beneficial<-Ftable[which(Ftable$fitted1>thben),]
F_beneficial<-F_beneficial[!duplicated(F_beneficial$site),] # to impede overlapping regions being counted twice

# .. shift values to have a zero class for the exponential fit
r <- R_beneficial$fitted1-thben
m <- T_beneficial$fitted1-thben
k <- F_beneficial$fitted1-thben

# .. maximum-likelihood fitting to an exponential distribution
fitR <- fitdistr(r, "exponential") 
fitM <- fitdistr(m, "exponential")
fitK <- fitdistr(k, "exponential")

# .. print p-values
print(ks.test(r, "pexp", fitR$estimate))
print(ks.test(m, "pexp", fitM$estimate))
print(ks.test(k, "pexp", fitK$estimate)) 

# .. initialize graphical output
png("ben_tail.png", width=(2400/2)*0.011, height=(820/2)*0.011, units = 'in', res = 300)
par(mfrow=c(1,3), family="sans", cex=1.35)
mai = c(1, 0.1, 0.1, 0.1)
opar <- par(lwd = 1.5)

# .. retrieve graphical parameters from base R hist() to build a customized, logarithmic histogram
expR<-hist(r, breaks=20, plot=FALSE)
expM<-hist(m, breaks=6, plot=FALSE)
expK<-hist(k, breaks=6, plot=FALSE)
range<-length(expR$density)

# .. plot a customized, logarithmic histogram
Rpoints<-barplot(log(expR$density+1), xlim=c(0,range+1),names.arg=seq(thben,thben+(range-1)*0.005,0.005), space=0, col=rgb(0,0,0,0.25), border='grey', cex.axis=1, cex.lab=2,  xaxt='n', ylab='log (density)', yaxt='n')
lines(Rpoints-0.5, log(expR$density+1), type='s', lwd=2, col='grey24')
lines(Rpoints, log(dexp(expR$breaks, rate = fitR$estimate))[1:length(Rpoints)], col = "black", lwd=2, lty=2)

# .. custom axes
axis(2, lwd = 1.5)
ticks<-1:length(expR$breaks)
axis(1, at = ticks, pos=-0.2, labels = FALSE, cex.axis = 1, tck=-0.03, lwd=1.5)
axis(1, at = c(2,7,12,17,22), pos=-0.2, labels = round((expR$breaks+0.020)[c(2,7,12,17,22)],3), lwd=1.5, cex.axis = 1, tck=-0.08)

# .. plot a customized, logarithmic histogram
Mpoints<-barplot(log(c(expM$density,rep(0,range-length(expM$density)))+1),  xlim=c(0,range+1),names.arg=seq(thben,thben+(range-1)*0.005,0.005), space=0, col=col2K, border=col2K, cex.axis=1, cex.names=1,  xaxt='n', yaxt='n')
lines(Rpoints-0.5, log(c(expM$density,rep(0,range-length(expM$density)))+1), type='s', lwd=2, col='grey24')
lines(Mpoints, log(dexp(expM$breaks, rate = fitM$estimate))[1:length(Mpoints)], col = "black", lwd=2, lty=2)

# .. custom axes
axis(side = 2, lwd = 1.5)
axis(1, at = ticks, pos=-0.2, labels = FALSE, cex.axis = 1, tck=-0.03, lwd=1.5)
axis(1, at = c(1,5,9), pos=-0.2, labels = round((expR$breaks+0.020)[c(1,5,9)],3), lwd=1.5, cex.axis = 1, tck=-0.08)

# .. plot a customized, logarithmic histogram
Kpoints<-barplot(log(c(expK$density,rep(0,range-length(expK$density)))+1),  xlim=c(0,range+1),names.arg=seq(thben,thben+(range-1)*0.005,0.005), space=0, col=col15K, border=col15K, cex.axis=1, cex.names=1,  xaxt='n', yaxt='n')
lines(Rpoints-0.5, log(c(expK$density,rep(0,range-length(expK$density)))+1), type='s', lwd=2, col='grey24')
lines(Kpoints, log(dexp(expK$breaks, rate = fitK$estimate))[1:length(Kpoints)], col = "black", lwd=2, lty=2)

# .. custom axes
axis(side = 2, lwd = 1.5)
axis(1, at = ticks, pos=-0.2, labels = FALSE, cex.axis = 1, tck=-0.03, lwd=1.5)
axis(1, at = c(1,5,9), pos=-0.2, labels = round((expR$breaks+0.020)[c(1,5,9)],3), lwd=1.5, cex.axis = 1, tck=-0.08)

# .. close graphical output
dev.off()
