# .. This script resamples datasets to assess the role of noise in statistical analyses of the shape of DFE.

# .. make condensed version of deleterious tails with only fitness effects and standard errors
Rtail<-Rtable[Rtable$fitted1>-0.5 & Rtable$fitted1<=-0.015, c('fitted1', 'sterr1')]
Ttail<-Ttable[Ttable$fitted1>-0.5 & Ttable$fitted1<=-0.015,c('fitted1', 'sterr1')]
Ftail<-Ftable[Ftable$fitted1>-0.5 & Ftable$fitted1<=-0.015,c('fitted1', 'sterr1')]

# .. matrix to store p-values from KS tests
delsam<-matrix(NA,200, 3)

# .. conduct resampling with replacement using the logarithm of the regression standard errors as weights
for (b in 1:200) {

	Rr<-sample(Rtail$fitted1, nrow(Rtail), replace = TRUE, prob= log(1/Rtail$sterr1))
	Rr2<-sample(Rtail$fitted1, nrow(Rtail), replace = TRUE, prob= log(1/Rtail$sterr1))
	Tr<-sample(Ttail$fitted1, nrow(Ttail), replace = TRUE, prob= log(1/Ttail$sterr1))
	Fr<-sample(Ftail$fitted1, nrow(Ftail), replace = TRUE, prob= log(1/Ftail$sterr1))

	KS1<-ks.test(Rr, Rr2,alternative = "two.sided")	# Anc resampled delterious tail compared to another realization of Anc resampled delterious tail
	KS2<-ks.test(Rr, Tr,alternative = "two.sided")	# Anc vs 2K
	KS3<-ks.test(Rr, Fr,alternative = "two.sided")	# Anc vs 15K
	
	delsam[b,1:3]<-c(KS1$p.value,KS2$p.value, KS3$p.value)
}

# .. print average p-values 
print(colSums(delsam)/200)

# .. initialize graphical output
png(file="resampling.png", width=900, height=1050, res=64)
par(mfrow=c(2,2), family="sans", cex=2)
par(lwd=1.5)

# .. set a lower limit for visualization purposes
delsam[delsam<1e-8]<-1e-8

# .. plot p-values as scatter plot for Anc, 2K and 15K data sets
plot(delsam[,1], log='y', ylim=c(1e-8,1), main=signif(median(delsam[,1]), 2), xlab='resample' , ylab='p-value', cex=1.5, lwd=0.2, pch=21, bg=rgb(0,0,0,0.5), yaxt="n")
points(delsam[,2], cex=1.5, lwd=0.2, pch=21, bg=rgb(0.09,0.45,0.95,0.5))
points(delsam[,3], cex=1.5, lwd=0.2, pch=21, bg=rgb(0.959,0.25,0.22,0.5))
abline(h=0.01, lwd=0.75, lty=3, col='black')

# .. custom y-axis
ylab=seq(0,-8,by=-1)
yat=10^ylab
axis(2, at=yat, labels=yat,las=1, tcl=-0.75, las=3, lwd=1.5, cex.axis=1, cex.lab=1)

ylab=seq(2,9,by=1)
for (i in seq(0,-8,-1)) {
yat=ylab*10^i
axis(2, at=yat, labels=NA,las=1, tcl=-0.25, lwd=0.5, cex.axis=1, cex.lab=1)
}

# .. plot p-values as boxplot 
boxplot(delsam, log='y', ylim=c(1e-8,1), col=c(rgb(0,0,0,0.7), rgb(0.09,0.45,0.95,0.75), rgb(0.959,0.25,0.22,0.75)), yaxt="n", cex=1)
abline(h=0.01, lwd=0.75, lty=3, col='black')

# .. custom y-axis
ylab=seq(0,-8,by=-1)
yat=10^ylab
axis(2, at=yat, labels=yat,las=1, tcl=-0.75, las=3, lwd=1.5, cex.axis=1, cex.lab=1)

ylab=seq(2,9,by=1)
for (i in seq(0,-8,-1)) {
yat=ylab*10^i
axis(2, at=yat, labels=NA,las=1, tcl=-0.25, lwd=0.5, cex.axis=1, cex.lab=1)
}


# .. now, let's make condensed versions of beneficial tails with only fitness effects and standard errors
Rtail<-Rtable[which(Rtable$fitted1>0.015), c('fitted1', 'sterr1')]
Ttail<-Ttable[which(Ttable$fitted1>0.015), c('fitted1', 'sterr1')]
Ftail<-Ftable[which(Ftable$fitted1>0.015), c('fitted1', 'sterr1')]

# .. load a useful library (to use fitdist(), which implements MLE and GoF analyses for all the desired distributions) 
#library(fitdistrplus)	# uncomment in case it was not already loaded

# .. matrix to store p-values from KS tests
bensam<-matrix(NA,200, 3)

# .. conduct resampling with replacement using the logarithm of the regression standard errors as weights
for (b in 1:200) {

	r<-sample(Rtail$fitted1, nrow(Rtail), replace = TRUE)-0.015	# shift values to have a zero class for the exponential fit
	m<-sample(Ttail$fitted1, nrow(Ttail), replace = TRUE)-0.015
	k<-sample(Ftail$fitted1, nrow(Ftail), replace = TRUE)-0.015
	
	expR <- fitdist(r, "exp")
	expM <- fitdist(m, "exp") 
	expK <- fitdist(k, "exp")
	
	bensam[b,1:3]<-c(ks.test(r, "pexp", expR$estimate)$p.value, ks.test(m, "pexp", expM$estimate)$p.value, ks.test(k, "pexp", expK$estimate)$p.value)
}

# .. print average p-values 
print(c(signif(mean(bensam[,1]), 2), round(mean(bensam[,2]),2),round(mean(bensam[,3]),2)))

# .. set a lower limit for visualization purposes
bensam[bensam<1e-8]<-1e-8

# .. plot p-values as scatter plot for Anc, 2K and 15K data sets
plot(bensam[,1], log='y', ylim=c(1e-8,1), main=signif(median(bensam[,1]), 2), xlab='resample' , ylab='p-value', cex=1.5, lwd=0.2, pch=21, bg=rgb(0,0,0,0.5), yaxt="n")
points(bensam[,2], cex=1.5, lwd=0.2, pch=21, bg=rgb(0.09,0.45,0.95,0.5))
points(bensam[,3], cex=1.5, lwd=0.2, pch=21, bg=rgb(0.959,0.25,0.22,0.5))
abline(h=0.01, lwd=0.75, lty=3, col='black')

# .. custom y-axis
ylab=seq(0,-10,by=-1)
yat=10^ylab
axis(2, at=yat, labels=NA,las=1, tcl=-0.75, las=3, lwd=1.5, cex.axis=1, cex.lab=1)
axis(2, at=10^seq(0,-10,by=-2), labels=expression(10^0, 10^-2, 10^-4, 10^-6, 10^-8, 10^-10),las=1, tcl=-0.75, las=3, lwd=1.5, cex.axis=1, cex.lab=1)

ylab=seq(2,9,by=1)
for (i in seq(0,-10,-1)) {
yat=ylab*10^i
axis(2, at=yat, labels=NA,las=1, tcl=-0.25, lwd=0.5, cex.axis=1, cex.lab=1)
}

# .. plot p-values as boxplot 
boxplot(bensam[,1:3], log='y', ylim=c(1e-8,1), col=c(rgb(0,0,0,0.7), rgb(0.09,0.45,0.95,0.75), rgb(0.959,0.25,0.22,0.75)), yaxt="n", cex=1)
abline(h=0.01, lwd=0.75, lty=3, col='black')

# .. custom y-axis
ylab=seq(0,-10,by=-1)
yat=10^ylab
axis(2, at=yat, labels=yat,las=1, tcl=-0.75, las=3, lwd=1.5, cex.axis=1, cex.lab=1)

ylab=seq(2,9,by=1)
for (i in seq(0,-10,-1)) {
yat=ylab*10^i
axis(2, at=yat, labels=NA,las=1, tcl=-0.25, lwd=0.5, cex.axis=1, cex.lab=1)
}

# .. close graphical output
dev.off()
