# .. This script plots the DFEs

# .. define a simple smoothing function
mav <- function(x,n){filter(x,rep(1/n,n), method = c("convolution"), sides=1)}

# .. load count table of fitness classes
R6<-read.table('Rhist.dat')
tK<-read.table('2Khist.dat')
fK<-read.table('15Khist.dat')

# .. initialize graphical output
png(file="DFE.png", width = 15.63, height = 8.54, units = 'in', res = 300)
par(mfrow=c(1,2), family="sans", cex=2, mgp=c(1.5,0.5,0))

# .. the plot
plot(R6[,1], R6[,2], log='y',xlim=c(-0.3,0.1),ylim=c(1,9000), las=1, type='l', lwd=4, pch=1, cex=1, xlab = "selection coefficient (s)",  ylab = "number of insertions", cex.lab=0.75,  cex.axis=0.67, yaxt="n", col= 'white')
lines(R6[,1]-1*step/2,mav(R6[,3],2),type='l',lwd=10, col='black')
lines(tK[,1]-1*step/2,mav(tK[,3],2),type='l',lwd=10, col=col2K)
lines(fK[,1]-1*step/2,mav(fK[,3],2),type='l',lwd=10, col=col15K)
legend("topleft", legend=c("Anc", "2K", "15K"),  col=c("black", col2K, col15K), lty = 1, lwd=8, cex=0.75, box.lty=0)
abline(v=0,lty=2,lwd=1)

# .. custom y-axis
ylab=seq(0,4,by=1)
yat=10^ylab
axis(2, at=yat, labels=yat,las=1, tcl=-0.5, las=3, lwd=1.5, cex.axis=0.67, cex.lab=0.75)

ylab=seq(2,9,by=1)
for (i in seq(0,3,1)) {
	yat=ylab*10^i
	axis(2, at=yat, labels=NA,las=1, tcl=-0.25, lwd=1.5, cex.axis=0.67, cex.lab=0.75)
}

# .. homogenize border lines
box(lwd=1.5)

# .. close graphical output
dev.off()
