# ..  This script plots the the predictive capacity of DFEs as a function of time in the LTEE, both for actual and simulated data. It also plots prevalence of loss-of-function mutations (if last bit is uncommented)

# .. load useful library (to use function std.error)
library(plotrix)

# .. define a simple smoothing function
mav <- function(x,n){filter(x,rep(1/n,n), method = c("convolution"), sides=2)}

# .. define custom colors
col2K<-'#cd6e6c'
col15K<-'#6c91cd'

# .. define generations at which the meta-genomic data was sampled (in case needed)
gens<-c(seq(1000,50000,500))

# .. load representativity results with real data
Rf<-read.table(file = 'Rrep.dat', sep = "\t", header = TRUE)
Tf<-read.table(file = 'Trep.dat', sep = "\t", header = TRUE)
Ff<-read.table(file = 'Frep.dat', sep = "\t", header = TRUE)

# .. load representativity results with simulated data (for neutral expectation)
Rsims<-read.table(file = 'Rrep_sims.dat', sep = "\t", header = TRUE)
Tsims<-read.table(file = 'Trep_sims.dat', sep = "\t", header = TRUE)
Fsims<-read.table(file = 'Frep_sims.dat', sep = "\t", header = TRUE)

# .. initialize graphical output
png(file="repre.png", width = (15.63*0.9)*(2/2.32), height = (8.54)/2.22, units = 'in', res = 300)
par(mfrow=c(1,3), family="sans", cex=1) #, mgp=c(1.5,0.5,0))

# .. create an empty plot with desired details and parameters
plot(gens, Rf[1:99,], xlim=c(1000,50000), ylim=c(0.015,0.6), type='o', col='white', cex.axis=1.35, cex.lab=1.3, xlab='generations', ylab='fraction of drivers')#, log='y')#,log='y')#, log='y')#, log='xy')

# .. add reference lines for 2K and 15K
abline(v=c(2000,15000), lty=2, col='darkgrey')

# .. smooth and plot ancestor data
Rf[2:99,]<-mav(Rf[2:99,],2)
lines(gens,Rf[1:99,1],type='l',lwd=1.5, col='black')
lines(gens,Rf[1:99,1],type='p',lwd=1.5, cex = 1.4, pch=21, bg='#BFBFBF',col='black')

# .. create and add neutral expectation with 90% confidence intervals 
Ravg<-apply(Rsims, 1, mean)
Rsd<-apply(Rsims, 1, std.error)
polygon(x = c(gens,rev(gens)), y = c(Ravg[1:99]+Rsd[1:99]*1.64, rev(Ravg[1:99]-Rsd[1:99]*1.64)), col= adjustcolor("black", alpha.f = 0.25), border=NA)

# .. smooth and plot 2K data
Tf[2:99,]<-mav(Tf[2:99,],2)
Tf[Tf==0]<-NA
lines(gens,Tf[1:99,1],type='l',lwd=1.5, col=col2K)
lines(gens,Tf[1:99,1],type='p',lwd=1.5, cex = 1.4, pch=21, bg='#f1d1e1',col=col2K)

# .. create and add neutral expectation with 90% confidence intervals 
Tavg<-apply(Tsims, 1, mean)
Tsd<-apply(Tsims, 1, std.error)
yvalues<-c(Tavg[1:99]+Tsd[1:99]*1.64, rev(Tavg[1:99]-Tsd[1:99]*1.64))
yvalues[is.na(yvalues)]<-0
polygon(x = c(gens,rev(gens)), y = yvalues, col= adjustcolor(col2K, alpha.f = 0.25), border=NA)

# .. smooth and plot 15K data
Ff[2:99,]<-mav(Ff[2:99,],2)
Ff[Ff==0]<-NA
lines(gens,Ff[1:99,1],type='l',lwd=1.5, col=col15K)
lines(gens,Ff[1:99,1],type='p',lwd=1.5, cex = 1.4, pch=21, bg='#d1e1f1',col=col15K)

# .. create and add neutral expectation with 90% confidence intervals 
Favg<-apply(Fsims, 1, mean)
Fsd<-apply(Fsims, 1, std.error)
yvalues<-c(Favg[1:99]+Fsd[1:99]*1.64, rev(Favg[1:99]-Fsd[1:99]*1.64))
yvalues[is.na(yvalues)]<-0
polygon(x = c(gens,rev(gens)), y = yvalues, col= adjustcolor(col15K, alpha.f = 0.25), border=NA)

# .. close graphical output
dev.off()

# uncomment all below to plot prevalence of loss-of-function mutations
#LOF<-read.table(file = 'LOF.dat', sep = "\t", header = TRUE)
#png(file="LOF.png",  width=1250/2.2, height=810/1.5, res=80)
#par(family="sans", cex=2, lwd=1.5)
#plot(gens, LOF[1:99,], xlim=c(1000,50000), ylim=c(0.02,0.8), type='o', col='white', cex.axis=1, cex.lab=1.3, xlab='generations', ylab='fraction of drivers')#, log='y')#,log='y')#, log='y')#, log='xy')
#lines(gens,mav(LOF[1:99,],2),type='l',lwd=1.3, col='black')
#lines(gens,mav(LOF[1:99,],2),type='p',lwd=1.3, cex = 1.4, pch=21, bg='#BFBFBF',col='black')
#lines(gens,c(Rf[1], mav(Rf[2:99],2)),type='l',lwd=1.5, col='#737373', lty=3)
#legend("topright", legend=c("LOF", 'average'),  col=c('#737373', '#737373'), lty = c(1,3), lwd= c(8,3), cex=0.75, box.lty=0)
#dev.off()
